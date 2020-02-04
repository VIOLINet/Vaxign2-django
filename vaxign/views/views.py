from django.shortcuts import render, redirect

# Create your views here.

import re
import urllib
import string
import random
import collections
import phpserialize
import xml.etree.ElementTree as ET

from django.conf import settings
from django.db.models import Q
from django.db.models import Sum
from vaxign.models import TVaxignQuery
from vaxign.models import TUserVaxignPrivilege
from vaxign.models import TVaxignProject
from vaxign.models import TVaxignAnalysisResults
from vaxign.models import TVaxignAlleleGroup

from vaxign.forms.results import ResultsForm
from vaxign.forms.runs import RunsForm

import logging
logger = logging.getLogger('console')

def index(request, name):
    if 'c_user_name' not in request.session.keys():
        try:
            php_session_str = open(str.format('{}/sess_{}', settings.SESSION_FILE_PATH, request.COOKIES.get('PHPSESSID'))).read().encode()
            php_session = phpserialize.loads(php_session_str)
            request.session['c_user_name'] = php_session['c_user_name'.encode()].decode()
            request.session['is_admin'] = php_session['is_admin'.encode()]
        except:
            pass
    else:
        try:
            php_session_str = open(str.format('{}/sess_{}', settings.SESSION_FILE_PATH, request.COOKIES.get('PHPSESSID'))).read().encode()
            php_session = phpserialize.loads(php_session_str)
            request.session['c_user_name'] = php_session['c_user_name'.encode()].decode()
            request.session['is_admin'] = php_session['is_admin'.encode()]
        except:
            if 'c_user_name' in request.session.keys():
                del request.session['c_user_name']
            if 'is_admin' in request.session.keys():
                del request.session['is_admin']
    
    if name == 'index' or name not in ['dynamic', 'precompute']:
        name = 'dynamic'
    context = {'vaxign_option': name}
    
    # Form 1: Dynamic Analysis
    if 'Submit1' in request.POST:
        tmpData = {}
        for key in request.POST:
            if key == 'basic_analysis':
                tmpData[key] = request.POST.getlist(key)
            else:
                tmpData[key] = request.POST[key]
        if 'sequence_file' in request.FILES:
            file = request.FILES['sequence_file']
            fileContent = ''
            for chunk in file.chunks():
                fileContent += chunk.decode()
            tmpData['sequence'] = fileContent.strip()
            tmpData['sequence_type'] = 'protein_fasta'
            
        # Process sequence data
        # Sequence input type: NCBI Protein GI or Refseq
        if tmpData['sequence_type'] == 'protein_gi':
            logger.debug("Selected NCBI Protein GI. Retrieving sequence from NCBI...")
            try:
                url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=text&rettype=fasta&id=' + re.sub('[,\s]+', ',', tmpData['sequence'])
                with urllib.request.urlopen(url) as file:
                    tmpData['sequence'] = file.read().decode()
            except:
                tmpData['sequence'] = None
        # Sequence input type: NCBI Protein Refseq
        if tmpData['sequence_type'] == 'protein_refseq':
            logger.debug( "Selected NCBI Protein Refseq. Retrieving sequence from NCBI...")
            try:
                url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=text&rettype=fasta&id=' + re.sub('[,\s]+', ',', tmpData['sequence'])
                with urllib.request.urlopen(url) as file:
                    tmpData['sequence'] = file.read().decode()
            except:
                tmpData['sequence'] = None
        # Sequence input type: FASTA URL
        if tmpData['sequence_type'] == 'protein_fasta_url':
            logger.debug("Selected custom protein FASTA file link(s). Retrieving sequence from the source(s)...")
            tmp = ''
            for url in re.split('[\r\n]+', tmpData['sequence']):
                with urllib.request.urlopen(url.strip()) as file:
                    tmp += file.read().decode() + '\n'
            tmpData['sequence'] = tmp
        # Sequence input type: NCBI Gene ID
        if tmpData['sequence_type'] == 'gene_id':
            logger.debug("Selected NCBI Gene ID. Retrieving sequence from NCBI...")
            url1 = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=protein&term=srcdb_refseq[PROP]&id=' + re.sub('[,\s]+', '', tmpData['sequence'])
            with urllib.request.urlopen(url1) as file1:
                geneIDs = list(set(re.findall('<Link>\s+<Id>(\d+)<\/Id>\s+<\/Link>', file1.read().decode())))
            url2 = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=text&rettype=fasta&id=' + ','.join(geneIDs)
            with urllib.request.urlopen(url2) as file2:
                tmpData['sequence'] = file2.read().decode()
        # Sequence input type: UniprotKB Protein ID
        if tmpData['sequence_type'] == 'protein_uniprotkb':
            logger.debug("Selected UniprotKB Protein ID. Retrieving sequence from Uniprot...")
            url = str.format('https://www.uniprot.org/uniprot/?query=id:{}&format=fasta', re.sub('[,\s]+', '+OR+id:', tmpData['sequence']))
            with urllib.request.urlopen(url) as file:
                tmpData['sequence'] = file.read().decode()
        # Sequence input type: NCBI Bioproject ID
        if tmpData['sequence_type'] == 'bioproject_id':
            logger.debug( "Selected NCBI Bioproject ID. Retrieving sequence from NCBI...")
            if tmpData['sequence'].strip().startswith( 'PRJNA' ):
                bioprojectID = tmpData['sequence'].strip().replace('PRJNA', '')
            else:
                bioprojectID = tmpData['sequence'].strip()
            url1 = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=bioproject&db=protein&id='+bioprojectID
            proteinIDs = []
            with urllib.request.urlopen(url1) as file1:
                tree1 = ET.fromstring(file1.read().decode())
                for proteinID in tree1.find('LinkSet').find('LinkSetDb').iter('Id'):
                    proteinIDs.append(proteinID.text)
            url2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            data2 = urllib.parse.urlencode({
                'db': 'protein',
                'rettype': 'fasta',
                'retmode': 'text',
                'retmax': '10000',
                'id': ','.join(proteinIDs),
            }).encode()
            request2 = urllib.request.Request(url=url2, data=data2)
            with urllib.request.urlopen(request2) as file2:
                tmpData['sequence'] = file2.read().decode()
            
            url3 = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=bioproject&id='+bioprojectID
            with urllib.request.urlopen(url3) as file3:
                tree3 = ET.fromstring(file3.read().decode())
                organism = tree3.find('DocumentSummarySet').find('DocumentSummary').find('Organism_Name').text
                strain = tree3.find('DocumentSummarySet').find('DocumentSummary').find('Organism_Strain').text
                if tmpData['group_name'] == '':
                    tmpData['group_name'] = organism
                if tmpData['genome_name'] == '':
                    tmpData['genome_name'] = str.format('{} {}', organism, strain)
        # Sequence input type: Uniprot Proteome ID
        if tmpData['sequence_type'] == 'protein_uniprot_proteome':
            logger.debug("Selected Uniprot Proteome ID. Retrieving sequence from Uniprot...")
            url = str.format('https://www.uniprot.org/uniprot/?query=proteome:{}&format=fasta', tmpData['sequence'].strip())
            with urllib.request.urlopen(url) as file:
                tmpData['sequence'] = file.read().decode()
            
        form1 = RunsForm(tmpData or None)
        context['vaxign_option'] = 'dynamic'
    else:
        form1 = RunsForm()
    if request.method == 'POST':
        if form1.is_valid():
            formData = form1.cleaned_data
            queryID = ''.join(random.choice('ABCDEFGHJKMNPQRSTUVWXYZ23456789') for i in range(10))
            formData['c_projectfk'] = 0
            formData['group_name'] = ''
            formData['genome_name'] = ''
            request.session['runs_form_data'] = formData
            return redirect('/vaxign2/run/'+queryID)
    context['form1'] = form1
    
    # Form 2: Precompute Query
    if 'Submit2' in request.POST:
        form2 = ResultsForm(request.POST or None)
        context['vaxign_option'] = 'precompute'
    else:
        form2 = ResultsForm()
    if request.method == 'POST':
        if form2.is_valid():
            if 'query_seq_ids_file' in request.FILES:
                file = request.FILES['query_seq_ids_file']
                fileContent = ''
                for chunk in file.chunks():
                    fileContent += chunk.decode()
                form2.cleaned_data['query_seq_ids'] = fileContent.strip()
            
            formData = form2.cleaned_data
            queryID = formData['query_id']
            request.session['results_form_data'] = formData
            return redirect('/vaxign2/query/'+queryID+'/results')
        
    context['form2'] = form2
    
    if 'c_user_name' in request.session:
        privilege = TUserVaxignPrivilege.objects.filter(
            Q(c_userfk=request.session['c_user_name'])
        )
        projects = TVaxignProject.objects.filter(
            Q(c_vaxign_projectid__in=privilege.values_list('c_projectfk', flat=True)) & Q(c_list_title=1)
        )
        queries = TVaxignQuery.objects.filter(
            Q(c_is_public=1) & ~Q(c_species_short='') |
            Q(c_projectfk__in=projects.values_list('c_vaxign_projectid', flat=True) ) & ~Q(c_species_short='')
            ).order_by('c_species_short', 'c_species_name')
    else:
        queries = TVaxignQuery.objects.filter(
            Q(c_is_public=1) & ~Q(c_species_short=''),
        ).order_by('c_species_short', 'c_species_name')
    
    genomeGroups = []
    genomeGroupsCount = []
    for row in queries:
        if row.c_species_short not in genomeGroups:
            genomeGroups.append(row.c_species_short)
            genomeGroupsCount.append(0)
        index = genomeGroups.index(row.c_species_short)
        genomeGroupsCount[index] += 1
    
    context['genome_groups_count'] = zip(genomeGroups, genomeGroupsCount)
    
    return render(request, 'home.html', context)

def aucs(request):
    
    return render(request, 'docs/aucs.html')

def disclaimer(request):
    
    return render(request, 'docs/disclaimer.html')

def docs(request):
    
    return render(request, 'docs/docs.html')

def faqs(request):
    
    return render(request, 'docs/faqs.html')

def tutorial(request):
    
    return render(request, 'docs/tutorial.html')

def updates(request):
    
    return render(request, 'docs/updates.html')

def stats(request):
    
    context = {}
    
    queries = TVaxignQuery.objects.filter(
        Q(c_is_public=1) & ~Q(c_species_short='')
    ).order_by('c_species_short', 'c_species_name')
    
    groupGenCount = {}
    groupSeqCount = {}
    genomeSeqCount = {}
    totalCount = 0
    for query in queries:
        count = TVaxignAnalysisResults.objects.filter(
            Q(c_query_id=query.c_query_id)
        ).count()
        genomeSeqCount[query.c_species_name] = count 
        if query.c_species_short not in groupSeqCount.keys():
            groupSeqCount[query.c_species_short] = 0
            groupGenCount[query.c_species_short] = 0
        groupSeqCount[query.c_species_short] += count
        groupGenCount[query.c_species_short] += 1
        totalCount += count
    
    current = ''
    groupIndex = 1
    genomeIndex = 1
    results = collections.OrderedDict()
    for row in queries.order_by('c_species_short', 'c_species_name'):
        if current != row.c_species_short:
            results[row.c_species_short] = [row.c_species_short, groupIndex, row.c_species_short, groupGenCount[row.c_species_short], groupSeqCount[row.c_species_short]] 
            current = row.c_species_short
            groupIndex += 1
            genomeIndex = 1
        if genomeIndex % 2 == 0:
            results[row.c_query_id] = [row.c_species_short, genomeIndex, row.c_species_name, row.c_ref_acc, genomeSeqCount[row.c_species_name], '#F8FAFA']
        else:
            results[row.c_query_id] = [row.c_species_short, genomeIndex, row.c_species_name, row.c_ref_acc, genomeSeqCount[row.c_species_name], '#EEEEEE']
        genomeIndex += 1
    context['results'] = results
    context['group_count'] = len(groupGenCount.keys())
    context['genome_count'] = len(genomeSeqCount.keys())
    context['total_count'] = totalCount
    
    return render(request, 'statistics.html', context)

def login(request):
    
    if 'c_user_name' not in request.session.keys() and 'is_admin' not in request.session.keys():
        try:
            php_session_str = open(str.format('{}/sess_{}', settings.SESSION_FILE_PATH, request.COOKIES.get('PHPSESSID'))).read().encode()
            php_session = phpserialize.loads(php_session_str)
            request.session['c_user_name'] = php_session['c_user_name'.encode()].decode()
            request.session['is_admin'] = php_session['is_admin'.encode()]
        except:
            pass
    
    return redirect('/vaxign2')
    
def logout(request):
    
    if 'c_user_name' in request.session.keys():
        del request.session['c_user_name']
    if 'is_admin' in request.session.keys():
        del request.session['is_admin']
    
    return redirect('/vaxign2')

def vaxignml(request):
    
    context = {}
    
    if 'Submit' in request.POST:
        tmpData = {}
        for key in request.POST:
            tmpData[key] = request.POST[key]
        tmpData['vaxignml_analysis'] = 'y'
        tmpData['vaxitop_analysis'] = 'n'
        if 'sequence_file' in request.FILES:
            file = request.FILES['sequence_file']
            fileContent = ''
            for chunk in file.chunks():
                fileContent += chunk.decode()
            tmpData['sequence'] = fileContent.strip()
            tmpData['sequence_type'] = 'protein_fasta'
            
        # Process sequence data
        # Sequence input type: NCBI Protein GI or Refseq
        if tmpData['sequence_type'] == 'protein_gi':
            logger.debug( "Selected NCBI Protein GI. Retrieving sequence from NCBI...")
            try:
                url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=text&rettype=fasta&id=' + re.sub('[,\s]+', ',', tmpData['sequence'])
                with urllib.request.urlopen(url) as file:
                    tmpData['sequence'] = file.read().decode()
            except:
                tmpData['sequence'] = None
        # Sequence input type: NCBI Protein Refseq
        if tmpData['sequence_type'] == 'protein_refseq':
            logger.debug( "Selected NCBI Protein Refseq. Retrieving sequence from NCBI...")
            try:
                url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=text&rettype=fasta&id=' + re.sub('[,\s]+', ',', tmpData['sequence'])
                with urllib.request.urlopen(url) as file:
                    tmpData['sequence'] = file.read().decode()
            except:
                tmpData['sequence'] = None
        # Sequence input type: FASTA URL
        if tmpData['sequence_type'] == 'protein_fasta_url':
            logger.debug( "Selected custom protein FASTA file link(s). Retrieving sequence from the source(s)...")
            tmp = ''
            for url in re.split('[\r\n]+', tmpData['sequence']):
                with urllib.request.urlopen(url.strip()) as file:
                    tmp += file.read().decode() + '\n'
            tmpData['sequence'] = tmp
        # Sequence input type: NCBI Gene ID
        if tmpData['sequence_type'] == 'gene_id':
            logger.debug( "Selected NCBI Gene ID. Retrieving sequence from NCBI...")
            url1 = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=protein&term=srcdb_refseq[PROP]&id=' + re.sub('[,\s]+', '', tmpData['sequence'])
            with urllib.request.urlopen(url1) as file1:
                geneIDs = list(set(re.findall('<Link>\s+<Id>(\d+)<\/Id>\s+<\/Link>', file1.read().decode())))
            url2 = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=text&rettype=fasta&id=' + ','.join(geneIDs)
            with urllib.request.urlopen(url2) as file2:
                tmpData['sequence'] = file2.read().decode()
        # Sequence input type: UniprotKB Protein ID
        if tmpData['sequence_type'] == 'protein_uniprotkb':
            logger.debug("Selected UniprotKB Protein ID. Retrieving sequence from Uniprot...")
            url = str.format('https://www.uniprot.org/uniprot/?query=id:{}&format=fasta', re.sub('[,\s]+', '+OR+id:', tmpData['sequence']))
            with urllib.request.urlopen(url) as file:
                tmpData['sequence'] = file.read().decode()
        # Sequence input type: NCBI Bioproject ID
        if tmpData['sequence_type'] == 'bioproject_id':
            logger.debug( "Selected NCBI Bioproject ID. Retrieving sequence from NCBI...")
            if tmpData['sequence'].strip().startswith( 'PRJNA' ):
                bioprojectID = tmpData['sequence'].strip().replace('PRJNA', '')
            else:
                bioprojectID = tmpData['sequence'].strip()
            url1 = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=bioproject&db=protein&id='+bioprojectID
            proteinIDs = []
            with urllib.request.urlopen(url1) as file1:
                tree1 = ET.fromstring(file1.read().decode())
                for proteinID in tree1.find('LinkSet').find('LinkSetDb').iter('Id'):
                    proteinIDs.append(proteinID.text)
            url2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            data2 = urllib.parse.urlencode({
                'db': 'protein',
                'rettype': 'fasta',
                'retmode': 'text',
                'retmax': '10000',
                'id': ','.join(proteinIDs),
            }).encode()
            request2 = urllib.request.Request(url=url2, data=data2)
            with urllib.request.urlopen(request2) as file2:
                tmpData['sequence'] = file2.read().decode()
            
            url3 = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=bioproject&id='+bioprojectID
            with urllib.request.urlopen(url3) as file3:
                tree3 = ET.fromstring(file3.read().decode())
                organism = tree3.find('DocumentSummarySet').find('DocumentSummary').find('Organism_Name').text
                strain = tree3.find('DocumentSummarySet').find('DocumentSummary').find('Organism_Strain').text
                if tmpData['group_name'] == '':
                    tmpData['group_name'] = organism
                if tmpData['genome_name'] == '':
                    tmpData['genome_name'] = str.format('{} {}', organism, strain)
        # Sequence input type: Uniprot Proteome ID
        if tmpData['sequence_type'] == 'protein_uniprot_proteome':
            logger.debug("Selected Uniprot Proteome ID. Retrieving sequence from Uniprot...")
            url = str.format('https://www.uniprot.org/uniprot/?query=proteome:{}&format=fasta', tmpData['sequence'].strip())
            with urllib.request.urlopen(url) as file:
                tmpData['sequence'] = file.read().decode()
            
        form = RunsForm(tmpData or None)
    else:
        form = RunsForm()
    if request.method == 'POST':
        if form.is_valid():
            formData = form.cleaned_data
            queryID = ''.join(random.choice('ABCDEFGHJKMNPQRSTUVWXYZ23456789') for i in range(10))
            formData['c_projectfk'] = 0
            formData['group_name'] = ''
            formData['genome_name'] = ''
            request.session['runs_form_data'] = formData
            return redirect('/vaxign2/run/'+queryID)
    context['form'] = form
    
    return render(request, 'vaxignml.html', context)

def vaxitop(request):
    
    context = {}
    
    if 'Submit' in request.POST:
        tmpData = {}
        for key in request.POST:
            tmpData[key] = request.POST[key]
        if 'sequence_file' in request.FILES:
            file = request.FILES['sequence_file']
            fileContent = ''
            for chunk in file.chunks():
                fileContent += chunk.decode()
            tmpData['sequence'] = fileContent.strip()
            tmpData['sequence_type'] = 'protein_fasta'
            
        # Process sequence data
        # Sequence input type: NCBI Protein GI or Refseq
        if tmpData['sequence_type'] == 'protein_gi':
            logger.debug( "Selected NCBI Protein GI. Retrieving sequence from NCBI...")
            try:
                url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=text&rettype=fasta&id=' + re.sub('[,\s]+', ',', tmpData['sequence'])
                with urllib.request.urlopen(url) as file:
                    tmpData['sequence'] = file.read().decode()
            except:
                tmpData['sequence'] = None
        # Sequence input type: NCBI Protein Refseq
        if tmpData['sequence_type'] == 'protein_refseq':
            logger.debug( "Selected NCBI Protein Refseq. Retrieving sequence from NCBI...")
            try:
                url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=text&rettype=fasta&id=' + re.sub('[,\s]+', ',', tmpData['sequence'])
                with urllib.request.urlopen(url) as file:
                    tmpData['sequence'] = file.read().decode()
            except:
                tmpData['sequence'] = None
        # Sequence input type: FASTA URL
        if tmpData['sequence_type'] == 'protein_fasta_url':
            logger.debug( "Selected custom protein FASTA file link(s). Retrieving sequence from the source(s)...")
            tmp = ''
            for url in re.split('[\r\n]+', tmpData['sequence']):
                with urllib.request.urlopen(url.strip()) as file:
                    tmp += file.read().decode() + '\n'
            tmpData['sequence'] = tmp
        # Sequence input type: NCBI Gene ID
        if tmpData['sequence_type'] == 'gene_id':
            logger.debug( "Selected NCBI Gene ID. Retrieving sequence from NCBI...")
            url1 = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=protein&term=srcdb_refseq[PROP]&id=' + re.sub('[,\s]+', '', tmpData['sequence'])
            with urllib.request.urlopen(url1) as file1:
                geneIDs = list(set(re.findall('<Link>\s+<Id>(\d+)<\/Id>\s+<\/Link>', file1.read().decode())))
            url2 = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&retmode=text&rettype=fasta&id=' + ','.join(geneIDs)
            with urllib.request.urlopen(url2) as file2:
                tmpData['sequence'] = file2.read().decode()
        # Sequence input type: UniprotKB Protein ID
        if tmpData['sequence_type'] == 'protein_uniprotkb':
            logger.debug("Selected UniprotKB Protein ID. Retrieving sequence from Uniprot...")
            url = str.format('https://www.uniprot.org/uniprot/?query=id:{}&format=fasta', re.sub('[,\s]+', '+OR+id:', tmpData['sequence']))
            with urllib.request.urlopen(url) as file:
                tmpData['sequence'] = file.read().decode()
        # Sequence input type: NCBI Bioproject ID
        if tmpData['sequence_type'] == 'bioproject_id':
            logger.debug( "Selected NCBI Bioproject ID. Retrieving sequence from NCBI...")
            if tmpData['sequence'].strip().startswith( 'PRJNA' ):
                bioprojectID = tmpData['sequence'].strip().replace('PRJNA', '')
            else:
                bioprojectID = tmpData['sequence'].strip()
            url1 = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=bioproject&db=protein&id='+bioprojectID
            proteinIDs = []
            with urllib.request.urlopen(url1) as file1:
                tree1 = ET.fromstring(file1.read().decode())
                for proteinID in tree1.find('LinkSet').find('LinkSetDb').iter('Id'):
                    proteinIDs.append(proteinID.text)
            url2 = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            data2 = urllib.parse.urlencode({
                'db': 'protein',
                'rettype': 'fasta',
                'retmode': 'text',
                'retmax': '10000',
                'id': ','.join(proteinIDs),
            }).encode()
            request2 = urllib.request.Request(url=url2, data=data2)
            with urllib.request.urlopen(request2) as file2:
                tmpData['sequence'] = file2.read().decode()
            
            url3 = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=bioproject&id='+bioprojectID
            with urllib.request.urlopen(url3) as file3:
                tree3 = ET.fromstring(file3.read().decode())
                organism = tree3.find('DocumentSummarySet').find('DocumentSummary').find('Organism_Name').text
                strain = tree3.find('DocumentSummarySet').find('DocumentSummary').find('Organism_Strain').text
                if tmpData['group_name'] == '':
                    tmpData['group_name'] = organism
                if tmpData['genome_name'] == '':
                    tmpData['genome_name'] = str.format('{} {}', organism, strain)
        # Sequence input type: Uniprot Proteome ID
        if tmpData['sequence_type'] == 'protein_uniprot_proteome':
            logger.debug("Selected Uniprot Proteome ID. Retrieving sequence from Uniprot...")
            url = str.format('https://www.uniprot.org/uniprot/?query=proteome:{}&format=fasta', tmpData['sequence'].strip())
            with urllib.request.urlopen(url) as file:
                tmpData['sequence'] = file.read().decode()
        
        if tmpData['sequence'] != '':
            queryID = ''.join(random.choice('ABCDEFGHJKMNPQRSTUVWXYZ23456789') for i in range(10))
            request.session['sequence'] = tmpData['sequence']

            return redirect('/vaxign2/run/'+queryID+'/vaxitop')
        else:
            form = RunsForm(tmpData or None)
    else:
        form = RunsForm()
    context['form'] = form
    
    return render(request, 'vaxitop.html', context)
