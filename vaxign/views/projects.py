from django.shortcuts import render, redirect

# Create your views here.

import re
import time
import json
import random
import urllib
import string
import collections
import xml.etree.ElementTree as ET

from django.db.models import Q
from django.db.models.functions import Now
from vaxign.models import TUserVaxignPrivilege
from vaxign.models import TVaxignProject
from vaxign.models import TUser
from vaxign.models import TVaxignQuery
from vaxign.models import TVaxignAnalysisResults

from django.forms import model_to_dict
from vaxign.forms.projects import ProjectsForm
from vaxign.forms.runs import RunsForm

from django.http import HttpResponse

import logging
logger = logging.getLogger('console')

def index(request):
    
    if 'c_user_name' not in request.session:
        return redirect('/users/login.php?redirect=/vaxign2/login')
    elif request.session['c_user_name'] not in TUser.objects.all().values_list('c_user_name', flat=True):
        return HttpResponse(status=403)
    
    context = {}
    
    user = request.session['c_user_name']
    if request.session['is_admin']:
        projectIDs = TVaxignProject.objects.all().values_list('c_vaxign_projectid', flat=True)
    else:
        projectIDs = TUserVaxignPrivilege.objects.filter(
            Q(c_userfk=user)
        ).values_list('c_projectfk', flat=True)
    
    curators = {}
    for row in TUserVaxignPrivilege.objects.filter(Q(c_projectfk__in=projectIDs)):
        if row.c_projectfk not in curators.keys():
            curators[row.c_projectfk] = []
        try:
            others = TUser.objects.get(c_user_name=row.c_userfk)
            curators[row.c_projectfk].append(str.format("{} {}", others.c_first_name, others.c_last_name).title())
        except:
            pass
    
    projects = []
    for row in TVaxignProject.objects.filter(
        Q(c_vaxign_projectid__in=projectIDs)
    ).order_by('-c_vaxign_projectid'):
        tmp = row.__dict__
        if row.c_vaxign_projectid in curators.keys():
            tmp['curators'] = ', '.join(curators[row.c_vaxign_projectid])
        else:
            tmp['curators'] = ''
        tmp['format_time'] = row.c_creation_time.strftime('%Y-%m-%d %H:%M')
        projects.append(tmp)
    context['projects'] = projects
    
    return render(request, 'projects/projects.html', context)

def add(request):
    
    if 'c_user_name' not in request.session:
        return redirect('/users/register.php?redirect=/vaxign2/login')
    elif request.session['c_user_name'] not in TUser.objects.all().values_list('c_user_name', flat=True):
        return HttpResponse(status=403)
    
    context = {}
    
    form = ProjectsForm(request.POST or None)
    
    if request.method == 'POST':
        if form.is_valid():
            formData = form.cleaned_data
            
            project = TVaxignProject(
                c_vaxign_projectname = formData['project_name'],
                c_description = formData['description'],
                c_note = formData['note'],
                c_lastupdate = Now(),
                c_isopen = formData['is_public'],
                c_authorized = 1,
                c_userfk = request.session['c_user_name'],
                c_institution = formData['institution'],
                c_investigator = formData['investigator'],
                c_list_title = formData['list_project'],
                c_creation_time = Now(),
                c_is_locked = 0,
                c_show_settings = 0,
                c_grant = formData['grant'],
            )
            project.save()
            
            permission = TUserVaxignPrivilege(
                c_userfk = request.session['c_user_name'],
                c_projectfk = project.c_vaxign_projectid,
                c_permission = 'update',
            )
            permission.save()
            
            return redirect('/vaxign2/project/'+str(project.c_vaxign_projectid))
            
    context['form'] = form
    
    return render(request, 'projects/add.html', context)

def open(request, projectID):
    
    if 'c_user_name' not in request.session:
        return redirect('/users/register.php?redirect=/vaxign2/login')
    elif request.session['c_user_name'] not in TUser.objects.all().values_list('c_user_name', flat=True):
        return HttpResponse(status=403)
    
    context = {}
    
    if int(projectID) in TVaxignProject.objects.all().values_list('c_vaxign_projectid', flat=True):
        project = TVaxignProject.objects.get(c_vaxign_projectid=projectID)
        context['project'] = project
        
        tmpQueries = TVaxignQuery.objects.filter(c_projectfk=projectID).order_by('c_submit_time')
        queries = []
        for tmpQuery in tmpQueries:
            query = tmpQuery.__dict__
            query['c_submit_time'] = tmpQuery.c_submit_time.strftime('%Y-%m-%d %H:%M')
            if tmpQuery.c_finish_time:
                query['c_finish_time'] = tmpQuery.c_finish_time.strftime('%Y-%m-%d %H:%M')
            else:
                query['c_finish_time'] = ''
            query['sequence_count'] = TVaxignAnalysisResults.objects.filter(c_query_id=tmpQuery.c_query_id).count()
            queries.append(query)
        context['queries'] = queries
        
        owner = TUser.objects.get(c_user_name=project.c_userfk)
        context['owner'] = owner
        
        curators = TUser.objects.filter(
            Q(c_user_name__in=TUserVaxignPrivilege.objects.filter(c_projectfk=projectID).values_list('c_userfk', flat=True)) &
            ~Q(c_user_name=owner.c_user_name)
        )
        context['curators'] = curators
        
        return render(request, 'projects/open.html', context)
    else:
        return redirect('/vaxign2/project')

def edit(request, projectID):
    
    if 'c_user_name' not in request.session:
        return redirect('/users/register.php?redirect=/vaxign2/login')
    
    context = {}
    
    if request.method == 'POST':
        if request.POST['project_id'] != projectID:
            return redirect('/vaxign2/project/'+str(projectID))
        
        form = ProjectsForm(request.POST)
        if form.is_valid():
            formData = form.cleaned_data
            
            project = TVaxignProject.objects.get(c_vaxign_projectid=projectID)
            
            project.c_vaxign_projectname = formData['project_name']
            project.c_description = formData['description']
            project.c_note = formData['note']
            project.c_lastupdate = Now()
            project.c_isopen = formData['is_public']
            project.c_institution = formData['institution']
            project.c_investigator = formData['investigator']
            project.c_list_title = formData['list_project']
            project.c_grant = formData['grant']
            project.save()
            
            return redirect('/vaxign2/project/'+str(project.c_vaxign_projectid))
    else:
        project = TVaxignProject.objects.get(c_vaxign_projectid=projectID)
        formData = {
            'project_id': projectID,
            'project_name': project.c_vaxign_projectname,
            'is_public': project.c_isopen,
            'list_project': project.c_list_title,
            'description': project.c_description,
            'institution': project.c_institution,
            'investigator': project.c_investigator,
            'grant': project.c_grant,
            'note': project.c_note,
        }
        
        form = ProjectsForm(formData)
        
    context['form'] = form
    
    return render(request, 'projects/edit.html', context)

def remove(request, projectID):
    
    if 'c_user_name' not in request.session:
        return redirect('/users/register.php?redirect=/vaxign2/login')
    elif request.session['c_user_name'] not in TUser.objects.all().values_list('c_user_name', flat=True):
        return HttpResponse(status=403)
    
    try:
        permission = TUserVaxignPrivilege.objects.filter(
            Q(c_projectfk=projectID)
        )
        project = TVaxignProject.objects.get(
            Q(c_userfk=request.session['c_user_name']) & Q(c_vaxign_projectid=projectID)
        )
    except:
        return HttpResponse(status=400)
    
    permission.delete()
    project.delete()
    return HttpResponse("Success.")

def curator(request, projectID, email, type):
    
    if 'c_user_name' not in request.session:
        return redirect('/users/register.php?redirect=/vaxign2/login')
    elif request.session['c_user_name'] not in TUser.objects.all().values_list('c_user_name', flat=True):
        return HttpResponse(status=403)
    
    if type == 'add':
        users = TUser.objects.all().values_list('c_user_name', flat=True)
        if email not in users:
            return HttpResponse(status=403)
        permission = TUserVaxignPrivilege(
            c_userfk = email,
            c_projectfk = projectID,
            c_permission = 'update',
        )
        permission.save()
        return HttpResponse("Success.")
    elif type == 'remove':
        TUserVaxignPrivilege.objects.get(
            Q(c_userfk=email) & Q(c_projectfk=projectID)
        ).delete()
        return HttpResponse("Success.")
    else:
        return HttpResponse(status=403)

def run(request, projectID):
    
    if 'c_user_name' not in request.session:
        return redirect('/users/register.php?redirect=/vaxign2/login')
    elif request.session['c_user_name'] not in TUser.objects.all().values_list('c_user_name', flat=True):
        return HttpResponse(status=403)
    
    context = {}
    
    project = TVaxignProject.objects.get(c_vaxign_projectid=projectID)
    context['project'] = project
    
    if request.method == 'POST':
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
            formData['c_projectfk'] = projectID
            request.session['runs_form_data'] = formData
            return redirect('/vaxign2/run/'+queryID)
    context['form'] = form
        
    return render(request, "projects/run.html", context)
    
def querySetting(request, projectID, queryID):
    
    if 'c_user_name' not in request.session:
        return redirect('/users/register.php?redirect=/vaxign2/login')
    elif request.session['c_user_name'] not in TUser.objects.all().values_list('c_user_name', flat=True):
        return HttpResponse(status=403)
    
    try:
        query = TVaxignQuery.objects.get(
            Q(c_query_id=queryID) & Q(c_projectfk=projectID)
        )
    except:
        return HttpResponse(status=400)
    
    context = {}
    
    project = TVaxignProject.objects.get(c_vaxign_projectid=projectID)
    context['project'] = project
    
    organism = ''
    if query.c_organism is None or query.c_organism == '':
        if query.c_gram == '--positive':
            organism = 'Gram Positive Bacteria'
        elif query.c_gram == '--negative':
            organism = 'Gram Negative Bacteria'
    else:
        if query.c_organism == 'bacteria':
            if query.c_gram == '--positive':
                organism = 'Gram Positive Bacteria'
            elif query.c_gram == '--negative':
                organism = 'Gram Negative Bacteria'
        else:
            organism = query.c_organism
    
    basic_analyses = []
    vaxign_ml = False
    mhc_i = False
    mhc_ii = False
    for analysis in query.c_included_analyses.split(','):
        if analysis == 'psort':
            basic_analyses.append('Subcellular Localization')
        if analysis == 'spaan':
            basic_analyses.append('Adhesion Probability')
        if analysis == 'tmhmm':
            basic_analyses.append('Transmembrane Helices')
        if analysis == 'blasth':
            basic_analyses.append('Similarity to Human Proteins')
        if analysis == 'blastm':
            basic_analyses.append('Similarity to Mouse Proteins')
        if analysis == 'blastp':
            basic_analyses.append('Similarity to Pig Proteins')
        if analysis == 'vaxign_ml':
            vaxign_ml = True
        if analysis == 'mast_I':
            mhc_i = True
        if analysis == 'mast_II':
            mhc_ii = True
    
    query_detail = {
        'c_query_id': query.c_query_id,
        'organism': organism,
        'basic_analyses': basic_analyses,
        'vaxign_ml' : vaxign_ml,
        'mhc_i': mhc_i,
        'mhc_ii': mhc_ii,
    }
    context['query'] = query_detail
    
    return render(request, "projects/setting.html", context)
    
def queryRemove(request, projectID, queryID):
    
    if 'c_user_name' not in request.session:
        return redirect('/users/register.php?redirect=/vaxign2/login')
    elif request.session['c_user_name'] not in TUser.objects.all().values_list('c_user_name', flat=True):
        return HttpResponse(status=403)
    
    try:
        query = TVaxignQuery.objects.get(
            Q(c_query_id=queryID) & Q(c_projectfk=projectID)
        )
    except:
        return HttpResponse(status=400)
    
    query.delete()
    
    return HttpResponse("Success.")

def queryPermission(request, projectID, queryID, type):
    
    if 'c_user_name' not in request.session:
        return redirect('/users/register.php?redirect=/vaxign2/login')
    elif request.session['c_user_name'] not in TUser.objects.all().values_list('c_user_name', flat=True):
        return HttpResponse(status=403)
    
    query = TVaxignQuery.objects.get(
        Q(c_query_id=queryID) & Q(c_projectfk=projectID)
    )
    query.c_is_public = type
    query.save()
    
    return HttpResponse("Success.")
