from django.shortcuts import render, redirect

# Create your views here.

import os
import re
import urllib
import subprocess
import pandas as pd
import plotly.graph_objects as go
from Bio import Entrez
from Bio.KEGG import REST
from goatools.base import get_godag
from plotly import offline

from django.conf import settings

from django.db.models import Q
from vaxign.models import TVaxignAnalysisResults
from vaxign.models import TVaxignEggnogFunctions
from vaxign.models import TVaxignEggnogOrtholog
from vaxign.models import TVaxignAlleleGroup
from vaxign.models import TVaxignMastResults
from vaxign.models import TVaxignPopulationCoverage

from vaxign.views.runs import BulkCreateManager

import logging
from django.core.exceptions import ValidationError
logger = logging.getLogger('console')

def protein_eggnog_function(request, queryID, seqID):
    
    context = {'query_id': queryID, 'sequence_id': seqID}
    
    try:
        sequence = TVaxignAnalysisResults.objects.get(c_sequence_id=seqID)
        context['sequence'] = sequence.c_sequence
    except:
        return HttpResponse("No eggNOG Prediction available")
    
    if sequence.c_eggnog_function_computed == None or sequence.c_eggnog_function_computed == 'n':
        if not os.path.exists(os.path.join(settings.WORKSPACE_DIR, queryID, queryID+'.fasta')):
            if not os.path.exists(settings.VAXIGN2_TMP_DIR):
                os.mkdir(settings.VAXIGN2_TMP_DIR)
            if not os.path.exists(os.path.join(settings.VAXIGN2_TMP_DIR, queryID)):
                os.mkdir(os.path.join(settings.VAXIGN2_TMP_DIR, queryID))
                
            queryPath = os.path.join(settings.VAXIGN2_TMP_DIR, queryID)
            open(os.path.join(queryPath, seqID+'.fasta'), 'w').write(str.format(""">{}
{}""", sequence.c_sequence_id, sequence.c_sequence))
            cmd = ['python2', os.path.join(settings.EGGNOG_PATH, 'emapper.py'),
                   '--cpu', '10',
                   '-i', os.path.join(queryPath, seqID+'.fasta'),
                   '--output', seqID,
                   '--output_dir', queryPath,
                   '--temp_dir', queryPath,
                   '--data_dir', settings.EGGNOG_DB_PATH,
                   '-m', 'diamond',
                   '-d', 'none',
                   '--tax_scope', 'auto',
                   '--go_evidence', 'non-electronic',
                   '--target_orthologs', 'all',
                   '--seed_ortholog_evalue', '0.001',
                   '--seed_ortholog_score', '60',
                   '--query-cover', '20',
                   '--subject-cover', '0',
                   '--override', '--no_file_comments', '--predict_ortho'
                   ]
            subprocess.call(cmd)
            for line in open(os.path.join(queryPath, seqID+'.emapper.annotations')).read().splitlines():
                if line == '':
                    continue
                tokens = line.split('\t')
                TVaxignEggnogFunctions(
                    c_sequence_id= tokens[0],
                    seed_eggnog_ortholog= tokens[1],
                    seed_ortholog_evalue= tokens[2],
                    seed_ortholog_score= tokens[3],
                    best_tax_level= tokens[4],
                    preferred_name= tokens[5],
                    gos=  tokens[6],
                    ec=  tokens[7],
                    kegg_ko=  tokens[8],
                    kegg_pathway=  tokens[9],
                    kegg_module=  tokens[10],
                    kegg_reaction=  tokens[11],
                    kegg_rclass=  tokens[12],
                    brite=  tokens[13],
                    kegg_tc=  tokens[14],
                    cazy=  tokens[15],
                    bigg_reaction=  tokens[16],
                    cogs=  tokens[-2],
                ).save()
            sequence.c_eggnog_function_computed = 'y'
            sequence.save()
                
        else:
            queryPath = os.path.join(settings.WORKSPACE_DIR, queryID)
            cmd = ['python2', os.path.join(settings.EGGNOG_PATH, 'emapper.py'),
                   '--cpu', '10',
                   '-i', os.path.join(queryPath, queryID+'.fasta'),
                   '--output', queryID,
                   '--output_dir', queryPath,
                   '--temp_dir', queryPath,
                   '--data_dir', settings.EGGNOG_DB_PATH,
                   '-m', 'diamond',
                   '-d', 'none',
                   '--tax_scope', 'auto',
                   '--go_evidence', 'non-electronic',
                   '--target_orthologs', 'all',
                   '--seed_ortholog_evalue', '0.001',
                   '--seed_ortholog_score', '60',
                   '--query-cover', '20',
                   '--subject-cover', '0',
                   '--override', '--no_file_comments', '--predict_ortho'
                   ]
            subprocess.call(cmd)
            for line in open(os.path.join(queryPath, queryID+'.emapper.annotations')).read().splitlines():
                if line == '':
                    continue
                tokens = line.split('\t')
                TVaxignEggnogFunctions(
                    c_sequence_id= tokens[0],
                    seed_eggnog_ortholog= tokens[1],
                    seed_ortholog_evalue= tokens[2],
                    seed_ortholog_score= tokens[3],
                    best_tax_level= tokens[4],
                    preferred_name= tokens[5],
                    gos=  tokens[6],
                    ec=  tokens[7],
                    kegg_ko=  tokens[8],
                    kegg_pathway=  tokens[9],
                    kegg_module=  tokens[10],
                    kegg_reaction=  tokens[11],
                    kegg_rclass=  tokens[12],
                    brite=  tokens[13],
                    kegg_tc=  tokens[14],
                    cazy=  tokens[15],
                    bigg_reaction=  tokens[16],
                    cogs=  tokens[-2],
                ).save()
            sequence.c_eggnog_function_computed = 'y'
            sequence.save()
    try:    
        eggnog = TVaxignEggnogFunctions.objects.get(c_sequence_id=seqID)
        
        godag = get_godag(settings.GO_BASIC_OBO_PATH)
        gos = {
            'biological_process': {},
            'molecular_function': {},
            'cellular_component': {},
        }
        for goID in eggnog.gos.split(','):
            if goID not in godag.keys():
                continue
            go = godag[goID]
            gos[go.namespace][goID] = go.name
        context['gos'] = gos
        
        cogs = {}
        for cog in eggnog.cogs:
            cogs[cog] = settings.COG_MAP[cog]
        context['cogs'] = cogs
        
        keggs = {
            'kegg_pathway': {},
            'kegg_ko': {},
        }
        for kegg_pathway in eggnog.kegg_pathway.split(','):
            keggs['kegg_pathway'][kegg_pathway] = re.split('\s+', list(REST.kegg_get(kegg_pathway))[1],maxsplit=1)[1].strip()
        for kegg_ko in eggnog.kegg_ko.split(','):
            keggs['kegg_ko'][kegg_ko] = re.split('\s+', list(REST.kegg_get(kegg_ko))[2],maxsplit=1)[1].strip()
        context['keggs'] = keggs
    except:
        pass
    
    return render(request, 'queries/tabs/eggnog_function.html', context)

def protein_eggnog_ortholog(request, queryID, seqID, display):
    
    context = {'query_id': queryID, 'sequence_id': seqID}
    
    try:
        sequence = TVaxignAnalysisResults.objects.get(c_sequence_id=seqID)
        context['sequence'] = sequence.c_sequence
    except:
        return HttpResponse("No eggNOG Prediction available")
    
    if sequence.c_eggnog_ortholog_computed == None or sequence.c_eggnog_ortholog_computed == 'n':
        if not os.path.exists(os.path.join(settings.WORKSPACE_DIR, queryID, queryID+'.fasta')):
            if not os.path.exists(settings.VAXIGN2_TMP_DIR):
                os.mkdir(settings.VAXIGN2_TMP_DIR)
            if not os.path.exists(os.path.join(settings.VAXIGN2_TMP_DIR, queryID)):
                os.mkdir(os.path.join(settings.VAXIGN2_TMP_DIR, queryID))
                
            queryPath = os.path.join(settings.VAXIGN2_TMP_DIR, queryID)
            open(os.path.join(queryPath, seqID+'.fasta'), 'w').write(str.format(""">{}
{}""", sequence.c_sequence_id, sequence.c_sequence))
            cmd = ['python2', os.path.join(settings.EGGNOG_PATH, 'emapper.py'),
               '--cpu', '10',
               '-i', os.path.join(queryPath, seqID+'.fasta'),
               '--output', seqID,
               '--output_dir', queryPath,
               '--temp_dir', queryPath,
               '--data_dir', settings.EGGNOG_DB_PATH,
               '-m', 'diamond',
               '-d', 'none',
               '--tax_scope', 'auto',
               '--go_evidence', 'non-electronic',
               '--target_orthologs', 'all',
               '--seed_ortholog_evalue', '0.001',
               '--seed_ortholog_score', '60',
               '--query-cover', '20',
               '--subject-cover', '0',
               '--override', '--no_file_comments', '--predict_ortho'
               ]
            subprocess.call(cmd)
            for line in open(os.path.join(queryPath, seqID+'.emapper.predict_orthologs')).read().splitlines():
                if line.startswith('#') or line == '':
                    continue
                tokens = line.split('\t')
                TVaxignEggnogOrtholog(
                    c_sequence_id = tokens[0],
                    c_ortholog_taxon = tokens[1],
                    c_ortholog_gene = tokens[2],
                ).save()
            sequence.c_eggnog_ortholog_computed = 'y'
            sequence.save()
        else:
            queryPath = os.path.join(settings.WORKSPACE_DIR, queryID)
            cmd = ['python2', os.path.join(settings.EGGNOG_PATH, 'emapper.py'),
                   '--cpu', '10',
                   '-i', os.path.join(queryPath, queryID+'.fasta'),
                   '--output', queryID,
                   '--output_dir', queryPath,
                   '--temp_dir', queryPath,
                   '--data_dir', settings.EGGNOG_DB_PATH,
                   '-m', 'diamond',
                   '-d', 'none',
                   '--tax_scope', 'auto',
                   '--go_evidence', 'non-electronic',
                   '--target_orthologs', 'all',
                   '--seed_ortholog_evalue', '0.001',
                   '--seed_ortholog_score', '60',
                   '--query-cover', '20',
                   '--subject-cover', '0',
                   '--override', '--no_file_comments', '--predict_ortho'
                   ]
            subprocess.call(cmd)
            for line in open(os.path.join(queryPath, queryID+'.emapper.predict_orthologs')).read().splitlines():
                if line.startswith('#') or line == '':
                    continue
                tokens = line.split('\t')
                TVaxignEggnogOrtholog(
                    c_sequence_id = tokens[0],
                    c_ortholog_taxon = tokens[1],
                    c_ortholog_gene = tokens[2],
                ).save()
            sequence.c_eggnog_ortholog_computed = 'y'
            sequence.save()
    
    eggnog = TVaxignEggnogOrtholog.objects.filter(c_sequence_id=seqID)
    
    if eggnog.count() != 0:
        Entrez.email = settings.ADMIN_EMAIL
        tops = ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses', 'Other', 'Unclassified']
        
        taxonMap = {}
        treeMap = {}
        for record in Entrez.read(Entrez.efetch(db='Taxonomy', id=[str(i) for i in list(eggnog.values_list('c_ortholog_taxon', flat=True))])):
            taxonMap[record['TaxId']] = record['ScientificName']
            topFound = False
            parent = None
            top = None
            for lineage in record['Lineage'].split('; '):
                if not topFound and lineage not in tops:
                    continue
                elif not topFound:
                    parent = lineage
                    topFound = True
                    continue
                if parent not in treeMap.keys():
                    treeMap[parent] = set([])
                treeMap[parent].add(lineage)
                parent = lineage
            treeMap[parent] = set([])
            treeMap[parent].add(record['TaxId'])
        
        rawGeneNames = []
        for rawGenes in eggnog.values_list('c_ortholog_gene', flat=True):
            for rawGene in rawGenes.split(','):
                rawGeneNames.append(rawGene.split('.')[1])
        url = 'https://www.uniprot.org/uploadlists/'
        data = urllib.parse.urlencode({
            'from': 'GENENAME',
            'to': 'ID',
            'format': 'tab',
            'query': ' '.join(rawGeneNames)
        }).encode()
        proteinMap = {}
        with urllib.request.urlopen(url, data=data) as file:
            for line in file.read().decode().splitlines()[1:]:
                tokens = line.split('\t')
                proteinMap[tokens[0]] = tokens[1]
    
    if display == 'tree':
        orthologs = {}
        for ortholog in eggnog:
            taxon = ortholog.c_ortholog_taxon
            if str(taxon) not in taxonMap.keys():
                taxonName = ''
            else:
                taxonName = taxonMap[str(taxon)]
            proteins = []
            for gene in ortholog.c_ortholog_gene.split(','):
                geneID = gene.split('.')[1]
                if geneID not in proteinMap:
                    proteins.append(geneID)
                else:
                    proteins.append(proteinMap[geneID])
            if str(taxon) not in orthologs.keys():
                orthologs[str(taxon)] = []
            orthologs[str(taxon)] = proteins
        
        html = "<ul>"
        for top in tops:
            if top not in treeMap.keys():
                continue
            html += "<li class='jstree-open' data-jstree='{\"icon\":\"cell\"}'>"+top+"<ul>"
            queue = list(treeMap[top])
            branch = list(treeMap[top])
            while len(queue) != 0:
                current = queue.pop()
                if current in treeMap.keys():
                    html += "<li class='jstree-open' data-jstree='{\"icon\":\"cell\"}'>"+current+"<ul>"
                    queue += list(treeMap[current])
                    branch += list(treeMap[current])
                    del treeMap[current]
                else:
                    html += str.format("<li class='jstree-open' data-jstree='{{\"icon\":\"cell\"}}'><a class='has_taxon'href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={}' target='_blank'>{} (Taxon:{})</a><ul>",
                                       current, taxonMap[current], current)
                    if current in orthologs.keys():
                        for protein in orthologs[current]:
                            html += str.format("<li class='jstree-open' data-jstree='{{\"icon\":\"protein\"}}'><a href='https://www.uniprot.org/uniprot/{}'>{}</a></li>", protein, protein)
                    html += "</ul></li>"
                    while len(branch) > 0 and branch[-1] not in queue:
                        branch.pop()
                        if len(branch) > 0 and branch[-1] not in queue:
                            html += "</ul></li>"
            html += "</ul></li>"
        html += "</ul>"
        context['orthologs'] = orthologs
        context['tree'] = html
        
        return render(request, 'queries/tabs/eggnog_ortholog_tree.html', context)
    else:
        orthologs = []
        for ortholog in eggnog:
            taxon = ortholog.c_ortholog_taxon
            if str(taxon) not in taxonMap.keys():
                taxonName = ''
            else:
                taxonName = taxonMap[str(taxon)]
            proteins = []
            for gene in ortholog.c_ortholog_gene.split(','):
                geneID = gene.split('.')[1]
                if geneID not in proteinMap:
                    proteins.append(geneID)
                else:
                    proteins.append(proteinMap[geneID])
            orthologs.append({
                'taxon_id': taxon,
                'taxon_name': taxonName,
                'proteins': proteins,
            })
        context['orthologs'] = orthologs
        return render(request, 'queries/tabs/eggnog_ortholog_table.html', context)

def protein_population_coverage(request, queryID, seqID, mhc_class, country_code=None):
    
    if country_code == None:
        country_code = settings.POPULATION_COVERAGE_ISO_CODE.values()
    else:
        country_code = country_code.split(',')
    
    context = {'query_id': queryID, 'sequence_id': seqID}
    
    try:
        sequence = TVaxignAnalysisResults.objects.get(c_sequence_id=seqID)
        context['sequence'] = sequence.c_sequence
    except:
        return HttpResponse("No Vaxitop Prediction available")
    
    if TVaxignMastResults.objects.filter(c_sequence_id=seqID).count() == 0:
        if not os.path.exists(os.path.join(settings.WORKSPACE_DIR, queryID, queryID+'.fasta')):
            if not os.path.exists(settings.VAXIGN2_TMP_DIR):
                os.mkdir(settings.VAXIGN2_TMP_DIR)
            if not os.path.exists(os.path.join(settings.VAXIGN2_TMP_DIR, queryID)):
                os.mkdir(os.path.join(settings.VAXIGN2_TMP_DIR, queryID))
                
            queryPath = os.path.join(settings.VAXIGN2_TMP_DIR, queryID)
            open(os.path.join(queryPath, seqID+'.fasta'), 'w').write(str.format(""">{}
{}""", sequence.c_sequence_id, sequence.c_sequence))
            alleleGroups = TVaxignAlleleGroup.objects.all().values_list('c_allele_group_id', flat=True)
            cmds = []
            for group in alleleGroups:
                cmd = str.format("{} {} {} -mt 0.1 -hit_list > {}",
                                    os.path.join(settings.VAXITOP_PATH, 'lib', 'meme', 'bin', 'mast'),
                                    os.path.join(settings.VAXITOP_PATH, 'pssm', str(group)+'.xml'),
                                    os.path.join(queryPath, seqID+'.fasta'),
                                    os.path.join(queryPath, str(group)+'.matching.txt'),
                                )
                cmds.append(cmd)
        else:
            queryPath = os.path.join(settings.WORKSPACE_DIR, queryID)
            cmds = []
            for group in alleleGroups:
                cmd = str.format("{} {} {} -mt 0.1 -hit_list > {}",
                                    os.path.join(settings.VAXITOP_PATH, 'lib', 'meme', 'bin', 'mast'),
                                    os.path.join(settings.VAXITOP_PATH, 'pssm', str(group)+'.xml'),
                                    os.path.join(queryPath, queryID+'.fasta'),
                                    os.path.join(queryPath, str(group)+'.matching.txt'),
                                )
                cmds.append(cmd)

        _multiprocess(cmds, 10)
        bulk = BulkCreateManager()
        for group in alleleGroups:
            for line in open(os.path.join(queryPath, str(group)+'.matching.txt')):
                if line.startswith('#'):
                    continue
                tokens = re.split('[ ]+', line)
                mast = TVaxignMastResults()
                mast.c_sequence_id = tokens[0]
                mast.c_allele_group_id = group
                mast.c_hit_start = tokens[4]
                mast.c_hit_end = tokens[5]
                mast.c_hit_p_value = tokens[-1]
                mast.c_hit_z_value = norm.ppf(1-float(tokens[-1]))
                bulk.add(mast)
        bulk.done()
    
    mhc_i_reference_alleles = ['HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:03', 'HLA-A*02:06', 'HLA-A*03:01', 'HLA-A*11:01', 'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*26:01', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*31:01', 'HLA-A*32:01', 'HLA-A*33:01', 'HLA-A*68:01', 'HLA-A*68:02', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*15:01', 'HLA-B*35:01', 'HLA-B*40:01', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*51:01', 'HLA-B*53:01', 'HLA-B*57:01', 'HLA-B*58:01'];

    mhc_ii_reference_alleles = ['HLA-DPA1*01:03/DPB1*02:01', 'HLA-DPA1*01:03/DPB1*04:01', 'HLA-DPA1*02:01/DPB1*01:01', 'HLA-DPA1*02:01/DPB1*05:01', 'HLA-DPA1*03:01/DPB1*04:02', 'HLA-DQA1*01:01/DQB1*05:01', 'HLA-DQA1*01:02/DQB1*06:02', 'HLA-DQA1*05:01/DQB1*02:01', 'HLA-DQA1*05:01/DQB1*03:01', 'HLA-DRB1*01:01', 'HLA-DRB1*03:01', 'HLA-DRB1*04:01', 'HLA-DRB1*04:05', 'HLA-DRB1*07:01', 'HLA-DRB1*08:02', 'HLA-DRB1*09:01', 'HLA-DRB1*11:01', 'HLA-DRB1*12:01', 'HLA-DRB1*13:02', 'HLA-DRB1*15:01', 'HLA-DRB3*01:01', 'HLA-DRB3*02:02', 'HLA-DRB4*01:01', 'HLA-DRB5*01:01'];

    groupMap = {}
    for row in TVaxignAlleleGroup.objects.filter(mhc_allele__in=mhc_i_reference_alleles+mhc_ii_reference_alleles):
        groupMap[row.c_allele_group_id] = row.mhc_allele
    
    mhc_i_epitope_restriction = {}
    mhc_ii_epitope_restriction = {}
    combine_epitope_restriction = {}
    for row in TVaxignMastResults.objects.filter(
        Q(c_allele_group_id__in=groupMap.keys()) & Q(c_sequence_id=seqID) & Q(c_hit_p_value__lte=0.05)
    ):
        epitope = sequence.c_sequence[row.c_hit_start-1:row.c_hit_end]
        if epitope == '':
            continue
        allele = groupMap[row.c_allele_group_id]
        if allele in mhc_i_reference_alleles:
            if epitope not in mhc_i_epitope_restriction.keys():
                mhc_i_epitope_restriction[epitope] = set([])
            mhc_i_epitope_restriction[epitope].add(allele)
        elif allele in mhc_ii_reference_alleles:
            if epitope not in mhc_ii_epitope_restriction.keys():
                mhc_ii_epitope_restriction[epitope] = set([])
            mhc_ii_epitope_restriction[epitope].add(allele)
        if epitope not in combine_epitope_restriction.keys():
            combine_epitope_restriction[epitope] = set([])
        combine_epitope_restriction[epitope].add(allele)
    
    if TVaxignPopulationCoverage.objects.filter(c_sequence_id=seqID).count() == 0:
        if not os.path.exists(os.path.join(settings.WORKSPACE_DIR, queryID, queryID+'.fasta')):
            if not os.path.exists(settings.VAXIGN2_TMP_DIR):
                os.mkdir(settings.VAXIGN2_TMP_DIR)
            if not os.path.exists(os.path.join(settings.VAXIGN2_TMP_DIR, queryID)):
                os.mkdir(os.path.join(settings.VAXIGN2_TMP_DIR, queryID))
            queryPath = os.path.join(settings.VAXIGN2_TMP_DIR, queryID)
        else:
            queryPath = os.path.join(settings.WORKSPACE_DIR, queryID)
            
        bulk = BulkCreateManager()
        
        # MHC-I
        output = []
        for epitope, alleles in mhc_i_epitope_restriction.items():
            output.append(str.format("{}\t{}", epitope, ','.join(alleles)))
        open(os.path.join(queryPath, seqID+'_mhc_i_popcov.input'), 'w').write('\n'.join(output))
        f = open(os.path.join(queryPath, seqID+'_mhc_i_popcov.output') , 'w' )
        subprocess.call( ['python2', os.path.join(settings.IEDB_PATH, 'population_coverage', 'calculate_population_coverage.py'), 
                          '-c', 'I', '-p', ','.join(list(settings.POPULATION_COVERAGE_ISO_CODE.keys())),
                          '-f', os.path.join(queryPath, seqID+'_mhc_i_popcov.input'),
                          ], stdout = f )
        f.close()
        for line in open(os.path.join(queryPath, seqID+'_mhc_i_popcov.output')).read().splitlines()[2:]:
            if line.startswith('average'):
                break
            tokens = line.split('\t')
            popcov = TVaxignPopulationCoverage(
                c_sequence_id=seqID,
                mhc_class='I',
                c_country=tokens[0],
                c_coverage=tokens[1].strip('%'),
            )
            bulk.add(popcov)
        
        # MHC-II
        output = []
        for epitope, alleles in mhc_ii_epitope_restriction.items():
            output.append(str.format("{}\t{}", epitope, ','.join(alleles)))
        open(os.path.join(queryPath, seqID+'_mhc_ii_popcov.input'), 'w').write('\n'.join(output))
        f = open(os.path.join(queryPath, seqID+'_mhc_ii_popcov.output') , 'w' )
        subprocess.call( ['python2', os.path.join(settings.IEDB_PATH, 'population_coverage', 'calculate_population_coverage.py'), 
                          '-c', 'II', '-p', ','.join(list(settings.POPULATION_COVERAGE_ISO_CODE.keys())),
                          '-f', os.path.join(queryPath, seqID+'_mhc_ii_popcov.input'),
                          ], stdout = f )
        f.close()
        for line in open(os.path.join(queryPath, seqID+'_mhc_ii_popcov.output')).read().splitlines()[2:]:
            if line.startswith('average'):
                break
            tokens = line.split('\t')
            popcov = TVaxignPopulationCoverage(
                c_sequence_id=seqID,
                mhc_class='II',
                c_country=tokens[0],
                c_coverage=tokens[1].strip('%'),
            )
            bulk.add(popcov)
        
        # Combined
        output = []
        for epitope, alleles in combine_epitope_restriction.items():
            output.append(str.format("{}\t{}", epitope, ','.join(alleles)))
        open(os.path.join(queryPath, seqID+'_combine_popcov.input'), 'w').write('\n'.join(output))
        f = open(os.path.join(queryPath, seqID+'_combine_popcov.output') , 'w' )
        subprocess.call( ['python2', os.path.join(settings.IEDB_PATH, 'population_coverage', 'calculate_population_coverage.py'), 
                          '-c', 'combined', '-p', ','.join(list(settings.POPULATION_COVERAGE_ISO_CODE.keys())),
                          '-f', os.path.join(queryPath, seqID+'_combine_popcov.input'),
                          ], stdout = f )
        f.close()
        for line in open(os.path.join(queryPath, seqID+'_combine_popcov.output')).read().splitlines()[2:]:
            if line.startswith('average'):
                break
            tokens = line.split('\t')
            popcov = TVaxignPopulationCoverage(
                c_sequence_id=seqID,
                mhc_class='I,II',
                c_country=tokens[0],
                c_coverage=tokens[1].strip('%'),
            )
            bulk.add(popcov)
        
        # Update MySQL
        bulk.done()
        
        # Clean up
        os.remove(os.path.join(queryPath, seqID+'_mhc_i_popcov.input'))
        os.remove(os.path.join(queryPath, seqID+'_mhc_i_popcov.output'))
        os.remove(os.path.join(queryPath, seqID+'_mhc_ii_popcov.input'))
        os.remove(os.path.join(queryPath, seqID+'_mhc_ii_popcov.output'))
        os.remove(os.path.join(queryPath, seqID+'_combine_popcov.input'))
        os.remove(os.path.join(queryPath, seqID+'_combine_popcov.output'))
    
    popcov = TVaxignPopulationCoverage.objects.filter(c_sequence_id=seqID)
    mhc_i_data = {
        'code':[],
        'country':[],
        'popcov':[],
    }
    mhc_ii_data = {
        'code':[],
        'country':[],
        'popcov':[],
    }
    combine_data = {
        'code':[],
        'country':[],
        'popcov':[],
    }
    for row in popcov:
        if settings.POPULATION_COVERAGE_ISO_CODE[row.c_country] not in country_code:
            continue
        if row.mhc_class == 'I':
            mhc_i_data['code'].append(settings.POPULATION_COVERAGE_ISO_CODE[row.c_country]),
            mhc_i_data['country'].append(row.c_country)
            mhc_i_data['popcov'].append(row.c_coverage)
        elif row.mhc_class == 'II':
            mhc_ii_data['code'].append(settings.POPULATION_COVERAGE_ISO_CODE[row.c_country]),
            mhc_ii_data['country'].append(row.c_country)
            mhc_ii_data['popcov'].append(row.c_coverage)
        elif row.mhc_class == 'I,II':
            combine_data['code'].append(settings.POPULATION_COVERAGE_ISO_CODE[row.c_country]),
            combine_data['country'].append(row.c_country)
            combine_data['popcov'].append(row.c_coverage)
    
    if mhc_class == 'I':
        mhc_i_data = pd.DataFrame(mhc_i_data)
        mhc_i_fig = go.Figure(data=go.Choropleth(
            locations = mhc_i_data['code'],
            z = mhc_i_data['popcov'],
            text = mhc_i_data['country'],
            colorscale = 'Portland',
            autocolorscale = False,
            reversescale = True,
            marker_line_color = 'darkgray',
            marker_line_width = 0.5,
            colorbar_ticksuffix = '%',
            colorbar_title = 'Population Coverage (%)',
            ),
            layout = {
                'width': 1000,
                'height': 600,
                'geo': {
                    'landcolor':'lightgray',
                    'countrycolor':'lightgray',
                    'showframe':False,
                    'showcoastlines':False,
                    'projection_type':'equirectangular',
                },
            },
        )
        context['mhc_i'] = offline.plot(mhc_i_fig, output_type='div')
    elif mhc_class == 'II':
        mhc_ii_data = pd.DataFrame(mhc_ii_data)
        mhc_ii_fig = go.Figure(data=go.Choropleth(
            locations = mhc_ii_data['code'],
            z = mhc_ii_data['popcov'],
            text = mhc_ii_data['country'],
            colorscale = 'Portland',
            autocolorscale = False,
            reversescale = True,
            marker_line_color = 'darkgray',
            marker_line_width = 0.5,
            colorbar_ticksuffix = '%',
            colorbar_title = 'Population Coverage (%)',
            ),
            layout = {
                'width': 1000,
                'height': 600,
                'geo': {
                    'landcolor':'lightgray',
                    'countrycolor':'lightgray',
                    'showframe':False,
                    'showcoastlines':False,
                    'projection_type':'equirectangular',
                },
            },
        )
        context['mhc_ii'] = offline.plot(mhc_ii_fig, output_type='div')
    else:
        combine_data = pd.DataFrame(combine_data)
        combine_fig = go.Figure(data=go.Choropleth(
            locations = combine_data['code'],
            z = combine_data['popcov'],
            text = combine_data['country'],
            colorscale = 'Portland',
            autocolorscale = False,
            reversescale = True,
            marker_line_color = 'darkgray',
            marker_line_width = 0.5,
            colorbar_ticksuffix = '%',
            colorbar_title = 'Population Coverage (%)',
            ),
            layout = {
                'width': 1000,
                'height': 600,
                'geo': {
                    'landcolor':'lightgray',
                    'countrycolor':'lightgray',
                    'showframe':False,
                    'showcoastlines':False,
                    'projection_type':'equirectangular',
                },
            },
        )
        context['combine'] = offline.plot(combine_fig, output_type='div')
    
    
    return render(request, 'queries/tabs/population_coverage.html', context)

def population_coverage(request, queryID, mhc_class, country_code=None):
    
    if country_code == None:
        country_code = settings.POPULATION_COVERAGE_ISO_CODE.values()
    else:
        country_code = country_code.split(',')
    
    context = {'query_id': queryID}
    
    mhc_i_reference_alleles = ['HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:03', 'HLA-A*02:06', 'HLA-A*03:01', 'HLA-A*11:01', 'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*26:01', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*31:01', 'HLA-A*32:01', 'HLA-A*33:01', 'HLA-A*68:01', 'HLA-A*68:02', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*15:01', 'HLA-B*35:01', 'HLA-B*40:01', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*51:01', 'HLA-B*53:01', 'HLA-B*57:01', 'HLA-B*58:01'];

    mhc_ii_reference_alleles = ['HLA-DPA1*01:03/DPB1*02:01', 'HLA-DPA1*01:03/DPB1*04:01', 'HLA-DPA1*02:01/DPB1*01:01', 'HLA-DPA1*02:01/DPB1*05:01', 'HLA-DPA1*03:01/DPB1*04:02', 'HLA-DQA1*01:01/DQB1*05:01', 'HLA-DQA1*01:02/DQB1*06:02', 'HLA-DQA1*05:01/DQB1*02:01', 'HLA-DQA1*05:01/DQB1*03:01', 'HLA-DRB1*01:01', 'HLA-DRB1*03:01', 'HLA-DRB1*04:01', 'HLA-DRB1*04:05', 'HLA-DRB1*07:01', 'HLA-DRB1*08:02', 'HLA-DRB1*09:01', 'HLA-DRB1*11:01', 'HLA-DRB1*12:01', 'HLA-DRB1*13:02', 'HLA-DRB1*15:01', 'HLA-DRB3*01:01', 'HLA-DRB3*02:02', 'HLA-DRB4*01:01', 'HLA-DRB5*01:01'];

    groupMap = {}
    for row in TVaxignAlleleGroup.objects.filter(mhc_allele__in=mhc_i_reference_alleles+mhc_ii_reference_alleles):
        groupMap[row.c_allele_group_id] = row.mhc_allele
    
    mhc_i_data = {
        'code':[],
        'country':[],
        'popcov':[],
    }
    mhc_ii_data = {
        'code':[],
        'country':[],
        'popcov':[],
    }
    combine_data = {
        'code':[],
        'country':[],
        'popcov':[],
    }
    for sequence in TVaxignAnalysisResults.objects.filter(c_query_id=queryID):
        seqID = str(sequence.c_sequence_id)
        mhc_i_epitope_restriction = {}
        mhc_ii_epitope_restriction = {}
        combine_epitope_restriction = {}
        for row in TVaxignMastResults.objects.filter(
            Q(c_allele_group_id__in=groupMap.keys()) & Q(c_sequence_id=seqID) & Q(c_hit_p_value__lte=0.05)
        ):
            epitope = sequence.c_sequence[row.c_hit_start-1:row.c_hit_end]
            if epitope == '':
                continue
            allele = groupMap[row.c_allele_group_id]
            if allele in mhc_i_reference_alleles:
                if epitope not in mhc_i_epitope_restriction.keys():
                    mhc_i_epitope_restriction[epitope] = set([])
                mhc_i_epitope_restriction[epitope].add(allele)
            elif allele in mhc_ii_reference_alleles:
                if epitope not in mhc_ii_epitope_restriction.keys():
                    mhc_ii_epitope_restriction[epitope] = set([])
                mhc_ii_epitope_restriction[epitope].add(allele)
            if epitope not in combine_epitope_restriction.keys():
                combine_epitope_restriction[epitope] = set([])
            combine_epitope_restriction[epitope].add(allele)
        
        if TVaxignPopulationCoverage.objects.filter(c_sequence_id=seqID).count() == 0:
            if not os.path.exists(os.path.join(settings.WORKSPACE_DIR, queryID, queryID+'.fasta')):
                if not os.path.exists(settings.VAXIGN2_TMP_DIR):
                    os.mkdir(settings.VAXIGN2_TMP_DIR)
                if not os.path.exists(os.path.join(settings.VAXIGN2_TMP_DIR, queryID)):
                    os.mkdir(os.path.join(settings.VAXIGN2_TMP_DIR, queryID))
                queryPath = os.path.join(settings.VAXIGN2_TMP_DIR, queryID)
            else:
                queryPath = os.path.join(settings.WORKSPACE_DIR, queryID)
                
            bulk = BulkCreateManager()
            
            # MHC-I
            output = []
            for epitope, alleles in mhc_i_epitope_restriction.items():
                output.append(str.format("{}\t{}", epitope, ','.join(alleles)))
            open(os.path.join(queryPath, seqID+'_mhc_i_popcov.input'), 'w').write('\n'.join(output))
            f = open(os.path.join(queryPath, seqID+'_mhc_i_popcov.output') , 'w' )
            subprocess.call( ['python2', os.path.join(settings.IEDB_PATH, 'population_coverage', 'calculate_population_coverage.py'), 
                              '-c', 'I', '-p', ','.join(list(settings.POPULATION_COVERAGE_ISO_CODE.keys())),
                              '-f', os.path.join(queryPath, seqID+'_mhc_i_popcov.input'),
                              ], stdout = f )
            f.close()
            for line in open(os.path.join(queryPath, seqID+'_mhc_i_popcov.output')).read().splitlines()[2:]:
                if line.startswith('average'):
                    break
                tokens = line.split('\t')
                popcov = TVaxignPopulationCoverage(
                    c_sequence_id=seqID,
                    mhc_class='I',
                    c_country=tokens[0],
                    c_coverage=tokens[1].strip('%'),
                )
                bulk.add(popcov)
            
            # MHC-II
            output = []
            for epitope, alleles in mhc_ii_epitope_restriction.items():
                output.append(str.format("{}\t{}", epitope, ','.join(alleles)))
            open(os.path.join(queryPath, seqID+'_mhc_ii_popcov.input'), 'w').write('\n'.join(output))
            f = open(os.path.join(queryPath, seqID+'_mhc_ii_popcov.output') , 'w' )
            subprocess.call( ['python2', os.path.join(settings.IEDB_PATH, 'population_coverage', 'calculate_population_coverage.py'), 
                              '-c', 'II', '-p', ','.join(list(settings.POPULATION_COVERAGE_ISO_CODE.keys())),
                              '-f', os.path.join(queryPath, seqID+'_mhc_ii_popcov.input'),
                              ], stdout = f )
            f.close()
            for line in open(os.path.join(queryPath, seqID+'_mhc_ii_popcov.output')).read().splitlines()[2:]:
                if line.startswith('average'):
                    break
                tokens = line.split('\t')
                popcov = TVaxignPopulationCoverage(
                    c_sequence_id=seqID,
                    mhc_class='II',
                    c_country=tokens[0],
                    c_coverage=tokens[1].strip('%'),
                )
                bulk.add(popcov)
            
            # Combined
            output = []
            for epitope, alleles in combine_epitope_restriction.items():
                output.append(str.format("{}\t{}", epitope, ','.join(alleles)))
            open(os.path.join(queryPath, seqID+'_combine_popcov.input'), 'w').write('\n'.join(output))
            f = open(os.path.join(queryPath, seqID+'_combine_popcov.output') , 'w' )
            subprocess.call( ['python2', os.path.join(settings.IEDB_PATH, 'population_coverage', 'calculate_population_coverage.py'), 
                              '-c', 'combined', '-p', ','.join(list(settings.POPULATION_COVERAGE_ISO_CODE.keys())),
                              '-f', os.path.join(queryPath, seqID+'_combine_popcov.input'),
                              ], stdout = f )
            f.close()
            for line in open(os.path.join(queryPath, seqID+'_combine_popcov.output')).read().splitlines()[2:]:
                if line.startswith('average'):
                    break
                tokens = line.split('\t')
                popcov = TVaxignPopulationCoverage(
                    c_sequence_id=seqID,
                    mhc_class='I,II',
                    c_country=tokens[0],
                    c_coverage=tokens[1].strip('%'),
                )
                bulk.add(popcov)
            
            # Update MySQL
            bulk.done()
            
            # Clean up
            os.remove(os.path.join(queryPath, seqID+'_mhc_i_popcov.input'))
            os.remove(os.path.join(queryPath, seqID+'_mhc_i_popcov.output'))
            os.remove(os.path.join(queryPath, seqID+'_mhc_ii_popcov.input'))
            os.remove(os.path.join(queryPath, seqID+'_mhc_ii_popcov.output'))
            os.remove(os.path.join(queryPath, seqID+'_combine_popcov.input'))
            os.remove(os.path.join(queryPath, seqID+'_combine_popcov.output'))
    

        popcov = TVaxignPopulationCoverage.objects.filter(c_sequence_id=seqID)
        for row in popcov:
            if settings.POPULATION_COVERAGE_ISO_CODE[row.c_country] not in country_code:
                continue
            if row.mhc_class == 'I':
                mhc_i_data['code'].append(settings.POPULATION_COVERAGE_ISO_CODE[row.c_country]),
                mhc_i_data['country'].append(row.c_country)
                mhc_i_data['popcov'].append(row.c_coverage)
            elif row.mhc_class == 'II':
                mhc_ii_data['code'].append(settings.POPULATION_COVERAGE_ISO_CODE[row.c_country]),
                mhc_ii_data['country'].append(row.c_country)
                mhc_ii_data['popcov'].append(row.c_coverage)
            elif row.mhc_class == 'I,II':
                combine_data['code'].append(settings.POPULATION_COVERAGE_ISO_CODE[row.c_country]),
                combine_data['country'].append(row.c_country)
                combine_data['popcov'].append(row.c_coverage)
    
    if mhc_class == 'I':
        mhc_i_data = pd.DataFrame(mhc_i_data)
        mhc_i_fig = go.Figure(data=go.Choropleth(
            locations = mhc_i_data['code'],
            z = mhc_i_data['popcov'],
            text = mhc_i_data['country'],
            colorscale = 'Portland',
            autocolorscale = False,
            reversescale = True,
            marker_line_color = 'darkgray',
            marker_line_width = 0.5,
            colorbar_ticksuffix = '%',
            colorbar_title = 'Population Coverage (%)',
            ),
            layout = {
                'width': 1200,
                'height': 650,
                'geo': {
                    'landcolor':'lightgray',
                    'countrycolor':'lightgray',
                    'showframe':False,
                    'showcoastlines':False,
                    'projection_type':'equirectangular',
                },
            },
        )
        context['mhc_i'] = offline.plot(mhc_i_fig, output_type='div')
    elif mhc_class == 'II':
        mhc_ii_data = pd.DataFrame(mhc_ii_data)
        mhc_ii_fig = go.Figure(data=go.Choropleth(
            locations = mhc_ii_data['code'],
            z = mhc_ii_data['popcov'],
            text = mhc_ii_data['country'],
            colorscale = 'Portland',
            autocolorscale = False,
            reversescale = True,
            marker_line_color = 'darkgray',
            marker_line_width = 0.5,
            colorbar_ticksuffix = '%',
            colorbar_title = 'Population Coverage (%)',
            ),
            layout = {
                'width': 1200,
                'height': 650,
                'geo': {
                    'landcolor':'lightgray',
                    'countrycolor':'lightgray',
                    'showframe':False,
                    'showcoastlines':False,
                    'projection_type':'equirectangular',
                },
            },
        )
        context['mhc_ii'] = offline.plot(mhc_ii_fig, output_type='div')
    else:
        combine_data = pd.DataFrame(combine_data)
        combine_fig = go.Figure(data=go.Choropleth(
            locations = combine_data['code'],
            z = combine_data['popcov'],
            text = combine_data['country'],
            colorscale = 'Portland',
            autocolorscale = False,
            reversescale = True,
            marker_line_color = 'darkgray',
            marker_line_width = 0.5,
            colorbar_ticksuffix = '%',
            colorbar_title = 'Population Coverage (%)',
            ),
            layout = {
                'width': 1200,
                'height': 650,
                'geo': {
                    'landcolor':'lightgray',
                    'countrycolor':'lightgray',
                    'showframe':False,
                    'showcoastlines':False,
                    'projection_type':'equirectangular',
                },
            },
        )
        context['combine'] = offline.plot(combine_fig, output_type='div')
    
    
    return render(request, 'queries/population_coverage.html', context)