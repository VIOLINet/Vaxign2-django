from django.shortcuts import render, redirect

# Create your views here.

import os
import re
import subprocess
from Bio.KEGG import REST
from goatools.base import get_godag

from django.conf import settings

from django.db.models import Q
from vaxign.models import TVaxignAnalysisResults
from vaxign.models import TVaxignEggnogFunctions

import logging
logger = logging.getLogger('console')

def protein_eggnog_function(request, queryID, seqID):
    
    context = {'query_id': queryID, 'sequence_id': seqID}
    
    try:
        sequence = TVaxignAnalysisResults.objects.get(c_sequence_id=seqID)
        context['sequence'] = sequence.c_sequence
    except:
        return HttpResponse("No eggNOG Prediction available")
    
    if TVaxignEggnogFunctions.objects.filter(c_sequence_id=seqID).count() == 0:
        if not os.path.exists(os.path.join(settings.WORKSPACE_DIR, queryID, queryID+'.fasta')):
            if not os.path.exists(settings.VAXIGN2_TMP_DIR):
                os.mkdir(settings.VAXIGN2_TMP_DIR)
            if not os.path.exists(os.path.join(settings.VAXIGN2_TMP_DIR, queryID)):
                os.mkdir(os.path.join(settings.VAXIGN2_TMP_DIR, queryID))
                
            queryPath = os.path.join(settings.VAXIGN2_TMP_DIR, queryID)
            open(os.path.join(queryPath, queryID+'.fasta'), 'w').write(str.format("""
>{}
{}
            """, sequence.c_sequence_id, sequence.c_sequence))
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
    eggnog = TVaxignEggnogFunctions.objects.get(c_sequence_id=seqID)
    
    godag = get_godag('go-basic.obo')
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
    
    return render(request, 'queries/tabs/eggnog_function.html', context)

