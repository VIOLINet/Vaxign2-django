from django.shortcuts import render, redirect

# Create your views here.

import os
import re
import json
import urllib
import pymysql
import collections
import numpy as np
import plotly.graph_objects as go

from plotly import offline

from django.conf import settings
from django.contrib import messages

from django.db.models import Q
from vaxign.models import TVaxignQuery
from vaxign.models import TVaxignAnalysisResults
from vaxign.models import OrthomclResultsSeqs
from vaxign.models import TVaxignBlastResults
from vaxign.models import TVaxignMastResults
from vaxign.models import TVaxignAlleleGroup

from vaxign.forms.results import ResultsForm
from vaxign.forms.ortholog import OrthologForm

from vaxign.views.runs import *
from vaxign.views.runs import _multiprocess

import logging
from django.core.exceptions import ValidationError
from django.http.response import HttpResponse
logger = logging.getLogger('console')

def index(request, queryID):
    
    context = {'query_id': queryID}
    
    return render(request, 'queries/vaxitop.html', context)

def heatmap(request, queryID):
    
    context = {}
    
    alleles = request.POST['alleles']
    pval = 0.05
    
    groupMap = {}
    groupIDMap = {}
    for row in TVaxignAlleleGroup.objects.all():
        groupMap[re.sub("[^a-zA-Z0-9]", '', row.mhc_allele)] = row.mhc_allele
        groupIDMap[row.c_allele_group_id] = row.mhc_allele
    
    
    selected = set([])
    for select in alleles.split(','):
        if select == '':
            continue
        tokens = select.split('_')
        if tokens[2] == 'any':
            matched = TVaxignAlleleGroup.objects.filter(
                Q(mhc_species=tokens[0]) & Q(mhc_class=tokens[1])
            ).values_list('c_allele_group_id', flat=True)
        elif tokens[3] == 'any':
            matched = TVaxignAlleleGroup.objects.filter(
                Q(mhc_species=tokens[0]) & Q(mhc_class=tokens[1]) & Q(mhc_allele=groupMap[tokens[2]])
            ).values_list('c_allele_group_id', flat=True)
        else:
            matched = TVaxignAlleleGroup.objects.filter(
                Q(mhc_species=tokens[0]) & Q(mhc_class=tokens[1]) & Q(mhc_allele=groupMap[tokens[2]]) & Q(epitope_length=tokens[3])
            ).values_list('c_allele_group_id', flat=True)
        selected.update(matched)
    sequences = collections.OrderedDict()
    for row in TVaxignAnalysisResults.objects.filter(c_query_id=queryID).order_by('c_sequence_id'):
        if row.c_protein_accession is not None and row.c_protein_accession != '':
            sequences[row.c_sequence_id] = {
                'sequence': row.c_sequence,
                'protein': row.c_protein_accession,
            }
        elif row.c_sequence_acc is not None and row.c_sequence_acc != '':
            sequences[row.c_sequence_id] = {
                'sequence': row.c_sequence,
                'protein': row.c_sequence_acc,
            }
        else:
            sequences[row.c_sequence_id] = {
                'sequence': row.c_sequence,
                'protein': row.c_sequence_id,
            }
    
    masts = {}
    for sequence in sequences.keys():
        masts[sequence] = {}
        for allele in selected:
            masts[sequence][allele] = 0
    
    for row in TVaxignMastResults.objects.filter(
        Q(c_sequence_id__in=sequences.keys()) & Q(c_allele_group_id__in=selected) & Q(c_hit_p_value__lte=pval)
    ):
        masts[row.c_sequence_id][row.c_allele_group_id] += 1 
    
    data = []
    
    for sequence in sequences.keys():
        row = []
        for allele in selected:
            row.append(masts[sequence][allele])
        data.append(row)
    data = np.array(data)
    
    fig = go.Figure(
        data = go.Heatmap(
            z = data,
            x = [groupIDMap[i] for i in selected],
            y = [sequences[i]['protein'] for i in sequences.keys()],
            colorscale = 'Hot',
            reversescale = True,
        ),
        layout = {
            'autosize': False,
            'height': 35*len(sequences),
            'width': 15*len(selected),
        },
    )
    
    context['heatmap'] = offline.plot(fig, output_type='div')
    
    return render(request, 'queries/vaxitop_heatmap.html', context)

def protein(request, queryID, seqID):
    
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
    
    return render(request, 'queries/tabs/vaxitop.html', context)

def iedb_search(request, queryID, seqID):
    
    context = {}
    
    sequence = TVaxignAnalysisResults.objects.get(c_sequence_id=seqID).c_sequence
    context['sequence'] = sequence
    
    heDB = pymysql.connect(
        host = settings.MYSQL_VIOLIN_HOST_IP,
        user = settings.MYSQL_VIOLIN_USER_NAME,
        password = settings.MYSQL_VIOLIN_USER_PWD,
    )
    
    # Search IEDB Linear T Cell Epitopes
    t_cell_epitopes = {}
    try:
        with heDB.cursor() as cursor:
            cursor.execute(str.format("""
SELECT iedb_public.epitope.epitope_id, iedb_public.epitope.linear_peptide_seq
FROM iedb_public.epitope
LEFT JOIN iedb_public.epitope_object ON iedb_public.epitope.epitope_id = iedb_public.epitope_object.epitope_id
LEFT JOIN iedb_public.curated_epitope ON iedb_public.epitope_object.object_id = iedb_public.curated_epitope.e_object_id
LEFT JOIN iedb_public.tcell ON iedb_public.curated_epitope.curated_epitope_id = iedb_public.tcell.curated_epitope_id
WHERE 
    '{}' LIKE CONCAT( '%', iedb_public.epitope.linear_peptide_seq, '%' )
AND
    CHAR_LENGTH( iedb_public.epitope.linear_peptide_seq ) >= 8
AND
    iedb_public.tcell.as_char_value != 'Negative'
GROUP BY iedb_public.epitope.epitope_id
            """, sequence))
            for row in cursor.fetchall():
                match = re.search(row[1], sequence)
                tmp = {
                    'epitope_id': row[0],
                    'linear_peptide_seq': row[1],
                    'alleles': {},
                    'start_pos': match.start,
                    'end_pos': match.end,
                }
                t_cell_epitopes[row[0]] = tmp
    finally:
        cursor.close()
    for epitope_id in t_cell_epitopes.keys():
        try:
            with heDB.cursor() as cursor:
                cursor.execute(str.format("""
SELECT iedb_public.mhc_allele_restriction.mhc_allele_restriction_id, iedb_public.mhc_allele_restriction.displayed_restriction FROM iedb_public.epitope
LEFT JOIN iedb_public.epitope_object ON iedb_public.epitope.epitope_id = iedb_public.epitope_object.epitope_id
LEFT JOIN iedb_public.curated_epitope ON iedb_public.epitope_object.object_id = iedb_public.curated_epitope.e_object_id
LEFT JOIN iedb_public.tcell ON iedb_public.curated_epitope.curated_epitope_id = iedb_public.tcell.curated_epitope_id
LEFT JOIN iedb_public.mhc_allele_restriction ON iedb_public.mhc_allele_restriction.mhc_allele_restriction_id = iedb_public.tcell.mhc_allele_restriction_id
WHERE
    iedb_public.epitope.epitope_id = {}
    AND
    iedb_public.mhc_allele_restriction.mhc_allele_restriction_id IS NOT NULL
    ;
                """, epitope_id))
                for row in cursor.fetchall():
                    t_cell_epitopes[epitope_id]['alleles'][row[0]] = row[1]
        finally:
            cursor.close()
    context['t_cell_epitopes'] = t_cell_epitopes.values()
    
    # Search IEDB Linear B Cell Epitopes
    b_cell_epitopes = {}
    try:
        with heDB.cursor() as cursor:
            cursor.execute(str.format("""
SELECT iedb_public.epitope.epitope_id, iedb_public.epitope.linear_peptide_seq
FROM iedb_public.epitope
LEFT JOIN iedb_public.epitope_object ON iedb_public.epitope.epitope_id = iedb_public.epitope_object.epitope_id
LEFT JOIN iedb_public.curated_epitope ON iedb_public.epitope_object.object_id = iedb_public.curated_epitope.e_object_id
LEFT JOIN iedb_public.bcell ON iedb_public.curated_epitope.curated_epitope_id = iedb_public.bcell.curated_epitope_id
WHERE 
    '{}' LIKE CONCAT( '%', iedb_public.epitope.linear_peptide_seq, '%' )
AND
    CHAR_LENGTH( iedb_public.epitope.linear_peptide_seq ) >= 5
AND
    iedb_public.bcell.as_char_value != 'Negative'
GROUP BY iedb_public.epitope.epitope_id
            """, sequence))
            for row in cursor.fetchall():
                match = re.search(row[1], sequence)
                tmp = {
                    'epitope_id': row[0],
                    'linear_peptide_seq': row[1],
                    'alleles': {},
                    'start_pos': match.start,
                    'end_pos': match.end,
                }
                b_cell_epitopes[row[0]] = tmp
    finally:
        cursor.close()
    context['b_cell_epitopes'] = b_cell_epitopes.values()
    
    return render(request, 'queries/tabs/iedb_epitope.html', context)