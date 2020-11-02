from django.shortcuts import render, redirect

# Create your views here.

import re
import json
import urllib
import pymysql
import collections

from django.conf import settings
from django.contrib import messages

from django.db.models import Q
from vaxign.models import TVaxignQuery
from vaxign.models import TVaxignAnalysisResults
from vaxign.models import OrthomclResultsSeqs
from vaxign.models import TVaxignBlastResults
from vaxign.models import TVaxignMastResults

from vaxign.forms.results import ResultsForm
from vaxign.forms.ortholog import OrthologForm

import logging
logger = logging.getLogger('console')


def results(request, queryID):
    
    if 'results_form_data' in request.session:
        if request.session['results_form_data']['query_id'] != queryID:
            request.session['results_form_data'] = {
                'query_id':queryID,
                'have_orthologs':[],
                'have_no_orthologs':[],
                'localization': ['Any'],
                'tmhmm_PredHel_check': False,
                'tmhmm_PredHel_value': 1,
                'tmhmm_PredHel_opt': '<=',
                'spaan_score_check': False,
                'spaan_score_value': 0.51,
                'spaan_score_opt': '>=',
                'human_alignment': '',
                'mouse_alignment': '',
                'pig_alignment': '',
            }
    else:
        request.session['results_form_data'] = {
                'query_id':queryID,
                'have_orthologs':[],
                'have_no_orthologs':[],
                'localization': ['Any'],
                'tmhmm_PredHel_check': False,
                'tmhmm_PredHel_value': 1,
                'tmhmm_PredHel_opt': '<=',
                'spaan_score_check': False,
                'spaan_score_value': 0.51,
                'spaan_score_opt': '>=',
                'human_alignment': '',
                'mouse_alignment': '',
                'pig_alignment': '',
        }
    
    context = {
        'query_id': queryID,
    }
    
    if bool(request.POST):
        for key in request.POST.keys():
            if key in request.session['results_form_data'].keys():
                if key in ['localization', 'have_orthologs', 'have_no_orthologs']:
                    request.session['results_form_data'][key] = request.POST.getlist(key)
                else:
                    request.session['results_form_data'][key] = request.POST[key]
    
    form = ResultsForm(request.session['results_form_data'])
    if form.is_valid():
        formData = form.cleaned_data
    else:
        messages.info(request, form.errors)
        return redirect(request.META.get('HTTP_REFERER', '/'))
    context['form'] = form
    
    query = TVaxignQuery.objects.get(
        Q(c_query_id=queryID)
    )
    context['query_detail'] = query
    
    query.c_included_analyses = query.c_included_analyses.split( ',' )
    
    if query.c_finish_time != '' and query.c_status == 'All selected steps finished':
        
        logger.debug('Query analysis finished. Gathering results...')
        
        # Initiate Q expression lists for model query
        Qs = Q(c_query_id=queryID)
        
        # Orthologs
        if 'have_orthologs' in formData.keys() and len(formData['have_orthologs']) != 0:
            logger.debug('Option [have orthologs] selected.')
            tmpOrtholog = formData['have_orthologs']
            if queryID in tmpOrtholog:
                tmpOrtholog.remove(queryID)
            formData['have_num_orthologs'] = len(tmpOrtholog)-1
            
            orthologResults = []
            for orthologObject in OrthomclResultsSeqs.objects.raw(str.format(
"""
SELECT result_id, count(distinct genome) as no_genomes FROM violin.orthomcl_results_seqs where genome in (
    '{}'
) group by result_id having no_genomes>={}
""", "','".join(tmpOrtholog), len(tmpOrtholog)
            )):
                orthologResults.append(orthologObject.result_id)
            
            orthologSeqs = OrthomclResultsSeqs.objects.filter(
                Q(result_id__in=orthologResults) & Q(genome=queryID)
            ).order_by('c_sequence_id').values_list('c_sequence_id',flat=True).distinct()
            
            Qs &= Q(c_sequence_id__in=orthologSeqs)
            
        # No orthologs
        if 'have_no_orthologs' in formData.keys() and len(formData['have_no_orthologs']) != 0:
            logger.debug('Option [exclude orthologs] selected.')
            tmpNoOrtholog = formData['have_no_orthologs']
            noOrthologResults = []
            for orthologObject in OrthomclResultsSeqs.objects.raw(str.format(
"""
SELECT result_id, count(distinct genome) as no_genomes FROM violin.orthomcl_results_seqs where genome in (
    '{}'
) group by result_id having no_genomes>={}
""", "','".join(tmpNoOrtholog), len(tmpNoOrtholog)
            )):
                noOrthologResults.append(orthologObject.result_id)
            
            noOrthologSeqs = OrthomclResultsSeqs.objects.filter(
                Q(result_id__in=noOrthologResults) & Q(genome=queryID)
            ).order_by('c_sequence_id').values_list('c_sequence_id',flat=True).distinct()
            
            Qs &= ~Q(c_sequence_id__in=noOrthologSeqs)
        
        # Localization
        if 'Any' not in formData['localization']:
            logger.debug('Option [localization] selected.')
            localizationQs = Q()
            for localization in formData['localization']:
                localizationQs |= Q(c_final_localization__icontains=localization)
            Qs &= (localizationQs)
        
        # Adhesin probability
        if formData['spaan_score_check']:
            logger.debug('Option [adhesin probability] selected.')
            if formData['spaan_score_opt'] == '>':
                Qs &= Q(c_spaan__gt=formData['spaan_score_value'])
            elif formData['spaan_score_opt'] == '>=':
                Qs &= Q(c_spaan__gte=formData['spaan_score_value'])
            elif formData['spaan_score_opt'] == '<':
                Qs &= Q(c_spaan__lt=formData['spaan_score_value'])
            elif formData['spaan_score_opt'] == '<=':
                Qs &= Q(c_spaan__lte=formData['spaan_score_value'])
        
        # Transmembrane helix
        if formData['tmhmm_PredHel_check']:
            logger.debug('Option [transmembrane helix] selected.')
            if formData['tmhmm_PredHel_opt'] == '>':
                Qs &= Q(c_tmhmm_predhel__gt=formData['tmhmm_PredHel_value'])
            elif formData['tmhmm_PredHel_opt'] == '>=':
                Qs &= Q(c_tmhmm_predhel__gte=formData['tmhmm_PredHel_value'])
            elif formData['tmhmm_PredHel_opt'] == '<':
                Qs &= Q(c_tmhmm_predhel__lt=formData['tmhmm_PredHel_value'])
            elif formData['tmhmm_PredHel_opt'] == '<=':
                Qs &= Q(c_tmhmm_predhel__lte=formData['tmhmm_PredHel_value'])
                
        # Human similarity
        logger.debug('Gathering human blast results...')
        humanSeqs = TVaxignBlastResults.objects.filter(
            Q(c_query_id=queryID) & Q(c_target_species='Human')
        ).distinct().values_list('c_sequence_id', flat=True)
        if formData['human_alignment'] == 'y':
            logger.debug('Option [similarity to human] selected.')
            Qs &= Q(c_sequence_id__in=humanSeqs)
        elif formData['human_alignment'] == 'n':
            Qs &= ~Q(c_sequence_id__in=humanSeqs)
        
        # Mouse similarity
        logger.debug('Gathering mouse blast results...')
        mouseSeqs = TVaxignBlastResults.objects.filter(
            Q(c_query_id=queryID) & Q(c_target_species='Mouse')
        ).distinct().values_list('c_sequence_id', flat=True)
        if formData['mouse_alignment'] == 'y':
            logger.debug('Option [similarity to mouse] selected.')
            Qs &= Q(c_sequence_id__in=mouseSeqs)
        elif formData['mouse_alignment'] == 'n':
            Qs &= ~Q(c_sequence_id__in=mouseSeqs)
        
        # Pig similarity
        logger.debug('Gathering pig blast results...')
        pigSeqs = TVaxignBlastResults.objects.filter(
            Q(c_query_id=queryID) & Q(c_target_species='Pig')
        ).distinct().values_list('c_sequence_id', flat=True)
        if formData['pig_alignment'] == 'y':
            logger.debug('Option [similarity to pig] selected.')
            Qs &= Q(c_sequence_id__in=pigSeqs)
        elif formData['pig_alignment'] == 'n':
            Qs &= ~Q(c_sequence_id__in=pigSeqs)
        
        # Keyword
        if formData['keywords'] != '':
            logger.debug('Option [keyword] selected.')
            if formData['keywords_type'] == 'c_gene_symbol':
                Qs &= Q(c_gene_symbol=formData['keywords'])
            elif formData['keywords_type'] == 'c_note':
                Qs &= Q(c_note__icontains=formData['keywords'])
                
        # Sequence IDs
        if formData['query_seq_ids'] != '':
            logger.debug('Option [sequence ids] selected.')
            seqIds = re.split(r'[\s,;:]+', formData['query_seq_ids'])
            if formData['query_seq_id_type'] == 'protein_gi':
                Qs &= Q(c_protein_gi__in=seqIds)
            elif formData['query_seq_id_type'] == 'protein_accession':
                Qs &= Q(c_protein_accession__in=seqIds)
            elif formData['query_seq_id_type'] == 'gene_id':
                Qs &= Q(c_gene_id__in=seqIds)
            elif formData['query_seq_id_type'] == 'locus_tag':
                Qs &= Q(c_locus_tag__in=seqIds)
        
        # Filter Q list and bind model result
        logger.debug('Query options processed. Gathering data...')
        results = TVaxignAnalysisResults.objects.filter(Qs).only(
            'c_query_id',
            'c_sequence_id'
            'c_protein_accession',
            'c_sequence_acc',
            'c_note',
            'c_gene_id',
            'c_gene_symbol',
            'c_locus_tag',
            'c_protein_length',
            'c_final_localization',
            'c_final_score',
            'c_spaan',
            'c_tmhmm_predhel',
            'c_vaxign_ml_score',
        ).order_by('-c_vaxign_ml_score', 'c_sequence_id').values()
        reactomeGenes = []
        for (i,row) in enumerate(results):
            if row['c_sequence_id'] in humanSeqs:
                results[i]['human_alignment'] = 'Yes'
            else:
                results[i]['human_alignment'] = None
            if row['c_sequence_id'] in mouseSeqs:
                results[i]['mouse_alignment'] = 'Yes'
            else:
                results[i]['mouse_alignment'] = None
            if row['c_sequence_id'] in pigSeqs:
                results[i]['pig_alignment'] = 'Yes'
            else:
                results[i]['pig_alignment'] = None
            if row['c_gene_symbol'] != None and row['c_gene_symbol'] != '':
                reactomeGenes.append(row['c_gene_symbol'])
            if row['c_spaan'] != None and row['c_spaan'] < 0.001:
                results[i]['c_spaan'] = 0
        context['results'] = results
        
        # Bind orthologs
        orthologs = collections.OrderedDict()
        for ortholog in TVaxignQuery.objects.filter(
            Q(c_ortholog_computed=1) & ~Q(c_query_id=queryID) & Q(c_species_short=formData['group_name'])
        ).order_by('c_species_name'): 
            orthologs[ortholog.c_query_id] = ortholog.c_species_name
        context['orthologs'] = orthologs
        
        return render(request, 'queries/results.html', context)
    
    else:
        return render(request, 'queries/refresh.html', context)

def ortholog(request, queryID):
    
    context = {
        'query_id': queryID,
    }
    
    query = TVaxignQuery.objects.get(
        Q(c_query_id=queryID)
    )
    query.c_included_analyses = query.c_included_analyses.split( ',' )
    context['query_detail'] = query
    
    # Bind all possible orthologs
    logger.debug("Gathering all orthologs available...")
    allOrthologs = collections.OrderedDict()
    for ortholog in TVaxignQuery.objects.filter(
        Q(c_ortholog_computed=1) & ~Q(c_query_id=queryID) & Q(c_species_short=query.c_species_short)
    ).order_by('c_species_name'): 
        allOrthologs[ortholog.c_query_id] = ortholog.c_species_name
    context['all_orthologs'] = allOrthologs
    
    if query.c_finish_time != '' and query.c_status == 'All selected steps finished':
        
        logger.debug('Query analysis finished. Gathering results...')
        
        if query.c_ortholog_computed != 1:
            logger.debug("Ortholog analysis not available. Redirecting...")
            messages.info(request, 'The ortholog result for the selected genome group is not available!')
            return redirect(request.META.get('HTTP_REFERER', '/'))
        
        if request.method == 'POST':
            form = OrthologForm(request.POST)
            
            if form.is_valid():
                formData = form.cleaned_data
                
                # Bind selected orthologs
                matchOrthologs = collections.OrderedDict()
                for ortholog in allOrthologs:
                    if ortholog in formData['have_orthologs']:
                        matchOrthologs[ortholog] = allOrthologs[ortholog]
                context['match_orthologs'] = matchOrthologs
                
                # Get a list of proteins and details in the query genome
                logger.debug('Gathering query protein detail...')
                querySeqs = TVaxignAnalysisResults.objects.filter(c_query_id=queryID).only(
                    'c_protein_accession',
                    'c_sequence_acc',
                    'c_note',
                    'c_gene_id',
                    'c_gene_symbol',
                ).order_by('c_sequence_id')
                
                # Map each query genome's protein to the corresponding ortholog's protein
                logger.debug("Gathering mapping from query protein to ortholog protein...")
                tmpResult = OrthomclResultsSeqs.objects.filter(
                    Q(c_sequence_id__in=querySeqs.values_list('c_sequence_id', flat=True)) & Q(genome=queryID)
                )
                mapQueryToTmp = {}
                for row in tmpResult:
                    mapQueryToTmp[row.c_sequence_id] = row.result_id
                
                orthologResult = OrthomclResultsSeqs.objects.filter(
                    Q(result_id__in=tmpResult) & Q(genome__in=formData['have_orthologs'])
                )
                mapTmpToOrtholog = {}
                for row in orthologResult:
                    if row.result_id not in mapTmpToOrtholog.keys():
                        mapTmpToOrtholog[row.result_id] = {}
                    mapTmpToOrtholog[row.result_id][row.genome] = row.c_sequence_id
                
                logger.debug('Gathering ortholog protein detail...')
                orthologDetail = {}
                for row in TVaxignAnalysisResults.objects.filter(c_sequence_id__in=orthologResult.values_list('c_sequence_id',flat=True)).only(
                    'c_protein_accession',
                    'c_sequence_acc',
                    'c_note',
                    'c_gene_id',
                    'c_gene_symbol',
                ):
                    orthologDetail[row.c_sequence_id] = row
                
                # Bind ortholog results
                results = collections.OrderedDict()
                for row in querySeqs:
                    results[row.c_sequence_id] = {
                        'c_protein_accession': row.c_protein_accession,
                        'c_sequence_acc': row.c_sequence_acc,
                        'c_note': row.c_note,
                        'c_gene_id': row.c_gene_id,
                        'c_gene_symbol': row.c_gene_symbol,
                        'orthologs':[]
                    }
                    for ortholog in matchOrthologs:
                        try:
                            results[row.c_sequence_id]['orthologs'].append(orthologDetail[mapTmpToOrtholog[mapQueryToTmp[row.c_sequence_id]][ortholog]])
                        except:
                            results[row.c_sequence_id]['orthologs'].append(None)
                context['results'] = results
                    
        else:
            formData = {
                'query_id':queryID,
                'group_name':query.c_species_short,
                'have_orthologs':[],
            }
            form = OrthologForm(formData)
        context['form'] = form
        
        return render(request, 'queries/ortholog.html', context)
            
    else:
        return render(request, 'queries/refresh.html', context)

def protein(request, queryID, seqID):
    
    context = {
        'query_id': queryID,
        'sequence_id': seqID,
    }
    
    logger.debug('Gathering query protein sequence and detail...')
    query = TVaxignQuery.objects.get(c_query_id=queryID)
    context['query_detail'] = query 
    
    results = TVaxignAnalysisResults.objects.get(
        Q(c_query_id=queryID) & Q(c_sequence_id=seqID)
    )
    results = results.__dict__
    if results['c_spaan'] != None and results['c_spaan'] < 0.001:
                    results['c_spaan'] = 0
    
    # Human similarity
    logger.debug('Gathering human blast results...')
    humanSeqs = TVaxignBlastResults.objects.filter(
        Q(c_sequence_id=seqID) & Q(c_query_id=queryID) & Q(c_target_species='Human')
    ).distinct().values_list('c_target_id', flat=True)
    results['human_alignment'] = list(humanSeqs)
    
    # Mouse similarity
    logger.debug('Gathering mouse blast results...')
    mouseSeqs = TVaxignBlastResults.objects.filter(
        Q(c_sequence_id=seqID) & Q(c_query_id=queryID) & Q(c_target_species='Mouse')
    ).distinct().values_list('c_target_id', flat=True)
    results['mouse_alignment'] = list(mouseSeqs)
    
    # Pig similarity
    logger.debug('Gathering pig blast results...')
    pigSeqs = TVaxignBlastResults.objects.filter(
        Q(c_sequence_id=seqID) & Q(c_query_id=queryID) & Q(c_target_species='Pig')
    ).distinct().values_list('c_target_id', flat=True)
    results['pig_alignment'] = list(pigSeqs)
    
    
    context['results'] = results
    
    return render(request, 'queries/sequence.html', context)

def blast(request, queryID, seqID, organism):
    
    context = {
        'query_id': queryID,
        'sequence_id': seqID,
        'organism': organism,
    }
    
    logger.debug('Gathering query protein sequence and detail...')
    query = TVaxignAnalysisResults.objects.get(
        Q(c_query_id=queryID) & Q(c_sequence_id=seqID)
    )
    if query.c_protein_accession != None and query.c_protein_accession != '':
        context['protein_accession'] = query.c_protein_accession
    else:
        context['protein_accession'] = query.c_sequence_acc
    context['protein_note'] = query.c_note
    
    if organism == 'human':
        blasts = TVaxignBlastResults.objects.filter(
            Q(c_sequence_id=seqID) & Q(c_query_id=queryID) & Q(c_target_species='Human')
        ).order_by('c_hit_from')
    elif organism == 'mouse':
        blasts = TVaxignBlastResults.objects.filter(
            Q(c_sequence_id=seqID) & Q(c_query_id=queryID) & Q(c_target_species='Mouse')
        ).order_by('c_hit_from')
    elif organism == 'pig':
        blasts = TVaxignBlastResults.objects.filter(
            Q(c_sequence_id=seqID) & Q(c_query_id=queryID) & Q(c_target_species='Pig')
        ).order_by('c_hit_from')
    else:
        blasts = []
    results = collections.OrderedDict()
    for row in blasts:
        if row.c_target_id not in results.keys():
            results[row.c_target_id] = []
        results[row.c_target_id].append({
            'c_hit_from':"%4s" % str(row.c_hit_from),
            'c_hseq':row.c_hseq,
            'c_hit_to':str(row.c_hit_to),
            'c_midline':row.c_midline,
            'c_query_from':"%4s" % str(row.c_query_from),
            'c_qseq':row.c_qseq,
            'c_query_to':str(row.c_query_to),
            
            })
    context['results'] = results
    
    return render(request, 'queries/blast.html', context)
    