from django.shortcuts import render, redirect

# Create your views here.

import os
import re
import sys
import subprocess
import multiprocessing
from pathlib import Path
from scipy.stats import norm
from collections import defaultdict

from django.conf import settings
from django.apps import apps

from django.db.models import Q
from django.db.models.functions import Now
from vaxign.models import TVaxignQuery
from vaxign.models import TVaxignAlleleGroup
from vaxign.models import TVaxignMastResults

from vaxign.views.email import email_submission

import logging
logger = logging.getLogger('console')

def index(request, queryID):
    
    logger.debug(str.format("Start running query: {}", queryID))
    
    context = {
        'query_id': queryID,
    }
    
    formData = request.session['runs_form_data']
    
    # Process analyses to be included
    included_analysis = formData['basic_analysis'].copy()
    if formData['vaxitop_analysis'] == 'y':
        included_analysis += ['mast_I','mast_II']
    if formData['vaxignml_analysis'] == 'y':
        included_analysis.append('vaxign_ml')
    
    logger.debug(str.format('Processed query {} parameters. Start running query...', queryID))
    
    # Clean up sequence special characters
    tmpSequence = []
    for line in formData['sequence'].split('\n'):
        if line.startswith( '>' ):
            tmpSequence.append(line) 
            continue
        line = line.upper()
        line = re.sub('[-\s]', '', line)
        line = re.sub('[^A-NP-TV-Z]', 'X', line)
        tmpSequence.append(line)
    formData['sequence'] = '\n'.join(tmpSequence)
    logger.debug('Cleaned up FASTA sequences.')
    
    os.mkdir(os.path.join(settings.WORKSPACE_DIR, queryID))
    logger.debug('Created workspace folder.')
    
    seqFile = str.format("{}.raw.fasta", queryID)
    open(os.path.join(settings.WORKSPACE_DIR, queryID, seqFile), 'w').write(formData['sequence'])
    logger.debug('Saved raw sequence in workspace.')
    
    if 'c_user_name' in request.session:
        user = request.session['c_user_name']
    else:
        user = ''
    query = TVaxignQuery(
        c_query_id = queryID,
        c_submit_time = Now(), 
        c_status = 'Loading Sequences',
        c_gram = formData['bacteria_strain'],
        c_species_name = formData['genome_name'],
        c_included_analyses = ','.join(included_analysis),
        c_user_name = user,
        c_is_public = '0',
        c_ortholog_computed = '0',
        c_species_short = formData['group_name'],
        c_note = formData['note'],
        c_projectfk = formData['c_projectfk'],
        c_has_z_value = 'n',
        c_organism = formData['organism'],
    )
    query.save()
    logger.debug('Inserted query record to database.')
    
    # Call vaxign-php to run the PHP version of Vaxign
    logger.debug('Started running vaxign-php in the background...')
    
    # Run load_sequence.php
    subprocess.Popen(['php',os.path.join(settings.VAXIGN_PHP_DIR, 'load_sequences.php'), queryID], cwd=settings.VAXIGN_PHP_DIR)
    
    # It includes the analyses of psortb, spaan, tmhmm, blasth, blastm, blasth, mast_I, and mast_II
    subprocess.Popen(['php',os.path.join(settings.VAXIGN_PHP_DIR, 'bg_run_query.php'), queryID], cwd=settings.VAXIGN_PHP_DIR)
    
    # Send email
    if formData['email'] != '':
        organism = ''
        if formData['organism'] is None or formData['organism'] == '':
            if formData['bacteria_strain'] == '--positive':
                organism = 'Gram Positive Bacteria'
            elif formData['bacteria_strain'] == '--negative':
                organism = 'Gram Negative Bacteria'
        else:
            if formData['organism'] == 'bacteria':
                if formData['bacteria_strain'] == '--positive':
                    organism = 'Gram Positive Bacteria'
                elif formData['bacteria_strain'] == '--negative':
                    organism = 'Gram Negative Bacteria'
            else:
                organism = formData['organism']
        analyses = []
        for analysis in included_analysis:
            if analysis == 'psort':
                analyses.append('Subcellular Localization')
            if analysis == 'spaan':
                analyses.append('Adhesion Probability')
            if analysis == 'tmhmm':
                analyses.append('Transmembrane Helices')
            if analysis == 'blasth':
                analyses.append('Similarity to Human Proteins')
            if analysis == 'blastm':
                analyses.append('Similarity to Mouse Proteins')
            if analysis == 'blastp':
                analyses.append('Similarity to Pig Proteins')
            if analysis == 'vaxign_ml':
                analyses.append('Vaxign-ML')
            if analysis == 'mast_I' or analysis == 'mast_II':
                analyses.append('Vaxitop')
        
        email_body = str.format("""
==================================================
Submission Note
==================================================
{}

==================================================
Organsim
==================================================
{}

==================================================
Included Analyes
==================================================
{}
        """, formData['note'], organism, analyses)
        
        email_submission(queryID, formData['email'], email_body)
    
    if formData['c_projectfk'] == 0:
        return redirect(str.format('/vaxign2/query/{}/results', queryID))
    else:
        return redirect(str.format('/vaxign2/project/{}', formData['c_projectfk']))
    

def vaxitop(request, queryID):
    
    logger.debug(str.format("Start running Vaxitop-only query: {}", queryID))
    
    context = {
        'query_id': queryID,
    }
    
    sequence = request.session['sequence']
    
    # Clean up sequence special characters
    tmpSequence = []
    for line in sequence.split('\n'):
        if line.startswith( '>' ):
            tmpSequence.append(line) 
            continue
        line = line.upper()
        line = re.sub('[-\s]', '', line)
        line = re.sub('[^A-NP-TV-Z]', 'X', line)
        tmpSequence.append(line)
    sequence = '\n'.join(tmpSequence)
    logger.debug('Cleaned up FASTA sequences.')
    
    os.mkdir(os.path.join(settings.WORKSPACE_DIR, queryID))
    logger.debug('Created workspace folder.')
    
    seqFile = str.format("{}.raw.fasta", queryID)
    open(os.path.join(settings.WORKSPACE_DIR, queryID, seqFile), 'w').write(sequence)
    logger.debug('Saved raw sequence in workspace.')
    
    query = TVaxignQuery(
        c_query_id = queryID,
        c_submit_time = Now(), 
        c_status = 'Loading Sequences',
        c_gram = '',
        c_species_name = '',
        c_included_analyses = ','.join(['mast_I','mast_II']),
        c_user_name = '',
        c_is_public = '0',
        c_ortholog_computed = '0',
        c_species_short = '',
        c_note = '',
        c_projectfk = 0,
        c_has_z_value = 'n',
        c_organism = '',
    )
    query.save()
    logger.debug('Inserted query record to database.')
    
    # Run load_sequence.php
    subprocess.Popen(['php',os.path.join(settings.VAXIGN_PHP_DIR, 'load_sequences.php'), queryID], cwd=settings.VAXIGN_PHP_DIR)
    
    alleleGroups = TVaxignAlleleGroup.objects.all().values_list('c_allele_group_id', flat=True)
    cmds = []
    for group in alleleGroups:
        cmd = str.format("{} {} {} -mt 0.1 -hit_list > {}",
                            os.path.join(settings.VAXITOP_PATH, 'lib', 'meme', 'bin', 'mast'),
                            os.path.join(settings.VAXITOP_PATH, 'pssm', str(group)+'.xml'),
                            os.path.join(settings.WORKSPACE_DIR, queryID, queryID+'.fasta'),
                            os.path.join(settings.WORKSPACE_DIR, queryID, str(group)+'.matching.txt'),
                        )
        cmds.append(cmd)
    _multiprocess(cmds, 10)
    logger.debug('Finished MAST.')
    bulk = BulkCreateManager()
    for group in alleleGroups:
        for line in open(os.path.join(settings.WORKSPACE_DIR, queryID, str(group)+'.matching.txt')):
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
    logger.debug('Inserting results to database...')
    bulk.done()
    
    query.c_status = 'All selected steps finished'
    query.save()
    
    return redirect(str.format('/vaxign2/query/{}/vaxitop', queryID))
    
def _multiprocess(cmds, process):
    pool = multiprocessing.Pool(process)
    complete = []
    pool_output = [pool.apply_async(subprocess.call,
                    (str(cmd), ),
                    {'stderr': subprocess.PIPE, 'shell': sys.platform != "win32"}) for cmd in cmds]
    pool.close()
    pool.join()

class BulkCreateManager(object):
    def __init__(self, chunk_size=100):
        self._create_queues = defaultdict(list)
        self.chunk_size = chunk_size
    def _commit(self, model_class):
        model_key = model_class._meta.label
        model_class.objects.bulk_create(self._create_queues[model_key])
        self._create_queues[model_key] = []
    def add(self, obj):
        model_class = type(obj)
        model_key = model_class._meta.label
        self._create_queues[model_key].append(obj)
        if len(self._create_queues[model_key]) >= self.chunk_size:
            self._commit(model_class)
    def done(self):
        for model_name, objs in self._create_queues.items():
            if len(objs) > 0:
                self._commit(apps.get_model(model_name))