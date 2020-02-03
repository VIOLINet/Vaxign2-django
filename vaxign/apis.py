import re
import pymysql

from rest_framework.decorators import api_view
from rest_framework.response import Response

from django.conf import settings

from django.db.models import Q
from vaxign.models import TVaxignQuery
from vaxign.models import TUserVaxignPrivilege
from vaxign.models import TVaxignProject
from vaxign.models import TUser
from vaxign.models import TVaxignAnalysisResults
from vaxign.models import TVaxignAlleleGroup
from vaxign.models import TVaxignMastResults

from vaxign.serializers import TVaxignQuerySerializer
from vaxign.serializers import TUserSerializer
from vaxign.serializers import IEDBEpitopeSerializer
from vaxign.serializers import TVaxignAlleleGroupSerializer
from vaxign.serializers import TVaxignMastResultsSerializer

from django.http import HttpResponse

import logging
logger = logging.getLogger('console')

@api_view()
def t_vaxign_query_group_all(request):
    if 'c_user_name' in request.session:
        privilege = TUserVaxignPrivilege.objects.filter(
            Q(c_userfk=request.session['c_user_name'])
        )
        projects = TVaxignProject.objects.filter(
            Q(c_vaxign_projectid__in=privilege.values_list('c_projectfk', flat=True)) & Q(c_list_title=1)
        )
        tmpQueries = {}
        for query in TVaxignQuery.objects.filter(
            Q(c_is_public=1) |
            Q(c_projectfk__in=projects.values_list('c_vaxign_projectid', flat=True) )
        ):
            if query.c_species_short == '':
                continue
            tmpQueries[query.c_species_short] = query.c_query_id
        queries = TVaxignQuery.objects.filter(c_query_id__in=tmpQueries.values()).order_by('c_species_short')
    else:
        tmpQueries = {}
        for query in TVaxignQuery.objects.filter(
            Q(c_is_public=1)
        ):
            if query.c_species_short == '':
                continue
            tmpQueries[query.c_species_short] = query.c_query_id
        queries = TVaxignQuery.objects.filter(c_query_id__in=tmpQueries.values()).order_by('c_species_short')
    
    return Response(TVaxignQuerySerializer(queries, many=True).data)

@api_view()
def t_vaxign_query_group(request, genomeGroup):
    
    if 'c_user_name' in request.session:
        privilege = TUserVaxignPrivilege.objects.filter(
            Q(c_userfk=request.session['c_user_name'])
        )
        projects = TVaxignProject.objects.filter(
            Q(c_vaxign_projectid__in=privilege.values_list('c_projectfk', flat=True)) & Q(c_list_title=1)
        )
        queries = TVaxignQuery.objects.filter(
            Q(c_is_public=1) & Q(c_species_short=genomeGroup) |
            Q(c_projectfk__in=projects.values_list('c_vaxign_projectid', flat=True) ) & Q(c_species_short=genomeGroup)
        ).order_by('c_species_short', 'c_species_name')
    else:
        queries = TVaxignQuery.objects.filter(
            Q(c_is_public=1) & Q(c_species_short=genomeGroup),
        ).order_by('c_species_short', 'c_species_name')
    
    return Response(TVaxignQuerySerializer(queries, many=True).data)

@api_view()
def t_vaxign_query_ortholog_exclude(request, queryID):
    
    try:
        genomeGroup = TVaxignQuery.objects.get(
            Q(c_query_id=queryID) & Q(c_ortholog_computed=1)
        ).c_species_short
    except:
        genomeGroup = ''
    
    queries = TVaxignQuery.objects.filter(
        Q(c_ortholog_computed=1) & ~Q(c_query_id=queryID) & Q(c_species_short=genomeGroup)
    ).order_by('c_species_name')
    
    return Response(TVaxignQuerySerializer(queries, many=True).data)

@api_view()
def t_user_query_all(request):
    
    if 'c_user_name' not in request.session:
        return HttpResponse(status=403)
    
    users = TUser.objects.all().order_by('c_last_name', 'c_first_name')
    
    return Response(TUserSerializer(users, many=True).data)

@api_view()
def t_vaxign_allele_group(request, species, mhc):
    
    group = TVaxignAlleleGroup.objects.filter(
        Q(mhc_species=species) & Q(mhc_class=mhc)
    ).order_by('mhc_allele','epitope_length')
    
    return Response(TVaxignAlleleGroupSerializer(group, many=True).data)

@api_view()
def t_vaxign_mast_results(request, queryID):
    
    groupMap = {}
    groupIDMap = {}
    for row in TVaxignAlleleGroup.objects.all():
        groupMap[re.sub("[^a-zA-Z0-9]", '', row.mhc_allele)] = row.mhc_allele
        groupIDMap[row.c_allele_group_id] = row.mhc_allele
    
    
    selected = set([])
    for select in request.GET['groups'].split(','):
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
    sequences = {}
    for row in TVaxignAnalysisResults.objects.filter(c_query_id=queryID):
        if row.c_protein_accession is not None and row.c_protein_accession != '':
            sequences[row.c_sequence_id] = {
                'sequence': row.c_sequence,
                'protein': row.c_protein_accession,
            }
        else:
            sequences[row.c_sequence_id] = {
                'sequence': row.c_sequence,
                'protein': row.c_sequence_acc,
            }
    
    masts = []
    for row in TVaxignMastResults.objects.filter(
        Q(c_sequence_id__in=sequences.keys()) & Q(c_allele_group_id__in=selected)
    ):
        masts.append({
            'c_allele_group_id': row.c_allele_group_id,
            'c_sequence_id': row.c_sequence_id,
            'c_hit_start': row.c_hit_start,
            'c_hit_end': row.c_hit_end,
            'c_hit_p_value': str.format('{:.3f}', float(row.c_hit_p_value)),
            'protein': sequences[row.c_sequence_id]['protein'],
            'epitope': sequences[row.c_sequence_id]['sequence'][row.c_hit_start-1:row.c_hit_end],
            'mhc_allele': groupIDMap[row.c_allele_group_id],
        })
    
    return Response(TVaxignMastResultsSerializer(masts, many=True).data)

@api_view()
def t_vaxign_mast_results_one_seq(request, queryID, seqID):
    
    groupMap = {}
    groupIDMap = {}
    for row in TVaxignAlleleGroup.objects.all():
        groupMap[re.sub("[^a-zA-Z0-9]", '', row.mhc_allele)] = row.mhc_allele
        groupIDMap[row.c_allele_group_id] = row.mhc_allele
    
    
    selected = set([])
    for select in request.GET['groups'].split(','):
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
    
    sequence = TVaxignAnalysisResults.objects.get(
        Q(c_query_id=queryID) & Q(c_sequence_id=seqID)
    )
    
    masts = []
    for row in TVaxignMastResults.objects.filter(
        Q(c_sequence_id=sequence.c_sequence_id) & Q(c_allele_group_id__in=selected)
    ):
        masts.append({
            'c_allele_group_id': row.c_allele_group_id,
            'c_sequence_id': row.c_sequence_id,
            'c_hit_start': row.c_hit_start,
            'c_hit_end': row.c_hit_end,
            'c_hit_p_value': str.format('{:.3f}', float(row.c_hit_p_value)),
            'protein': sequence.c_sequence_id,
            'epitope': sequence.c_sequence[row.c_hit_start-1:row.c_hit_end],
            'mhc_allele': groupIDMap[row.c_allele_group_id],
        })
    
    return Response(TVaxignMastResultsSerializer(masts, many=True).data)

@api_view()
def iedb_epitope_t_cell_linear(request, seqID):
    
    sequence = TVaxignAnalysisResults.objects.get(c_sequence_id=seqID).c_sequence
    
    heDB = pymysql.connect(
        host = settings.MYSQL_VIOLIN_HOST_IP,
        user = settings.MYSQL_VIOLIN_USER_NAME,
        password = settings.MYSQL_VIOLIN_USER_PWD,
    )
    
    epitopes = {}
    
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
                tmp = {
                    'epitope_id': row[0],
                    'linear_peptide_seq': row[1],
                    'alleles': {},
                }
                epitopes[row[0]] = tmp
    finally:
        cursor.close()
    
    for epitope_id in epitopes.keys():
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
                    epitopes[epitope_id]['alleles'][row[0]] = row[1]
        finally:
            cursor.close()
    
    return Response(IEDBEpitopeSerializer(epitopes.values(), many=True).data)
