from rest_framework import serializers

from vaxign.models import TVaxignQuery
from vaxign.models import TUser
from vaxign.models import TVaxignAlleleGroup
from vaxign.models import TVaxignMastResults

class TVaxignQuerySerializer(serializers.ModelSerializer):
    class Meta:
        model = TVaxignQuery
        fields = ['c_query_id', 'c_species_short', 'c_species_name']
        
class TUserSerializer(serializers.ModelSerializer):
    class Meta:
        model = TUser
        fields = ['c_user_name', 'c_first_name', 'c_last_name']
        
class TVaxignAlleleGroupSerializer(serializers.ModelSerializer):
    class Meta:
        model = TVaxignAlleleGroup
        fields = ['mhc_species', 'mhc_class', 'mhc_allele', 'epitope_length']

class TVaxignMastResultsSerializer(serializers.Serializer):
    c_allele_group_id = serializers.IntegerField()
    c_sequence_id = serializers.CharField()
    c_hit_start = serializers.IntegerField()
    c_hit_end = serializers.IntegerField()
    c_hit_p_value = serializers.FloatField()
    protein = serializers.CharField()
    epitope = serializers.CharField()
    mhc_allele = serializers.CharField()

class IEDBEpitopeSerializer(serializers.Serializer):
    epitope_id = serializers.IntegerField()
    linear_peptide_seq = serializers.CharField()
    alleles = serializers.DictField()
