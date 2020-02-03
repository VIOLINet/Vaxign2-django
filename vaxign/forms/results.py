from django import forms
from django.forms import MultipleChoiceField

import logging
logger = logging.getLogger('console')

class MultipleChoiceFieldNoValidation(forms.Field):
    def validate(self, value):
        pass

class ResultsForm(forms.Form):

    query_id = forms.CharField(label="Query/Genome ID", max_length=11)
    
    group_name = forms.CharField(max_length=50, required=False)
    
    query_seq_id_type = forms.ChoiceField(
        choices=(('protein_accession','NCBI Protein Accession'),('protein_gi','NCBI Protein GI'),('gene_id','NCBI Gene ID'),('locus_tag','Locus Tag')),
        required=False,
    )
    query_seq_ids = forms.CharField(
        widget=forms.Textarea(attrs={
            'cols':80,
            'rows':2,
        }),
        max_length=102400,
        required=False,
    )
    query_seq_ids_file = forms.FileField(required=False)
    
    keywords_type = forms.ChoiceField(
        choices=(
            ('c_gene_symbol','Gene Symbol'),
            ('c_note','Protein Note'),
        ),
        required=False,
    )
    keywords = forms.CharField(
        widget=forms.TextInput(attrs={
            'size':30,
            'maxlength':60,
        }),
        max_length=102400,
        required=False,
    )
    
    localization = forms.MultipleChoiceField(
        label="Subcellular Localization", 
        choices=(
            ('Any','Any Localization'),
            ('Cellwall',' Cellwall '),
            ('Cytoplasmic',' Cytoplasmic '),
            ('CytoplasmicMembrane',' Cytoplasmic Membrane '),
            ('Extracellular',' Extracellular '),
            ('OuterMembrane',' Outer Membrane '),
            ('Periplasmic',' Periplasmic '),
            ('Unknown',' Unknown '),
        ),
        initial=['Any'],
        required=False,
    )
    
    tmhmm_PredHel_opt = forms.ChoiceField(
        choices=(
            ('>','>'),
            ('>=','>='),
            ('<','<'),
            ('<=','<='),
        ),
        initial=['<='],       
        required=False,
    )
        
    tmhmm_PredHel_value = forms.IntegerField(
        label="Number of Transmembrane Helices", 
        initial=1, 
        min_value=0, 
        max_value=100, 
        required=False,
    )
    tmhmm_PredHel_check = forms.BooleanField(required=False)
    
    spaan_score_opt = forms.ChoiceField(
        choices=(
            ('>','>'),
            ('>=','>='),
            ('<','<'),
            ('<=','<='),
        ),
        initial=['>='],       
        required=False,
    )
    spaan_score_value = forms.FloatField(
        label="Adhesin Probability", 
        initial=0.51, 
        min_value=0.0, 
        max_value=1.0, 
        required=False
    )
    spaan_score_check = forms.BooleanField(required=False)

    have_orthologs = MultipleChoiceFieldNoValidation(
        widget=forms.SelectMultiple(attrs={
            'onChange':"chk_opt('id_have_orthologs', 'id_have_no_orthologs'); show_selected_no('id_have_orthologs', 'have_num_orthologs_label');",
            'disabled':'disabled',
            'size':3,
            }),
        required=False,
    )
    have_no_orthologs = MultipleChoiceFieldNoValidation(
        widget=forms.SelectMultiple(attrs={
            'onChange':"chk_opt('id_have_no_orthologs', 'id_have_orthologs'); show_selected_no('id_have_no_orthologs', 'have_num_no_orthologs_label');",
            'disabled':'disabled',
            'size':3,
            }),
        required=False,
    )
    
    human_alignment = forms.ChoiceField(
        choices=(
            ('y', 'Yes'),
            ('n', 'No'),
            ('','Do not use this option'),
        ),
        initial='',
        widget=forms.RadioSelect(),
        required=False,
    )
    mouse_alignment = forms.ChoiceField(
        choices=(
            ('y', 'Yes'),
            ('n', 'No'),
            ('','Do not use this option'),
        ),
        initial='',
        widget=forms.RadioSelect(),
        required=False,
    )
    pig_alignment = forms.ChoiceField(
        choices=(
            ('y', 'Yes'),
            ('n', 'No'),
            ('','Do not use this option'),
        ),
        initial='',
        widget=forms.RadioSelect(),
        required=False,
    )
    
    def clean(self):
        cleaned_data = super(ResultsForm, self).clean()
        query_id = cleaned_data.get('query_id')
        
        group_name = cleaned_data.get('group_name')
        
        query_seq_id_type = cleaned_data.get('query_seq_id_type')
        query_seq_ids = cleaned_data.get('query_seq_ids')
        query_seq_ids_file = cleaned_data.get('query_seq_ids_file')
        
        keywords_type = cleaned_data.get('keywords_type')
        keywords = cleaned_data.get('keywords')
        
        localization = cleaned_data.get('localization')
        
        tmhmm_PredHel_opt = cleaned_data.get('tmhmm_PredHel_opt')
        tmhmm_PredHel_value = cleaned_data.get('tmhmm_PredHel_value')
        tmhmm_PredHel_check = cleaned_data.get('tmhmm_PredHel_check')
        
        spaan_score_opt = cleaned_data.get('spaan_score_opt')
        spaan_score_value = cleaned_data.get('spaan_score_value')
        spaan_score_check = cleaned_data.get('spaan_score_check')
    
        have_orthologs = cleaned_data.get('have_orthologs')
        have_no_orthologs = cleaned_data.get('have_no_orthologs')
        
        human_alignment = cleaned_data.get('human_alignment')
        mouse_alignment = cleaned_data.get('mouse_alignment')
        pig_alignment = cleaned_data.get('pig_alignment')
        
        if not query_id:
            raise forms.ValidationError('Please select a genome.')