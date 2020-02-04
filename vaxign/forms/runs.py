from django import forms

import logging
logger = logging.getLogger('console')

class RunsForm(forms.Form):
    
    group_name = forms.CharField(
        widget=forms.TextInput(attrs={
            'class':'dropdown-input',
            'size':50,
        }),
        max_length=128,
        required=False,
    )
    
    genome_name = forms.CharField(
        widget=forms.TextInput(attrs={
            'size':50,
        }),
        max_length=128,
        required=False,
    )
    
    sequence_type = forms.ChoiceField(
        choices=(
            ('protein_fasta', 'Protein Sequence (FASTA Format)'),
            ('protein_uniprotkb', 'UniprotKB Protein ID'),
            ('protein_uniprot_proteome', 'Uniprot Proteome ID'),
            ('protein_gi','NCBI Protein ID'),
            ('protein_refseq','NCBI Protein Refseq'),
            ('gene_id','NCBI Gene ID'),
            ('bioproject_id','NCBI Bioproject ID'),
            ('protein_fasta_url', 'Protein Sequence (FASTA File Link)'),
        ),
    )
    
    sequence = forms.CharField(
        widget=forms.Textarea(attrs={
            'cols':120,
            'rows':8,
        }),
        max_length=8*1024*1024,
    )
    
    sequence_file = forms.FileField(required=False)
    
    organism = forms.ChoiceField(
        choices = (
            ('bacteria', 'Bacterium'),
            ('virus', 'Virus'),
            ('parasite','Parasite'),
        ),
    )
    
    bacteria_strain = forms.ChoiceField(
        choices = (
            ('--negative','Gram negative bacterium'),
            ('--positive','Gram positive bacterium'),
        ),
    )
    
    basic_analysis = forms.MultipleChoiceField(
        choices=(
            ('psort','Subcellular Localization'),
            ('tmhmm','Transmembrane Helix'),
            ('spaan','Adhesion Probability'),
            ('blasth','Similarity to Human Proteins'),
            ('blastm',' Similarity to Mouse Proteins'),
            ('blastp',' Similarity to Pig Proteins'),
        ),
        initial=['vaxign_ml','vaxitop','psort','tmhmm','spaan'],
        widget=forms.CheckboxSelectMultiple(attrs={
            'class':'col2',
        }),
        required=False,
    )
    
    vaxignml_analysis = forms.ChoiceField(
        choices = (
            ('y', 'Yes'),
            ('n', 'No'),
        ),
        initial='y',
        widget=forms.RadioSelect(),
    )
    
    vaxitop_analysis = forms.ChoiceField(
        choices = (
            ('y', 'Yes'),
            ('n', 'No'),
        ),
        initial='y',
        widget=forms.RadioSelect(),
    )
    
    note = forms.CharField(
        widget=forms.Textarea(attrs={
            'cols':120,
            'rows':3,
        }),
        max_length=128,
        required=False,
    )
    
    email = forms.EmailField(
        widget=forms.EmailInput(attrs={
            'size':50,
            'maxlength':60,
        }),
        max_length=128,
        required=False,
    )
    
    def clean(self):
        cleaned_data = super(RunsForm, self).clean()
        sequence = cleaned_data.get('sequence')
        
        sequence_type = cleaned_data.get('sequence_type')
        sequence_file = cleaned_data.get('sequence_file')
        
        organism = cleaned_data.get('organism')
        bacteria_strain = cleaned_data.get('bacteria_strain')
        
        basic_analysis = cleaned_data.get('basic_analysis')
        vaxignml_analysis = cleaned_data.get('vaxignml_analysis')
        vaxitop_analysis = cleaned_data.get('vaxitop_analysis')
        
        note = cleaned_data.get('note')
        email = cleaned_data.get('email')
        
        if not sequence or (sequence_type == 'protein_fasta' and not sequence.startswith('>')):
            raise forms.ValidationError('Please provide valid protein sequence FASTA or identifiers.')
        
        if not organism:
            raise forms.ValidationError('Please select an organism.')
        else:
            if organism == 'bacteria':
                if not bacteria_strain:
                    raise forms.ValidationError('Please select gram for bacterium.')
        
        if not basic_analysis and vaxignml_analysis == 'n'  and vaxitop_analysis == 'n':
            raise forms.ValidationError('Please select at least one analysis.')
        