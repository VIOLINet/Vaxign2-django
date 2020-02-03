from django import forms

import logging
from Bio.SubsMat.MatrixInfo import grant
logger = logging.getLogger('console')

class ProjectsForm(forms.Form):
    
    project_id = forms.IntegerField(
        widget=forms.HiddenInput(), 
        min_value=0, 
        max_value=10**10,
        required=False,
        )
    
    project_name = forms.CharField(
        label='Project Title*',
        min_length=3, 
        max_length=255,
    )
    
    is_public = forms.BooleanField(
        initial=False,
        label='Share your project to the public?',
        help_text='Anyone can select this project in the Vaxign Precompute Query.',
        widget=forms.CheckboxInput(attrs={
            'class':'form-check form-check-inline',
            'style':'display:inline-block;width:5%',
        }),
        required=False,
    )
    
    list_project = forms.BooleanField(
        initial=True,
        label='List your project in the Vaxign Precompute Query?',
        help_text='After login, you can select this project in the Vaxign Precompute Query.',
        widget=forms.CheckboxInput(attrs={
            'class':'form-check form-check-inline',
            'style':'display:inline-block;width:5%',
        }),
        required=False,
    )
    
    description = forms.CharField(
        label='Project Description',
        widget=forms.Textarea(attrs={
            'cols':80,
            'rows':6,
        }),
        max_length=5120,
        required=False,
    )
    
    institution = forms.CharField(
        initial='',
        label='Institution(s)',
        widget=forms.TextInput(attrs={
            'size':50,
            'maxlength':60,
        }),
        max_length=100,
        required=False,
    )
    
    investigator = forms.CharField(
        initial='',
        label='Investigator(s)',
        widget=forms.TextInput(attrs={
            'size':50,
            'maxlength':60,
        }),
        max_length=100,
        required=False,
    )
    
    grant = forms.CharField(
        initial='',
        label='Grant Support',
        widget=forms.TextInput(attrs={
            'size':50,
            'maxlength':60,
        }),
        max_length=100,
        required=False,
    )
    
    note = forms.CharField(
        initial='',
        label='Project Note',
        widget=forms.Textarea(attrs={
            'cols':80,
            'rows':3,
        }),
        max_length=2048,
        required=False,
    )
    
    def clean(self):
        cleaned_data = super(ProjectsForm, self).clean()
        
        project_id = cleaned_data.get('project_id')
        project_name = cleaned_data.get('project_name')
        
        is_public = cleaned_data.get('is_public')
        list_project = cleaned_data.get('list_project')
        
        description = cleaned_data.get('description')
        institution = cleaned_data.get('institution')
        investigator = cleaned_data.get('investigator')
        grant = cleaned_data.get('grant')
        note = cleaned_data.get('note')
        
        if not project_name:
            raise forms.ValidationError('Please provide a project name.')
        