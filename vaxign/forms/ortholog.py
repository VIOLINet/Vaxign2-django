from django import forms
from django.forms import MultipleChoiceField

import logging
logger = logging.getLogger('console')

class MultipleChoiceFieldNoValidation(forms.Field):
    def validate(self, value):
        pass

class OrthologForm(forms.Form):

    query_id = forms.CharField(widget=forms.HiddenInput(), max_length=11)
    
    group_name = forms.CharField(widget=forms.HiddenInput(), max_length=50)

    have_orthologs = MultipleChoiceFieldNoValidation(
        widget=forms.SelectMultiple(attrs={
            'onChange':"show_selected_no('id_have_orthologs', 'have_num_orthologs_label');",
            'disabled':'disabled',
            'size':3,
            }),
        required=False,
    )
    
    def clean(self):
        cleaned_data = super(OrthologForm, self).clean()
        query_id = cleaned_data.get('query_id')
        
        group_name = cleaned_data.get('group_name')
        
        have_orthologs = cleaned_data.get('have_orthologs')
        
        if not query_id:
            raise forms.ValidationError('Please select a genome.')
        elif not group_name:
            raise forms.ValidationError('Please select a genome group.')