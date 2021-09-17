from django import forms
from django.utils import timezone
from django.core.validators import MinValueValidator, FileExtensionValidator
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _
from django.urls import reverse


ovhngs = (
    (0, '15 overhangs - 98.5% fidelity'),
    (1, '20 overhangs - 98.1% fidelity'),
    (2, '25 overhangs - 95.8% fidelity'),
    (3, '30 overhangs - 91.7% fidelity'),
)

def fasta_validation(fa_file):
    if fa_file.content_type != 'application/octet-stream':
        raise ValidationError(_('Please choose a valid .fasta file.'), code='invalid')


class BundleForm(forms.Form):
    title = forms.CharField(max_length=250)
    description = forms.CharField(max_length=1000)
    multipart = forms.BooleanField(required=False)
    addgene = forms.BooleanField(required=False)
    igem = forms.BooleanField(required=False)
    dnasu = forms.BooleanField(required=False)
    min_blast = forms.IntegerField(validators=[MinValueValidator(50)])
    max_blast = forms.IntegerField(validators=[MinValueValidator(50)])
    min_synth = forms.IntegerField(validators=[MinValueValidator(50)])
    max_synth = forms.IntegerField(validators=[MinValueValidator(50)])
    mv_conc = forms.FloatField()
    dv_conc = forms.FloatField()
    dntp_conc = forms.FloatField()
    dna_conc = forms.FloatField()
    tm = forms.FloatField()
    backbone_file = forms.FileField( 
                        validators=[
                            FileExtensionValidator(allowed_extensions=['fasta', 'fa', 'faa']),
                            fasta_validation
                        ])
    insert_file = forms.FileField(
                        validators=[
                            FileExtensionValidator(allowed_extensions=['fasta', 'fa', 'faa']),
                            fasta_validation
                        ])
    multi_query = forms.BooleanField(required=False)
    overhangs = forms.ChoiceField(choices=ovhngs, required=False)
    scarless = forms.BooleanField(required=False) 
    overlap = forms.IntegerField(required=False)
    gibson = forms.BooleanField(required=False)
    goldengate = forms.BooleanField(required=False)
    slic = forms.BooleanField(required=False)
    pcr = forms.BooleanField(required=False)
    biobricks = forms.BooleanField(required=False)