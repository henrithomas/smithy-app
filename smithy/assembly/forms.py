from django import forms
from django.utils import timezone
from django.core.validators import MinValueValidator, FileExtensionValidator
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _
from django.urls import reverse
from pkg_resources import require


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
    gib_overlap = forms.IntegerField(required=False)
    gib_exonuclease_cost = forms.FloatField(required=False)
    gib_ligase_cost = forms.FloatField(required=False)
    gib_ligase_n_reacts = forms.IntegerField(required=False)
    gib_exonuclease_n_reacts = forms.IntegerField(required=False)
    gib_polymerase_cost = forms.FloatField(required=False)
    gib_polymerase_n_reacts = forms.IntegerField(required=False)
    gib_assembly_ps = forms.FloatField(required=False)
    bb_EcoRI_cost = forms.FloatField(required=False)
    bb_XbaI_cost = forms.FloatField(required=False)
    bb_SpeI_cost = forms.FloatField(required=False)
    bb_PstI_cost = forms.FloatField(required=False)
    bb_ligase_cost = forms.FloatField(required=False)
    bb_EcoRI_n_reacts = forms.IntegerField(required=False)
    bb_XbaI_n_reacts = forms.IntegerField(required=False)
    bb_SpeI_n_reacts = forms.IntegerField(required=False)
    bb_PstI_n_reacts = forms.IntegerField(required=False)
    bb_ligase_n_reacts = forms.IntegerField(required=False)
    bb_digestion_ps = forms.FloatField(required=False)
    bb_ligation_ps = forms.FloatField(required=False)
    gg_re_cost = forms.FloatField(required=False)
    gg_ligase_cost = forms.FloatField(required=False)
    gg_re_n_reacts = forms.IntegerField(required=False)
    gg_ligase_n_reacts = forms.IntegerField(required=False)
    gg_assembly_ps = forms.FloatField(required=False)
    pcr_overlap = forms.IntegerField(required=False)
    slic_overlap = forms.IntegerField(required=False)
    # slic_exonuclease_cost = forms.FloatField(required=False)
    # slic_ligase_cost = forms.FloatField(required=False)
    # slic_exonuclease_n_reacts = forms.IntegerField(required=False)
    # slic_ligase_n_reacts = forms.IntegerField(required=False)
    slic_chewback_ps = forms.FloatField(required=False)
    # slic_ligation_ps = forms.FloatField(required=False)
    primer_cost = forms.FloatField(required=False)
    part_cost = forms.FloatField(required=False)
    gene_cost = forms.FloatField(required=False)
    parts_pref = forms.FloatField(required=False)
    cost_pref = forms.FloatField(required=False)
    pcr_polymerase_cost = forms.FloatField(required=False)
    pcr_polymerase_n_reacts = forms.IntegerField(required=False)
    pcr_ps = forms.FloatField(required=False)
    gib_mmix_cost = forms.FloatField(required=False)
    gg_mmix_cost = forms.FloatField(required=False)
    pcr_mmix_cost = forms.FloatField(required=False)
    gib_mmix_n_reacts = forms.IntegerField(required=False)
    gg_mmix_n_reacts = forms.IntegerField(required=False)
    pcr_mmix_n_reacts = forms.IntegerField(required=False)