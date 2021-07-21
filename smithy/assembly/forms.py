from django import forms
from .models import GibsonAssembly
class GibsonForm(forms.ModelForm):
    class Meta:
        model = GibsonAssembly
        fields = [
            'title',
            'backbone_file',
            'insert_file',
            'addgene',
            'igem',
            'dnasu',
            'min_blast',
            'max_blast',
            'min_synth',
            'max_synth',
            'mv_conc',
            'dv_conc',
            'dntp_conc',
            'dntp_conc',
            'dna_conc', 
            'tm',
            'overlap']