from django import forms

class AssemblyForm(forms.Form):
    insert = forms.CharField(label='insert', max_length=100, widget=forms.Textarea)
    backbone = forms.CharField(label='backbone', max_length=100, widget=forms.Textarea)
    multipart = forms.BooleanField(label='multipart')
    addgene = forms.BooleanField(label='AddGene', required=False)
    dnasu = forms.BooleanField(label='DNASU', required=False)
    igem = forms.BooleanField(label='iGEM', required=False)
    min_blast = forms.IntegerField(label='min blast seq')
    max_blast = forms.IntegerField(label='max blast seq')
    min_synth = forms.IntegerField(label='min synthetic seq')
    max_synth = forms.IntegerField(label='max synthetic seq')
    mv_conc = forms.FloatField(label='monovalent ion concentration')
    dv_conc = forms.FloatField(label='divalent ion concentration')
    dntp_conc = forms.FloatField(label='dntp concentration')
    dna_conc = forms.FloatField(label='dna concentration')
    tm = forms.FloatField(label='melting temperature')

class GoldenGateForm(AssemblyForm):
    ovhngs = (
        ('15', '15-overhangs'),
        ('20', '20-overhangs'),
        ('25', '25-overhangs'),
        ('30', '30-overhangs'),
    )

    overhangs = forms.ChoiceField(label='overhangs', choices=ovhngs)

class GibsonForm(AssemblyForm):
    overlap = forms.IntegerField(label='overlap')
