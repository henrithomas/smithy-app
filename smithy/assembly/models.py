from django.db import models
from django.utils import timezone
from django.urls import reverse

ovhngs = (
        (15, '15 overhangs'),
        (20, '20 overhangs'),
        (25, '25 overhangs'),
        (30, '30 overhangs'),
    )

class Assembly(models.Model):
    title = models.CharField(max_length=250)
    date_created = models.DateTimeField(default=timezone.now)
    backbone = models.CharField(max_length=10000)
    insert = models.CharField(max_length=10000)
    multipart = models.BooleanField(default=False)
    addgene = models.BooleanField(default=False)
    igem = models.BooleanField(default=False)
    dnasu = models.BooleanField(default=False)
    min_blast = models.PositiveIntegerField(verbose_name='min BLAST seq size (nt)')
    max_blast = models.PositiveIntegerField(verbose_name='max BLAST seq size (nt)')
    min_synth = models.PositiveIntegerField(verbose_name='min synthetic seq size (nt)')
    max_synth = models.PositiveIntegerField(verbose_name='max synthetic seq size (nt)')
    mv_conc = models.FloatField(verbose_name='monovalent ion concentration (mM)')
    dv_conc = models.FloatField(verbose_name='divalent ion concentration (mM)')
    dntp_conc = models.FloatField(verbose_name='dNTP concentration (mM)')
    dna_conc = models.FloatField(verbose_name='dna concentration (nM)')
    tm = models.FloatField(verbose_name='melting temperature (C)')

    class Meta:
        abstract = True

    def __str__(self):
        return self.title


class GoldenGateAssembly(Assembly):
    overhangs = models.IntegerField(verbose_name='overhang count', choices=ovhngs)

    def get_absolute_url(self):
        return reverse('goldengate-detail', kwargs={'pk': self.pk})

class GibsonAssembly(Assembly):
    overlap = models.PositiveIntegerField()

    def get_absolute_url(self):
        return reverse('gibson-detail', kwargs={'pk': self.pk})