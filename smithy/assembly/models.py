from django.db import models
from django.db.models.base import Model
from django.db.models.fields import BooleanField, FloatField
from django.utils import timezone
from django.urls import reverse

ovhngs = (
    (0, '15 overhangs'),
    (1, '20 overhangs'),
    (2, '25 overhangs'),
    (3, '30 overhangs'),
)

# abstract base classes for assemblies, parts, and primers
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
    dna_conc = models.FloatField(verbose_name='DNA concentration (nM)')
    tm = models.FloatField(verbose_name='melting temperature (C)')
    backbone_file = models.FileField(upload_to='fasta/backbones/', blank=False)
    insert_file = models.FileField(upload_to='fasta/queries/', blank=False)
    
    class Meta:
        abstract = True

    def __str__(self):
        return self.title


class AssemblyPart(models.Model):
    name = models.CharField(max_length=250, default='no_name')
    date_created = models.DateTimeField(default=timezone.now)
    database = models.CharField(max_length=250, default='none')
    length = models.PositiveIntegerField()
    length_extended = models.PositiveIntegerField()
    seq = models.CharField(max_length=10000)
    seq_extended = models.CharField(max_length=10000)
    position = models.PositiveIntegerField(default=0)

    class Meta:
        abstract = True

    def __str__(self):
        return self.name


class AssemblyPrimer(models.Model):
    name = models.CharField(max_length=250, default='no_name')
    date_created = models.DateTimeField(default=timezone.now)
    primer_type = models.CharField(max_length=100)
    sequence = models.CharField(max_length=1000)
    footprint = models.CharField(max_length=1000)
    tail = models.CharField(max_length=1000)
    tm_total = FloatField()
    tm_footprint = FloatField()
    gc = FloatField()
    hairpin = BooleanField(default=False)
    hairpin_tm = FloatField()
    hairpin_dg = FloatField()
    hairpin_dh = FloatField()
    hairpin_ds = FloatField()
    homodimer = BooleanField(default=False)
    homodimer_tm = FloatField()
    homodimer_dg = FloatField()
    homodimer_dh = FloatField()
    homodimer_ds = FloatField()

    class Meta:
        abstract = True

    def __str__(self):
        return self.name


class AssemblySolution(models.Model):
    name = models.CharField(max_length=250, default='no_name')
    date_created = models.DateTimeField(default=timezone.now)
    backbone = models.CharField(max_length=100000)
    query = models.CharField(max_length=100000)
    solution = models.CharField(max_length=100000)
    parts_count = models.PositiveIntegerField()
    primers_count = models.PositiveIntegerField()
    match = models.FloatField()

    class Meta:
        abstract = True

    def __str__(self):
        return self.name


# specific classes/models in the db for supported assemblies, parts, and primers
class GoldenGateAssembly(Assembly):
    overhangs = models.IntegerField(verbose_name='overhang count', choices=ovhngs)

    def get_absolute_url(self):
        return reverse('goldengate-detail', kwargs={'pk': self.pk})


class GibsonAssembly(Assembly):
    overlap = models.PositiveIntegerField()

    def get_absolute_url(self):
        return reverse('gibson-detail', kwargs={'pk': self.pk})   


class GibsonSolution(AssemblySolution):
    assembly = models.ForeignKey(GibsonAssembly, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('gibson-solution-detail', kwargs={'pk': self.pk})


class GoldenGateSolution(AssemblySolution):
    assembly = models.ForeignKey(GoldenGateAssembly, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('goldengate-solution-detail', kwargs={'pk': self.pk})


class GoldenGatePart(AssemblyPart):
    solution = models.ForeignKey(GoldenGateSolution, on_delete=models.CASCADE)
    
    def get_absolute_url(self):
        return reverse('goldengate-part-detail', kwargs={'pk': self.pk})


class GibsonPart(AssemblyPart):
    solution = models.ForeignKey(GibsonSolution, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('gibson-part-detail', kwargs={'pk': self.pk})


class GoldenGatePrimer(AssemblyPrimer):
    part = models.ForeignKey(GoldenGatePart, on_delete=models.CASCADE)
    
    def get_absolute_url(self):
        return reverse('goldengate-primer-detail', kwargs={'pk': self.pk})


class GibsonPrimer(AssemblyPrimer):
    part = models.ForeignKey(GibsonPart, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('gibson-primer-detail', kwargs={'pk': self.pk})

