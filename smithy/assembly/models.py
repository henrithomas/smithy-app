from django.core import validators
from django.db import models
from django.db.models.base import Model
from django.db.models.fields import BooleanField, FloatField
from django.core.validators import MinValueValidator, FileExtensionValidator
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _
from django.utils import timezone
from django.urls import reverse

from tests import cut_locations

ovhngs = (
    (0, '15 overhangs - 98.5% fidelity'),
    (1, '20 overhangs - 98.1% fidelity'),
    (2, '25 overhangs - 95.8% fidelity'),
    (3, '30 overhangs - 91.7% fidelity'),
)

def fasta_validation(fa_file):
    if fa_file.file.content_type != 'application/octet-stream':
        raise ValidationError(_('Please choose a valid .fasta file.'), code='invalid')

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
    min_blast = models.PositiveIntegerField(verbose_name='min BLAST seq size (nt)',validators=[MinValueValidator(50)])
    max_blast = models.PositiveIntegerField(verbose_name='max BLAST seq size (nt)',validators=[MinValueValidator(50)])
    min_synth = models.PositiveIntegerField(verbose_name='min synthetic seq size (nt)',validators=[MinValueValidator(50)])
    max_synth = models.PositiveIntegerField(verbose_name='max synthetic seq size (nt)',validators=[MinValueValidator(50)])
    mv_conc = models.FloatField(verbose_name='monovalent ion concentration (mM)')
    dv_conc = models.FloatField(verbose_name='divalent ion concentration (mM)')
    dntp_conc = models.FloatField(verbose_name='dNTP concentration (mM)')
    dna_conc = models.FloatField(verbose_name='DNA concentration (nM)')
    tm = models.FloatField(verbose_name='melting temperature (C)')
    backbone_file = models.FileField(
                        upload_to='fasta/backbones/', 
                        validators=[
                            FileExtensionValidator(allowed_extensions=['fasta', 'fa', 'faa']),
                            fasta_validation
                        ],
                        blank=False)
    insert_file = models.FileField(
                        upload_to='fasta/queries/',
                        validators=[
                            FileExtensionValidator(allowed_extensions=['fasta', 'fa', 'faa']),
                            fasta_validation
                        ],
                        blank=False)
    multi_query = models.BooleanField(default=False)
    
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
    query_start = models.PositiveIntegerField(default=0)
    query_end = models.PositiveIntegerField(default=0)
    subject_start = models.PositiveIntegerField(default=0)
    subject_end = models.PositiveIntegerField(default=0)
    part_map = models.ImageField(default='default-part.png', upload_to='part-maps')

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
    plasmid_map = models.ImageField(default='default-plasmid.png', upload_to='assembly-maps')

    class Meta:
        abstract = True

    def __str__(self):
        return self.name


# specific classes/models in the db for supported assemblies, parts, and primers
class GoldenGateAssembly(Assembly):
    overhangs = models.IntegerField(verbose_name='overhang count', choices=ovhngs)
    scarless = models.BooleanField(default=False)
   
    def get_absolute_url(self):
        return reverse('goldengate-detail', kwargs={'pk': self.pk})


class GibsonAssembly(Assembly):
    overlap = models.PositiveIntegerField()

    def get_absolute_url(self):
        return reverse('gibson-detail', kwargs={'pk': self.pk})   


class BioBricksAssembly(Assembly):

    def get_absolute_url(self):
        return reverse('biobricks-detail', kwargs={'pk': self.pk})


class PCRAssembly(Assembly):
    overlap = models.PositiveIntegerField()

    def get_absolute_url(self):
        return reverse('pcr-detail', kwargs={'pk': self.pk})


class SLICAssembly(Assembly):
    overlap = models.PositiveIntegerField()

    def get_absolute_url(self):
        return reverse('slic-detail', kwargs={'pk': self.pk})


# class YeastAssembly(Assembly):
#     overlap = models.PositiveIntegerField()

#     def get_absolute_url(self):
#         return reverse('yeast-detail', kwargs={'pk': self.pk}) 


class GibsonSolution(AssemblySolution):
    assembly = models.ForeignKey(GibsonAssembly, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('gibson-solution-detail', kwargs={'pk': self.pk})


class GoldenGateSolution(AssemblySolution):
    assembly = models.ForeignKey(GoldenGateAssembly, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('goldengate-solution-detail', kwargs={'pk': self.pk})


class BioBricksSolution(AssemblySolution):
    assembly = models.ForeignKey(BioBricksAssembly, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('biobricks-solution-detail', kwargs={'pk': self.pk})


class PCRSolution(AssemblySolution):
    assembly = models.ForeignKey(PCRAssembly, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('pcr-solution-detail', kwargs={'pk': self.pk})


class SLICSolution(AssemblySolution):
    assembly = models.ForeignKey(SLICAssembly, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('slic-solution-detail', kwargs={'pk': self.pk})


class GoldenGatePart(AssemblyPart):
    solution = models.ForeignKey(GoldenGateSolution, on_delete=models.CASCADE)
    cuts = models.PositiveIntegerField()
    cut_locations = models.CharField(max_length=10000, default='{}')
    
    def get_absolute_url(self):
        return reverse('goldengate-part-detail', kwargs={'pk': self.pk})


class GibsonPart(AssemblyPart):
    solution = models.ForeignKey(GibsonSolution, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('gibson-part-detail', kwargs={'pk': self.pk})


class BioBricksPart(AssemblyPart):
    solution = models.ForeignKey(BioBricksSolution, on_delete=models.CASCADE)
    cuts = models.PositiveIntegerField()
    cut_locations = models.CharField(max_length=10000, default='{}')

    def get_absolute_url(self):
        return reverse('biobricks-part-detail', kwargs={'pk': self.pk})


class PCRPart(AssemblyPart):
    solution = models.ForeignKey(PCRSolution, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('pcr-part-detail', kwargs={'pk': self.pk})


class SLICPart(AssemblyPart):
    solution = models.ForeignKey(SLICSolution, on_delete=models.CASCADE)

    def get_absolute_url(self):
         return reverse('slic-part-detail', kwargs={'pk': self.pk})


class GoldenGatePrimer(AssemblyPrimer):
    part = models.ForeignKey(GoldenGatePart, on_delete=models.CASCADE)
    
    def get_absolute_url(self):
        return reverse('goldengate-primer-detail', kwargs={'pk': self.pk})


class GibsonPrimer(AssemblyPrimer):
    part = models.ForeignKey(GibsonPart, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('gibson-primer-detail', kwargs={'pk': self.pk})


class BioBricksPrimer(AssemblyPrimer):
    part = models.ForeignKey(BioBricksPart, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('biobricks-primer-detail', kwargs={'pk': self.pk})


class PCRPrimer(AssemblyPrimer):
    part = models.ForeignKey(PCRPart, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('pcr-primer-detail', kwargs={'pk': self.pk})


class SLICPrimer(AssemblyPrimer):
    part = models.ForeignKey(SLICPart, on_delete=models.CASCADE)

    def get_absolute_url(self):
        return reverse('slic-primer-detail', kwargs={'pk': self.pk})