from django.core import validators
from django.db import models
from django.db.models.base import Model
from django.db.models.fields import BooleanField, FloatField
from django.core.validators import MinValueValidator, FileExtensionValidator
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _
from django.utils import timezone
from django.urls import reverse

from django.db.models.signals import post_delete, pre_save
from django.dispatch import receiver
from django.db import models

ovhngs = (
    (0, '15 overhangs - 98.5% fidelity'),
    (1, '20 overhangs - 98.1% fidelity'),
    (2, '25 overhangs - 95.8% fidelity'),
    (3, '30 overhangs - 91.7% fidelity'),
)

enzyme_options = (
    (0.0, 'None'),
    (100.0, 'Small - 10x'),
    (500.0, 'Large - 50x')
)

from django.db.models.signals import post_delete
from django.dispatch import receiver
from django.db import models
 

"""
    model file deletion code from: https://cmljnelson.blog/2020/06/22/delete-files-when-deleting-models-in-django/  
""" 
""" 
    Whenever ANY model is deleted, if it has a file field on it, delete the associated file too
"""
@receiver(post_delete)
def delete_files_when_row_deleted_from_db(sender, instance, **kwargs):
    for field in sender._meta.concrete_fields:
        if isinstance(field, models.FileField) or isinstance(field, models.ImageField):
            instance_file_field = getattr(instance, field.name)
            delete_file_if_unused(sender, instance, field, instance_file_field)
            
""" Delete the file if something else get uploaded in its place"""
@receiver(pre_save)
def delete_files_when_file_changed(sender, instance, **kwargs):
    # Don't run on initial save
    if not instance.pk:
        return
    for field in sender._meta.concrete_fields:
        if isinstance(field, models.FileField) or isinstance(field, models.ImageField):
            #its got a file field. Let's see if it changed
            try:
                instance_in_db = sender.objects.get(pk=instance.pk)
            except sender.DoesNotExist:
                # We are probably in a transaction and the PK is just temporary
                # Don't worry about deleting attachments if they aren't actually saved yet.
                return
            instance_in_db_file_field = getattr(instance_in_db, field.name)
            instance_file_field = getattr(instance, field.name)
            if instance_in_db_file_field.name != instance_file_field.name:
                delete_file_if_unused(sender, instance, field, instance_in_db_file_field)


""" 
    Only delete the file if no other instances of that model are using it
"""    
def delete_file_if_unused(model, instance, field, instance_file_field):
    dynamic_field = {}
    dynamic_field[field.name] = instance_file_field.name
    other_refs_exist = model.objects.filter(**dynamic_field).exclude(pk=instance.pk).exists()
    if not other_refs_exist:
        instance_file_field.delete(False)

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
    primer_cost = models.FloatField(verbose_name='oligo cost ($)', default=0.0)
    part_cost = models.FloatField(verbose_name='block cost ($)', default=0.0)
    gene_cost = models.FloatField(verbose_name='mega block cost ($)', default=0.0)
    plasmid_cost = models.FloatField(verbose_name='plasmid cost ($)', default=0.0)
    pcr_polymerase_cost = models.FloatField(verbose_name='PCR polymerase cost ($)', default=0.0)
    pcr_polymerase_n_reacts = models.PositiveIntegerField(verbose_name='PCR polymerase number of reactions', default=1, validators=[MinValueValidator(1)])
    pcr_ps = models.FloatField(verbose_name='PCR probability of success', default=0.0)
    cost_pref = models.FloatField(verbose_name='cost preference', default=0.0)
    parts_pref = models.FloatField(verbose_name='parts count preference', default=0.0)

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
    parts_file = models.FileField(default='default-parts.csv', upload_to='csv/parts/')
    primers_file = models.FileField(default='default-primers.csv', upload_to='csv/primers/')
    synth_amount = models.FloatField(default=0.0)
    re_enzymes = models.BooleanField(default=False)
    part_length_average = models.PositiveIntegerField(default=0)
    primer_length_average = models.PositiveIntegerField(default=0)
    longest_part = models.PositiveIntegerField(default=0)
    shortest_part = models.PositiveIntegerField(default=0)
    tm_average = models.FloatField(default=0.0)
    hours = models.FloatField(default=0.0)
    experiments = models.PositiveIntegerField(default=0)
    solution_length = models.PositiveIntegerField(default=0)
    db_parts = models.PositiveIntegerField(default=0)
    synth_parts = models.PositiveIntegerField(default=0)
    cost_summary = models.CharField(max_length=10000, default='{}')
    time_summary = models.CharField(max_length=10000, default='{}')
    risk_summary = models.CharField(max_length=10000, default='{}')
    order_file = models.FileField(default='default-order.csv', upload_to='csv/orders/')

    class Meta:
        abstract = True

    def __str__(self):
        return self.name


# specific classes/models in the db for supported assemblies, parts, and primers
class GoldenGateAssembly(Assembly):
    overhangs = models.IntegerField(verbose_name='overhang count', choices=ovhngs)
    scarless = models.BooleanField(default=False)
    re_cost = models.FloatField(verbose_name='Type2s RE cost ($)', default=0.0)
    ligase_cost = models.FloatField(verbose_name='ligase cost ($)', default=0.0)
    re_n_reacts = models.PositiveIntegerField(verbose_name='Type2s RE number of reactions', default=1, validators=[MinValueValidator(1)])
    ligase_n_reacts = models.PositiveIntegerField(verbose_name='ligase number of reactions', default=1, validators=[MinValueValidator(1)])
    assembly_ps = models.FloatField(verbose_name='one-pot digestion-ligation probability of success', default=0.0)

    def get_absolute_url(self):
        return reverse('goldengate-detail', kwargs={'pk': self.pk})


class GibsonAssembly(Assembly):
    overlap = models.PositiveIntegerField()
    exonuclease_cost = models.FloatField(verbose_name='exonuclease cost ($)', default=0.0)
    ligase_cost = models.FloatField(verbose_name='ligase cost ($)', default=0.0)
    polymerase_cost = models.FloatField(verbose_name='polymerase cost ($)', default=0.0)
    exonuclease_n_reacts = models.PositiveIntegerField(verbose_name='exonuclease number of reactions', default=1, validators=[MinValueValidator(1)])
    ligase_n_reacts = models.PositiveIntegerField(verbose_name='ligase number of reactions', default=1, validators=[MinValueValidator(1)])
    polymerase_n_reacts = models.PositiveIntegerField(verbose_name='polymerase number of reactions', default=1, validators=[MinValueValidator(1)])
    assembly_ps = models.FloatField(verbose_name='gibson assembly probability of success', default=0.0)

    def get_absolute_url(self):
        return reverse('gibson-detail', kwargs={'pk': self.pk})   


class BioBricksAssembly(Assembly):
    EcoRI_cost = models.FloatField(verbose_name='EcoRI cost ($)', default=0.0)
    XbaI_cost = models.FloatField(verbose_name='XbaI cost ($)', default=0.0)
    SpeI_cost = models.FloatField(verbose_name='SpeI cost ($)', default=0.0)
    PstI_cost = models.FloatField(verbose_name='PstI cost ($)', default=0.0)
    ligase_cost = models.FloatField(verbose_name='ligase cost ($)', default=0.0)
    EcoRI_n_reacts = models.PositiveIntegerField(verbose_name='EcoRI number of reactions', default=1, validators=[MinValueValidator(1)])
    XbaI_n_reacts = models.PositiveIntegerField(verbose_name='XbaI number of reactions', default=1, validators=[MinValueValidator(1)])
    SpeI_n_reacts = models.PositiveIntegerField(verbose_name='SpeI number of reactions', default=1, validators=[MinValueValidator(1)])
    PstI_n_reacts = models.PositiveIntegerField(verbose_name='PstI number of reactions', default=1, validators=[MinValueValidator(1)])
    ligase_n_reacts = models.PositiveIntegerField(verbose_name='ligase number of reactions', default=1, validators=[MinValueValidator(1)])
    digestion_ps = models.FloatField(verbose_name='digestion probability of success', default=0.0)
    ligation_ps = models.FloatField(verbose_name='ligation probability of success', default=0.0)

    def get_absolute_url(self):
        return reverse('biobricks-detail', kwargs={'pk': self.pk})


class PCRAssembly(Assembly):
    overlap = models.PositiveIntegerField()
    polymerase_cost = models.FloatField(verbose_name='polymerase cost ($)', default=0.0)

    def get_absolute_url(self):
        return reverse('pcr-detail', kwargs={'pk': self.pk})


class SLICAssembly(Assembly):
    overlap = models.PositiveIntegerField()
    exonuclease_cost = models.FloatField(verbose_name='exonuclease cost ($)', default=0.0)
    ligase_cost = models.FloatField(verbose_name='ligase cost ($)', default=0.0)
    exonuclease_n_reacts = models.PositiveIntegerField(verbose_name='exonuclease number of reactions', default=1, validators=[MinValueValidator(1)])
    ligase_n_reacts = models.PositiveIntegerField(verbose_name='ligase number of reactions', default=1, validators=[MinValueValidator(1)])
    chewback_ps = models.FloatField(verbose_name='chewback probability of success', default=0.0)
    ligation_ps = models.FloatField(verbose_name='ligation probability of success', default=0.0)

    def get_absolute_url(self):
        return reverse('slic-detail', kwargs={'pk': self.pk})


class AssemblyBundle(models.Model):
    title = models.CharField(max_length=250)
    date_created = models.DateTimeField(default=timezone.now)
    description = models.CharField(max_length=1000)
    gibson = models.ManyToManyField(GibsonAssembly)
    goldengate = models.ManyToManyField(GoldenGateAssembly)
    slic = models.ManyToManyField(SLICAssembly)
    pcr = models.ManyToManyField(PCRAssembly)
    biobricks = models.ManyToManyField(BioBricksAssembly)


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
