from typing import List
from Bio.Blast import Record
from django.utils import timezone
from django.db.models.deletion import Collector
from django.shortcuts import render,redirect
from .service import (
    gibson_create_service,
    goldengate_create_service,
    biobricks_create_service, 
    pcr_create_service,
    slic_create_service,
    bundle_create_service
)
from .models import (
    AssemblyBundle,
    BioBricksAssembly,
    BioBricksPart,
    BioBricksPrimer,
    BioBricksSolution,
    GibsonAssembly, 
    GibsonPart,
    GibsonPrimer,
    GoldenGateAssembly,
    GoldenGatePart,
    GoldenGatePrimer,
    GibsonSolution,
    GoldenGateSolution,
    PCRAssembly,
    PCRSolution,
    PCRPart,
    PCRPrimer,
    SLICAssembly,
    SLICSolution,
    SLICPart,
    SLICPrimer
)
from django.views.generic import (
    DetailView,
    CreateView,
    ListView
)
from django.contrib.messages.views import SuccessMessageMixin 
import os
from .forms import BundleForm
from itertools import chain


def home(request):
    return render(request, 'assembly/home.html')

def about(request):
    return render(request, 'assembly/about.html', {'title': 'Assembly About'}) 

# Gibson Assembly
class GibsonDetailView(DetailView):
    model = GibsonAssembly
    context_object_name = 'gibson'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # TODO change how the queries happen here, maybe change to code in the GET request
        context['title'] = self.object.title 
        context['solutions'] = self.object.gibsonsolution_set.all()
        return context


class GibsonCreateView(SuccessMessageMixin, CreateView):
    model = GibsonAssembly
    success_message = 'View your new Gibson assembly below...'
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
        'overlap',
        'multi_query',
        'primer_cost',
        'part_cost',
        'gene_cost',
        'exonuclease_cost',
        'ligase_cost',
        'polymerase_cost',
        'pcr_polymerase_cost',
        'pcr_polymerase_n_reacts',
        'pcr_ps',
        'cost_pref',
        'parts_pref',
        'ligase_n_reacts',
        'exonuclease_n_reacts',
        'polymerase_n_reacts',
        'assembly_ps'
    ]

    def form_valid(self, form):
        self.object = form.save()
        gibson_create_service(self.object)
        return super().form_valid(form)  

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'New Gibson Assembly'
        return context


class GibsonSolutionDetailView(DetailView):
    model = GibsonSolution
    context_object_name = 'gibson_solution'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        context['exonuclease_cost'] = self.object.assembly.exonuclease_cost
        context['ligase_cost'] = self.object.assembly.ligase_cost
        context['polymerase_cost'] = self.object.assembly.polymerase_cost
        context['parts'] = self.object.gibsonpart_set.all()
        context['primer_sets'] =  [part.gibsonprimer_set.all() for part in context['parts']]
        return context


class GibsonPartDetailView(DetailView):
    model = GibsonPart
    context_object_name = 'gibson_part'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        context['primers'] = self.object.gibsonprimer_set.all()
        return context


class GibsonPrimerDetailView(DetailView):
    model = GibsonPrimer
    context_object_name = 'gibson_primer'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        return context


class GibsonListView(ListView):
    model = GibsonAssembly
    context_object_name = 'gibsons'

    def get_queryset(self):
        qs = super().get_queryset()
        return qs.all() \
            .only('title', 'date_created') \
            .order_by('-date_created')
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'Gibson Assemblies'
        return context

# Golden Gate Assembly
class GoldenGateDetailView(DetailView):
    model = GoldenGateAssembly
    context_object_name = 'goldengate'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.title
        context['solutions'] = self.object.goldengatesolution_set.all()
        return context


class GoldenGateCreateView(SuccessMessageMixin, CreateView):
    model = GoldenGateAssembly
    success_message = 'View your new Golden Gate assembly below...'
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
        'overhangs',
        'multi_query',
        'scarless',
        'primer_cost',
        'part_cost',
        'gene_cost',
        're_cost',
        'ligase_cost',
        'pcr_polymerase_cost',
        'pcr_polymerase_n_reacts',
        'pcr_ps',
        'cost_pref',
        'parts_pref',
        're_n_reacts',
        'ligase_n_reacts',
        'assembly_ps'
    ]

    def form_valid(self, form):
        self.object = form.save()
        goldengate_create_service(self.object)
        return super().form_valid(form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'New Golden Gate Assembly'
        return context


class GoldenGateSolutionDetailView(DetailView):
    model = GoldenGateSolution  
    context_object_name = 'goldengate_solution'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        context['re_cost'] = self.object.assembly.re_cost
        context['ligase_cost'] = self.object.assembly.ligase_cost
        context['parts'] = self.object.goldengatepart_set.all()
        context['primer_sets'] =  [part.goldengateprimer_set.all() for part in context['parts']]
        return context   


class GoldenGatePartDetailView(DetailView):
    model = GoldenGatePart
    context_object_name = 'goldengate_part'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        context['primers'] = self.object.goldengateprimer_set.all()
        return context


class GoldenGatePrimerDetailView(DetailView):
    model = GoldenGatePrimer
    context_object_name = 'goldengate_primer'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        return context


class GoldenGateListView(ListView):
    model = GoldenGateAssembly
    context_object_name = 'goldengates'

    def get_queryset(self):
        qs = super().get_queryset()
        return qs.all() \
            .only('title', 'date_created') \
            .order_by('-date_created')
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'Golden Gate Assemblies'
        return context

# Bio Bricks Assembly
class BioBricksDetailView(DetailView):
    model = BioBricksAssembly
    context_object_name = 'biobrick'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.title
        context['solutions'] = self.object.biobrickssolution_set.all()
        return context


class BioBricksCreateView(SuccessMessageMixin, CreateView):
    model = BioBricksAssembly
    success_message = 'View your new BioBricks assembly below...'
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
        'multi_query',
        'primer_cost',
        'part_cost',
        'gene_cost',
        'EcoRI_cost',
        'XbaI_cost',
        'SpeI_cost',
        'PstI_cost',
        'pcr_polymerase_cost',
        'pcr_polymerase_n_reacts',
        'pcr_ps',
        'cost_pref',
        'parts_pref',
        'EcoRI_n_reacts',
        'XbaI_n_reacts',
        'SpeI_n_reacts',
        'PstI_n_reacts',
        'digestion_ps',
        'ligation_ps',
        'ligase_cost',
        'ligase_n_reacts'
    ]

    def form_valid(self, form):
        self.object = form.save()
        biobricks_create_service(self.object)
        return super().form_valid(form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'New BioBricks Assembly'
        return context


class BioBricksSolutionDetailView(DetailView):
    model = BioBricksSolution
    context_object_name = 'biobricks_solution'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        context['EcoRI_cost'] = self.object.assembly.EcoRI_cost
        context['XbaI_cost'] = self.object.assembly.XbaI_cost
        context['SpeI_cost'] = self.object.assembly.SpeI_cost
        context['PstI_cost'] = self.object.assembly.PstI_cost
        context['parts'] = self.object.biobrickspart_set.all()
        context['primer_sets'] =  [part.biobricksprimer_set.all() for part in context['parts']]
        return context


class BioBricksPartDetailView(DetailView):
    model = BioBricksPart
    context_object_name = 'biobricks_part'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        context['primers'] = self.object.biobricksprimer_set.all()
        return context


class BioBricksPrimerDetailView(DetailView):
    model = BioBricksPrimer
    context_object_name = 'biobricks_primer'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        return context


class BioBricksListView(ListView):
    model = BioBricksAssembly
    context_object_name = 'biobricks'

    def get_queryset(self):
        qs = super().get_queryset()
        return qs.all() \
            .only('title', 'date_created') \
            .order_by('-date_created')
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'BioBricks Assemblies'
        return context

# PCR Assembly
class PCRDetailView(DetailView):
    model = PCRAssembly
    context_object_name = 'pcr'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.title
        context['solutions'] = self.object.pcrsolution_set.all()
        return context


class PCRCreateView(SuccessMessageMixin, CreateView):
    model = PCRAssembly
    success_message = 'View your new PCR-SOE assembly below...'
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
        'overlap',
        'multi_query',
        'primer_cost',
        'part_cost',
        'gene_cost',
        'pcr_polymerase_cost',
        'pcr_polymerase_n_reacts',
        'pcr_ps',
        'cost_pref',
        'parts_pref'
    ]

    def form_valid(self, form):
        self.object = form.save()
        pcr_create_service(self.object)
        return super().form_valid(form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'New PCR-SOE Assembly'
        return context


class PCRSolutionDetailView(DetailView):
    model = PCRSolution
    context_object_name = 'pcr_solution'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        context['polymerase_cost'] = self.object.assembly.polymerase_cost
        context['parts'] = self.object.pcrpart_set.all()
        context['primer_sets'] =  [part.pcrprimer_set.all() for part in context['parts']]
        return context


class PCRPartDetailView(DetailView):
    model = PCRPart
    context_object_name = 'pcr_part'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        context['primers'] = self.object.pcrprimer_set.all()
        return context


class PCRPrimerDetailView(DetailView):
    model = PCRPrimer
    context_object_name = 'pcr_primer'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        return context


class PCRListView(ListView):
    model = PCRAssembly
    context_object_name = 'pcrs'

    def get_queryset(self):
        qs = super().get_queryset()
        return qs.all() \
            .only('title', 'date_created') \
            .order_by('-date_created')
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'PCR-SOE Assemblies'
        return context

# SLIC Assembly
class SLICDetailView(DetailView):
    model = SLICAssembly
    context_object_name = 'slic'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.title
        context['solutions'] = self.object.slicsolution_set.all()
        return context


class SLICCreateView(SuccessMessageMixin, CreateView):
    model = SLICAssembly
    success_message = 'View your new SLIC assembly below...'
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
        'overlap',
        'multi_query',
        'primer_cost',
        'part_cost',
        'gene_cost',
        'exonuclease_cost',
        'ligase_cost',
        'pcr_polymerase_cost',
        'pcr_polymerase_n_reacts',
        'pcr_ps',
        'cost_pref',
        'parts_pref',
        'exonuclease_n_reacts',
        'ligase_n_reacts',
        'chewback_ps',
        'ligation_ps'
    ]

    def form_valid(self, form):
        self.object = form.save()
        slic_create_service(self.object)
        return super().form_valid(form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'New SLIC Assembly'
        return context


class SLICSolutionDetailView(DetailView):
    model = SLICSolution
    context_object_name = 'slic_solution'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        context['exonuclease_cost'] = self.object.assembly.exonuclease_cost
        context['ligase_cost'] = self.object.assembly.ligase_cost
        context['parts'] = self.object.slicpart_set.all()
        context['primer_sets'] =  [part.slicprimer_set.all() for part in context['parts']]
        return context


class SLICPartDetailView(DetailView):
    model = SLICPart
    context_object_name = 'slic_part'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        context['primers'] = self.object.slicprimer_set.all()
        return context


class SLICPrimerDetailView(DetailView):
    model = SLICPrimer
    context_object_name = 'slic_primer'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        return context


class SLICListView(ListView):
    model = SLICAssembly
    context_object_name = 'slics'

    def get_queryset(self):
        qs = super().get_queryset()
        return qs.all() \
            .only('title', 'date_created') \
            .order_by('-date_created')
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'SLIC Assemblies'
        return context

# Assembly Bundles
def assembly_bundle(request):
    if request.method == 'POST':
        bundle_form = BundleForm(request.POST, request.FILES)

        if bundle_form.is_valid():
            bundle_pk = bundle_create_service(bundle_form.cleaned_data)
            return redirect('bundle-detail', bundle_pk)
    else:
        bundle_form = BundleForm()
    return render(request, 'assembly/assemblybundle_form.html', {'bundle_form': bundle_form})


class AssemblyBundleDetailView(DetailView):
    model = AssemblyBundle
    context_object_name = 'assembly_bundle'

    def get_context_data(self, **kwargs):

        context = super().get_context_data(**kwargs)
        # TODO change how the queries happen here, maybe change to code in the GET request
        context['title'] = self.object.title 
        
        context['gibson'] = self.object.gibson.filter().first()
        if context['gibson']:
            context['gibson_solutions'] = context['gibson'].gibsonsolution_set.all()
        
        context['goldengate'] = self.object.goldengate.filter().first()
        if context['goldengate']:
            context['goldengate_solutions'] = context['goldengate'].goldengatesolution_set.all()

        context['slic'] = self.object.slic.filter().first()
        if context['slic']:
            context['slic_solutions'] = context['slic'].slicsolution_set.all()
        
        context['pcr'] = self.object.pcr.filter().first()
        if context['pcr']:
            context['pcr_solutions'] = context['pcr'].pcrsolution_set.all()

        context['biobrick'] = self.object.biobricks.filter().first()
        if context['biobrick']:
            context['biobrick_solutions'] = context['biobrick'].biobrickssolution_set.all()

        return context


class AssemblyBundleListView(ListView):
    model = AssemblyBundle
    context_object_name = 'bundles'

    def get_queryset(self):
        qs = super().get_queryset()
        return qs.all() \
            .only('title', 'date_created', 'description') \
            .order_by('-date_created')
    
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'Assembly Bundles'
        return context

def assemblies_list(request):
    gibsons = GibsonAssembly.objects.all()
    goldengates = GoldenGateAssembly.objects.all()
    pcrs = PCRAssembly.objects.all()
    slics = SLICAssembly.objects.all()
    biobricks = BioBricksAssembly.objects.all()

    assemblies = sorted(
        chain(gibsons, goldengates, pcrs, slics, biobricks),
        key=lambda assembly: assembly.date_created, reverse=True
    )

    return render(request, 'assembly/assemblies_list.html', { 'assemblies': assemblies})
    
