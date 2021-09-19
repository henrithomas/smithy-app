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
    CreateView
)
from django.contrib.messages.views import SuccessMessageMixin 
import os
from .forms import BundleForm


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
        'multi_query'
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
        'scarless'
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
        'multi_query'
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
        'multi_query'
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
        'multi_query'
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


def assembly_bundle(request):
    if request.method == 'POST':
        bundle_form = BundleForm(request.POST, request.FILES)

        if bundle_form.is_valid():
            # TODO add bundle service and assembly building here
            bundle_pk = bundle_create_service(bundle_form.cleaned_data)
            redirect('bundle-detail', bundle_pk)
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
        # context['solutions'] = self.object.gibsonsolution_set.all()
        return context