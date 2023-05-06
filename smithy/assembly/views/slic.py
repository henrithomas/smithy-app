from ..service import slic_create_service
from ..models import (
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
        'parts_pref',
        'chewback_ps',
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