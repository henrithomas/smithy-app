from ..service import pcr_create_service
from ..models import (
    PCRAssembly,
    PCRSolution,
    PCRPart,
    PCRPrimer
)
from django.views.generic import (
    DetailView,
    CreateView,
    ListView
)
from django.contrib.messages.views import SuccessMessageMixin 

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
        'parts_pref',
        'mastermix_cost',
        'mastermix_n_reacts'
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