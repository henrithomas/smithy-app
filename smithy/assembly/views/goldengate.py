from ..services.goldengate import run_goldengate
from ..models import (
    GoldenGateAssembly,
    GoldenGatePart,
    GoldenGatePrimer,
    GoldenGateSolution
)
from django.views.generic import (
    DetailView,
    CreateView,
    ListView
)
from django.contrib.messages.views import SuccessMessageMixin

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
        'assembly_ps',
        'mastermix_cost',
        'mastermix_n_reacts'
    ]

    def form_valid(self, form):
        self.object = form.save()
        run_goldengate(self.object)
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