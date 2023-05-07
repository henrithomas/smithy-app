from ..services.gibson import run_gibson
from ..models import (
    GibsonAssembly, 
    GibsonPart,
    GibsonPrimer,
    GibsonSolution
)
from django.views.generic import (
    DetailView,
    CreateView,
    ListView
)
from django.contrib.messages.views import SuccessMessageMixin 

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
        'assembly_ps',
        'mastermix_cost',
        'mastermix_n_reacts'
    ]

    def form_valid(self, form):
        self.object = form.save()
        run_gibson(self.object)
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