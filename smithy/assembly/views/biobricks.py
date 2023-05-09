from ..services.biobricks import run_biobricks
from ..models import (
    BioBricksAssembly,
    BioBricksPart,
    BioBricksPrimer,
    BioBricksSolution
)
from django.views.generic import (
    DetailView,
    CreateView,
    ListView
)
from django.contrib.messages.views import SuccessMessageMixin 

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
        run_biobricks(self.object)
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