from django.db.models.deletion import Collector
from django.shortcuts import render,redirect
from django.views.generic.edit import CreateView
from .forms import GibsonForm
from .models import (
    GibsonAssembly, 
    GibsonPart,
    GibsonPrimer,
    GoldenGateAssembly,
    GoldenGatePart,
    GoldenGatePrimer,
    GibsonSolution
)
from django.views.generic import (
    DetailView,
    CreateView
)
from django.contrib.messages.views import SuccessMessageMixin 

from assemblies.gibson import GibsonAssembler
from assemblies.goldengate import GoldenGateAssembler
from time import sleep

def db_list(addgene, igem, dnasu):
    l = []
    if addgene:
        l.append('addgene')
    if igem: 
        l.append('igem')
    if dnasu:
        l.append('dnasu')
    return l

def home(request):
    return render(request, 'assembly/home.html')

def about(request):
    return render(request, 'assembly/about.html', {'title': 'Assembly About'}) 


class GibsonDetailView(DetailView):
    model = GibsonAssembly
    context_object_name = 'gibson'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # TODO change how the queries happen here, maybe change to code in the GET request
        context['title'] = self.object.title 
        context['parts'] = self.object.gibsonpart_set.all()
        context['primer_sets'] =  [part.gibsonprimer_set.all() for part in context['parts']]
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
        'overlap'
    ]

    def form_valid(self, form):
        self.object = form.save()
        gib_assembler = GibsonAssembler(
                            self.object.mv_conc, 
                            self.object.dv_conc, 
                            self.object.dna_conc,
                            self.object.dntp_conc, 
                            self.object.tm, 
                            self.object.backbone_file.path, 
                            self.object.insert_file.path, 
                            db_list(self.object.addgene, self.object.igem, self.object.dnasu), 
                            min_frag=self.object.min_blast, 
                            max_frag=self.object.max_blast, 
                            min_synth=self.object.min_synth, 
                            max_synth=self.object.max_synth,
                            overlap=self.object.overlap
        )
        results, error = gib_assembler.query()
        gib_assembler.solution_building(results)
        gib_assembly, gib_fragments = gib_assembler.design(solution=0)
        # save assembly parts with meta/annotations and their primers here
        # TODO update to have a match % and BLAST solution sequence
        # TODO add a foreach solution in for the assembly
        gibson_solution = GibsonSolution(
            name='New Solution',
            backbone=gib_assembler.backbone.seq,
            query=gib_assembler.query_record.seq,
            solution='',
            parts_count=len(gib_fragments),
            primers_count=len(gib_fragments) * 2,
            match=0.0,
            assembly=self.object
        )
        gibson_solution.save()

        for i, part in enumerate(gib_assembly):
            gibson_part_entry = GibsonPart(
                name=f'{part.name}-{i}',
                database=part.annotations['db'],
                length=part.template.seq.length, 
                length_extended=part.seq.length,
                seq=part.template.seq,
                seq_extended=part.seq,
                position=i,
                solution=gibson_solution            
            )
            gibson_part_entry.save()

            forward_primer = GibsonPrimer(
                name= (part.name + ' forward primer'),
                primer_type='fwd',
                sequence=part.forward_primer.seq,
                footprint=part.forward_primer.footprint,
                tail=part.forward_primer.tail,
                tm_total=0.0,
                tm_footprint=0.0,
                gc=0.0,
                hairpin=False,
                hairpin_tm=0.0,
                hairpin_dg=0.0,
                hairpin_dh=0.0,
                hairpin_ds=0.0,
                homodimer=False,
                homodimer_tm=0.0,
                homodimer_dg=0.0,
                homodimer_dh=0.0,
                homodimer_ds=0.0,
                part=gibson_part_entry
            )
            forward_primer.save()

            reverse_primer = GibsonPrimer(
                name= (part.name + ' reverse primer'),
                primer_type='fwd',
                sequence=part.reverse_primer.seq,
                footprint=part.reverse_primer.footprint,
                tail=part.reverse_primer.tail,
                tm_total=0.0,
                tm_footprint=0.0,
                gc=0.0,
                hairpin=False,
                hairpin_tm=0.0,
                hairpin_dg=0.0,
                hairpin_dh=0.0,
                hairpin_ds=0.0,
                homodimer=False,
                homodimer_tm=0.0,
                homodimer_dg=0.0,
                homodimer_dh=0.0,
                homodimer_ds=0.0,
                part=gibson_part_entry
            )
            reverse_primer.save()

        return super().form_valid(form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'New Gibson Assembly'
        return context


class GibsonPartDetailView(DetailView):
    model = GibsonPart
    context_object_name = 'gibson_part'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        return context


class GibsonPrimerDetailView(DetailView):
    model = GibsonPrimer
    context_object_name = 'gibson_primer'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        return context


class GoldenGateDetailView(DetailView):
    model = GoldenGateAssembly
    context_object_name = 'goldengate'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.title
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
        'overhangs'
    ]

    def form_valid(self, form):
        self.object = form.save()
        gg_assembler = GoldenGateAssembler(
                            self.object.mv_conc, 
                            self.object.dv_conc, 
                            self.object.dna_conc,
                            self.object.dntp_conc, 
                            self.object.tm, 
                            self.object.backbone_file.path, 
                            self.object.insert_file.path, 
                            db_list(self.object.addgene, self.object.igem, self.object.dnasu), 
                            min_frag=self.object.min_blast, 
                            max_frag=self.object.max_blast, 
                            min_synth=self.object.min_synth, 
                            max_synth=self.object.max_synth,
                            ovhngs=self.object.overhangs
        )

        return super().form_valid(form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = 'New Golden Gate Assembly'
        return context


class GoldenGatePartDetailView(DetailView):
    model = GoldenGatePart
    context_object_name = 'goldengate_part'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        return context


class GoldenGatePrimerDetailView(DetailView):
    model = GoldenGatePrimer
    context_object_name = 'goldengate_primer'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        return context


# def submission(request):
#     return render(request, 'assembly/submission.html', {'title': 'Assembly Submission'}) 

# def gibson_form(request):
#     return render(request, 'assembly/gibson_form.html', {'title': 'Gibson Assembly'})
    
# def goldengate_form(request):
#     return render(request, 'assembly/goldengate_form.html', {'title': 'Golden Gate Assembly'})

# def goldengate_dev(request):
#     if request.method == 'POST':
#         form = GoldenGateForm(request.POST)
#         if form.is_valid():
#             return redirect('assembly-submit') 
#     else:
#         form = GoldenGateForm()

#     return render(request, 'assembly/goldengate_dev.html', {'title': 'Golden Gate Assembly', 'form': form})

# def gibson_dev(request):
#     if request.method == 'POST':
#         form = GibsonForm(request.POST)
#         if form.is_valid():
#             return redirect('assembly-submit') 
#     else:
#         form = GibsonForm()

#     return render(request, 'assembly/gibson_dev.html', {'title': 'Gibson Assembly', 'form': form})
