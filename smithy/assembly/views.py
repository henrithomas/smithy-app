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
    GibsonSolution,
    GoldenGateSolution
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
        gib_assembly, gib_fragments = gib_assembler.design(solution=3)
        # save assembly parts with meta/annotations and their primers here
        # TODO update to have a match % and BLAST solution sequence
        # TODO add a foreach solution in for the assembly
        gibson_solution = GibsonSolution(
            name=f'Solution - {self.object.title}',
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
                name=part.name,
                database=part.annotations['db'],
                length=part.template.seq.length, 
                length_extended=part.seq.length,
                seq=part.template.seq,
                seq_extended=part.seq,
                position=i,
                solution=gibson_solution,
                query_start = part.annotations['query_start'],
                query_end = part.annotations['query_end'],
                subject_start = part.annotations['subject_start'],
                subject_end = part.annotations['subject_end'] 
            )
            gibson_part_entry.save()

            forward_primer = GibsonPrimer(
                name= f'{gibson_part_entry.name} forward primer',
                primer_type='fwd',
                sequence=part.forward_primer.seq,
                footprint=part.forward_primer.footprint,
                tail=part.forward_primer.tail,
                tm_total=part.annotations['forward_primer']['tm_total'],
                tm_footprint=part.annotations['forward_primer']['tm_footprint'],
                gc=part.annotations['forward_primer']['gc'],
                hairpin=part.annotations['forward_primer']['hairpin'],
                hairpin_tm=part.annotations['forward_primer']['hairpin_tm'],
                hairpin_dg=part.annotations['forward_primer']['hairpin_dg'],
                hairpin_dh=part.annotations['forward_primer']['hairpin_dh'],
                hairpin_ds=part.annotations['forward_primer']['hairpin_ds'],
                homodimer=part.annotations['forward_primer']['homodimer'],
                homodimer_tm=part.annotations['forward_primer']['homodimer_tm'],
                homodimer_dg=part.annotations['forward_primer']['homodimer_dg'],
                homodimer_dh=part.annotations['forward_primer']['homodimer_dh'],
                homodimer_ds=part.annotations['forward_primer']['homodimer_ds'],
                part=gibson_part_entry
            )
            forward_primer.save()

            reverse_primer = GibsonPrimer(
                name= f'{gibson_part_entry.name} reverse primer ',
                primer_type='rvs',
                sequence=part.reverse_primer.seq,
                footprint=part.reverse_primer.footprint,
                tail=part.reverse_primer.tail,
                tm_total=part.annotations['reverse_primer']['tm_total'],
                tm_footprint=part.annotations['reverse_primer']['tm_footprint'],
                gc=part.annotations['reverse_primer']['gc'],
                hairpin=part.annotations['reverse_primer']['hairpin'],
                hairpin_tm=part.annotations['reverse_primer']['hairpin_tm'],
                hairpin_dg=part.annotations['reverse_primer']['hairpin_dg'],
                hairpin_dh=part.annotations['reverse_primer']['hairpin_dh'],
                hairpin_ds=part.annotations['reverse_primer']['hairpin_ds'],
                homodimer=part.annotations['reverse_primer']['homodimer'],
                homodimer_tm=part.annotations['reverse_primer']['homodimer_tm'],
                homodimer_dg=part.annotations['reverse_primer']['homodimer_dg'],
                homodimer_dh=part.annotations['reverse_primer']['homodimer_dh'],
                homodimer_ds=part.annotations['reverse_primer']['homodimer_ds'],
                part=gibson_part_entry 
            )
            reverse_primer.save()

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
        
        results, error = gg_assembler.query()
        gg_assembler.solution_building(results)
        gg_assembly, gg_fragments = gg_assembler.design(solution=3)

        # TODO update to have a match % and BLAST solution sequence
        # TODO add a foreach solution in for the assembly
        goldengate_solution = GoldenGateSolution(
            name=f'Solution - {self.object.title}',
            backbone=gg_assembler.backbone.seq,
            query=gg_assembler.query_record.seq,
            solution='',
            parts_count=len(gg_fragments),
            primers_count=len(gg_fragments) * 2,
            match=0.0,
            assembly=self.object
        )
        goldengate_solution.save()

        for i, part in enumerate(gg_assembly):
            goldengate_part_entry = GoldenGatePart(
                name=part.name,
                database=part.annotations['db'],
                length=part.template.seq.length, 
                length_extended=part.seq.length,
                seq=part.template.seq,
                seq_extended=part.seq,
                position=i,
                solution=goldengate_solution,
                query_start = part.annotations['query_start'],
                query_end = part.annotations['query_end'],
                subject_start = part.annotations['subject_start'],
                subject_end = part.annotations['subject_end']             
            )
            goldengate_part_entry.save()

            forward_primer = GoldenGatePrimer(
                name= f'{goldengate_part_entry.name} forward primer',
                primer_type='fwd',
                sequence=part.forward_primer.seq,
                footprint=part.forward_primer.footprint,
                tail=part.forward_primer.tail,
                tm_total=part.annotations['forward_primer']['tm_total'],
                tm_footprint=part.annotations['forward_primer']['tm_footprint'],
                gc=part.annotations['forward_primer']['gc'],
                hairpin=part.annotations['forward_primer']['hairpin'],
                hairpin_tm=part.annotations['forward_primer']['hairpin_tm'],
                hairpin_dg=part.annotations['forward_primer']['hairpin_dg'],
                hairpin_dh=part.annotations['forward_primer']['hairpin_dh'],
                hairpin_ds=part.annotations['forward_primer']['hairpin_ds'],
                homodimer=part.annotations['forward_primer']['homodimer'],
                homodimer_tm=part.annotations['forward_primer']['homodimer_tm'],
                homodimer_dg=part.annotations['forward_primer']['homodimer_dg'],
                homodimer_dh=part.annotations['forward_primer']['homodimer_dh'],
                homodimer_ds=part.annotations['forward_primer']['homodimer_ds'],
                part=goldengate_part_entry
            )
            forward_primer.save()

            reverse_primer = GoldenGatePrimer(
                name= f'{goldengate_part_entry.name} reverse primer ',
                primer_type='rvs',
                sequence=part.reverse_primer.seq,
                footprint=part.reverse_primer.footprint,
                tail=part.reverse_primer.tail,
                tm_total=part.annotations['reverse_primer']['tm_total'],
                tm_footprint=part.annotations['reverse_primer']['tm_footprint'],
                gc=part.annotations['reverse_primer']['gc'],
                hairpin=part.annotations['reverse_primer']['hairpin'],
                hairpin_tm=part.annotations['reverse_primer']['hairpin_tm'],
                hairpin_dg=part.annotations['reverse_primer']['hairpin_dg'],
                hairpin_dh=part.annotations['reverse_primer']['hairpin_dh'],
                hairpin_ds=part.annotations['reverse_primer']['hairpin_ds'],
                homodimer=part.annotations['reverse_primer']['homodimer'],
                homodimer_tm=part.annotations['reverse_primer']['homodimer_tm'],
                homodimer_dg=part.annotations['reverse_primer']['homodimer_dg'],
                homodimer_dh=part.annotations['reverse_primer']['homodimer_dh'],
                homodimer_ds=part.annotations['reverse_primer']['homodimer_ds'],
                part=goldengate_part_entry 
            )
            reverse_primer.save()

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
        return context


class GoldenGatePrimerDetailView(DetailView):
    model = GoldenGatePrimer
    context_object_name = 'goldengate_primer'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        return context

