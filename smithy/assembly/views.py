from Bio.Blast import Record
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
import os
from dna_features_viewer import (
    GraphicFeature, 
    GraphicRecord, 
    CircularGraphicRecord
)
from django.core.files import File


def db_list(addgene, igem, dnasu):
    l = []
    if addgene:
        l.append('addgene')
    if igem: 
        l.append('igem')
    if dnasu:
        l.append('dnasu')
    return l

def primer_map():
    pass

def part_map(part_model, part, left, right, name, space):
    part_plot_name = f'{name}-map.png'
    temp_plot = f'/home/hthoma/projects/smithy-app/smithy/media/images/{part_plot_name}'

    display_len = int(part.template.seq.length * 0.1)
    seq_len = 2 * display_len + 2 * space + part.template.seq.length
    
    part_start = display_len + space
    part_end = part_start + part.template.seq.length

    right_start = part_end + space
    right_end = right_start + display_len

    fwd_start = part_start - len(part.forward_primer.tail)
    fwd_end = part_start + len(part.forward_primer.footprint)

    rvs_start = part_end - len(part.reverse_primer.tail)
    rvs_end = part_end + len(part.reverse_primer.tail)

    features = [
        GraphicFeature(
            start=0,
            end=display_len, 
            strand=+1, 
            open_left=True, 
            label=f'{left.name}: ({left.annotations["query_start"]}, {left.annotations["query_end"]})',
            color='lightslategrey'
        ),
        GraphicFeature(
            start=part_start, 
            end=part_end, 
            strand=+1,  
            label=f'{part.name}: ({part.annotations["query_start"]}, {part.annotations["query_end"]})',
            color='darkseagreen'
        ),
        GraphicFeature(
            start=right_start, 
            end=right_end, 
            strand=+1,  
            label=f'{right.name}: ({right.annotations["query_start"]}, {right.annotations["query_end"]})',
            open_right=True,
            color='lightslategrey'
        ),
        GraphicFeature(
            thickness=10, 
            start= fwd_start,
            end=fwd_end,
            strand=+1,
            color="#ccccff",
            label=f'{part.name} fwd'
        ),
        GraphicFeature(
            thickness=10, 
            start= rvs_start,
            end=rvs_end,
            strand=-1,
            color="#ccccff",
            label=f'{part.name} rvs'
        )
    ] 

    record = GraphicRecord(sequence_length=seq_len, features=features)
    ax, _ = record.plot(figure_width=10)
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.figure.savefig(temp_plot)

    part_model.part_map.save(part_plot_name, File(open(temp_plot, 'rb')))
    part_model.save()
    os.remove(temp_plot)

def plasmid_map(solution_model, assembly, assembly_name, space, total_len):
    plot_name = f'{assembly_name}_map.png'
    temp_plot = f'/home/hthoma/projects/smithy-app/smithy/media/images/{plot_name}'   
    features = []
    part_start = 0
    part_end = 0

    # set features and colorings for the insert
    for i, part in enumerate(assembly[:-1]):
        part_end = part_start + part.template.seq.length
        if i % 2 == 0:
            part_color = 'darkseagreen'
        else:
            part_color = 'lightslategrey'
        features.append(
            GraphicFeature(
                start=part_start,
                end=part_end,
                strand=+1,
                label=part.name,
                color=part_color
            )
        )
        part_start = part_end + space

    # set the feature and coloring for the backbone
    part_end = part_start + assembly[-1].template.seq.length
    features.append(
        GraphicFeature(
            start=part_start,
            end=part_end,
            strand=+1,
            label=assembly[-1].name,
            color='lightblue'
        )
    )

    record = CircularGraphicRecord(sequence_length=total_len, features=features)
    ax, _ = record.plot(figure_width=10)
    ax.figure.savefig(temp_plot)

    solution_model.plasmid_map.save(plot_name, File(open(temp_plot, 'rb')))
    solution_model.save()
    os.remove(temp_plot)

def gibson_create_service(gibson_obj):
    gib_assembler = GibsonAssembler(
                        gibson_obj.mv_conc, 
                        gibson_obj.dv_conc, 
                        gibson_obj.dna_conc,
                        gibson_obj.dntp_conc, 
                        gibson_obj.tm, 
                        gibson_obj.backbone_file.path, 
                        gibson_obj.insert_file.path, 
                        db_list(gibson_obj.addgene, gibson_obj.igem, gibson_obj.dnasu), 
                        min_frag=gibson_obj.min_blast, 
                        max_frag=gibson_obj.max_blast, 
                        min_synth=gibson_obj.min_synth, 
                        max_synth=gibson_obj.max_synth,
                        overlap=gibson_obj.overlap,
                        multi_query=gibson_obj.multi_query
    )

    if gibson_obj.multi_query:
        results, error = gib_assembler.run_multi_query()
        gib_assembler.multi_query_solution_building(results)
    else:
        results, error = gib_assembler.query()
        gib_assembler.solution_building(results)
    gib_assembly, gib_fragments = gib_assembler.design(solution=0)

    total_len = gib_assembler.backbone.seq.length + gib_assembler.query_record.seq.length

    # save assembly parts with meta/annotations and their primers here
    # TODO update to have a match % and BLAST solution sequence
    # TODO add a foreach solution in for the assembly
    gibson_solution = GibsonSolution(
        name=f'Solution - {gibson_obj.title}',
        backbone=gib_assembler.backbone.seq,
        query=gib_assembler.query_record.seq,
        solution='',
        parts_count=len(gib_fragments),
        primers_count=len(gib_fragments) * 2,
        match=0.0,
        assembly=gibson_obj
    )
    gibson_solution.save()

    plasmid_map(gibson_solution, gib_assembly, gibson_obj.title, 0, total_len)


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

        if i == 0:
            left_index = len(gib_assembly) - 1
        else:
            left_index = i - 1

        if i == len(gib_assembly) - 1:
            right_index = 0
        else:
            right_index = i + 1

        part_map(
            gibson_part_entry, 
            part, 
            gib_assembly[left_index], 
            gib_assembly[right_index], 
            f'{part.name}-{i}', 
            0
        )

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
    pass

def goldengate_create_service(goldengate_obj):
    pass
    gg_assembler = GoldenGateAssembler(
                        goldengate_obj.mv_conc, 
                        goldengate_obj.dv_conc, 
                        goldengate_obj.dna_conc,
                        goldengate_obj.dntp_conc, 
                        goldengate_obj.tm, 
                        goldengate_obj.backbone_file.path, 
                        goldengate_obj.insert_file.path, 
                        db_list(goldengate_obj.addgene,goldengate_obj.igem, goldengate_obj.dnasu), 
                        min_frag=goldengate_obj.min_blast, 
                        max_frag=goldengate_obj.max_blast, 
                        min_synth=goldengate_obj.min_synth, 
                        max_synth=goldengate_obj.max_synth,
                        ovhngs=goldengate_obj.overhangs,
                        multi_query=goldengate_obj.multi_query
    )
    
    if goldengate_obj.multi_query:
        results, error = gg_assembler.run_multi_query()
        gg_assembler.multi_query_solution_building(results)
    else:
        results, error = gg_assembler.query()
        gg_assembler.solution_building(results)
    gg_assembly, gg_fragments = gg_assembler.design(solution=0)

    total_len = gg_assembler.backbone.seq.length + gg_assembler.query_record.seq.length + 4

    # TODO update to have a match % and BLAST solution sequence
    # TODO add a foreach solution in for the assembly
    goldengate_solution = GoldenGateSolution(
        name=f'Solution - {goldengate_obj.title}',
        backbone=gg_assembler.backbone.seq,
        query=gg_assembler.query_record.seq,
        solution='',
        parts_count=len(gg_fragments),
        primers_count=len(gg_fragments) * 2,
        match=0.0,
        assembly=goldengate_obj
    )
    goldengate_solution.save()

    plasmid_map(goldengate_solution, gg_assembly, goldengate_obj.title, 4, total_len)

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

        if i == 0:
            left_index = len(gg_assembly) - 1
        else:
            left_index = i - 1

        if i == len(gg_assembly) - 1:
            right_index = 0
        else:
            right_index = i + 1

        part_map(
            goldengate_part_entry, 
            part, 
            gg_assembly[left_index], 
            gg_assembly[right_index], 
            f'{part.name}-{i}', 
            4
        )

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
        'overhangs',
        'multi_query'
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
        return context


class GoldenGatePrimerDetailView(DetailView):
    model = GoldenGatePrimer
    context_object_name = 'goldengate_primer'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['title'] = self.object.name
        return context

