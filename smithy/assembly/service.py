from django.db.models.deletion import Collector
from tests import cut_locations
from .models import (
    AssemblyBundle,
    GibsonAssembly, 
    GibsonPart,
    GibsonPrimer,
    GoldenGateAssembly,
    GoldenGatePart,
    GoldenGatePrimer,
    GibsonSolution,
    GoldenGateSolution,
    PCRAssembly,
    PCRPart,
    PCRPrimer,
    PCRSolution,
    SLICAssembly,
    SLICPart,
    SLICPrimer,
    SLICSolution,
    BioBricksAssembly,
    BioBricksPart,
    BioBricksPrimer,
    BioBricksSolution
)
from assemblies.gibson import (
    GibsonAssembler,
    SLICAssembler,
    PCRAssembler
)
from assemblies.goldengate import GoldenGateAssembler
from assemblies.traditional_re import BioBrickAssembler
from assemblies.assembler import Assembler
from dna_features_viewer import (
    GraphicFeature, 
    GraphicRecord, 
    CircularGraphicRecord
)
from django.core.files import File
import os
import json
from datetime import datetime
import csv
import matplotlib.pyplot as plt
from statistics import mean


def parts_csv(solution_model, parts):
    """
    Exports the assembly parts data to the chosen file in csv format.


    Parameters
    ----------
    solution_model : AssemblySolution

    parts : list
        List of Amplicon objects of the assembly
        


    Returns
    -------
    None
    """
    file_name = f'{solution_model.name.replace(" ", "-")}-{solution_model.pk}-parts.csv'
    temp_file = f'/home/hthoma/projects/smithy-app/smithy/media/csv/{file_name}'

    fields = ['id', 'length', 'length_ext', 'seq', 'seq_ext', 'query_start', 'query_end', 'subject_start', 'subject_end']
    csv_list = [
        {
            'id': part.name, 
            'length': part.template.seq.length, 
            'length_ext': part.seq.length, 
            'seq': part.template.seq.watson, 
            'seq_ext': part.seq.watson,
            'query_start': part.annotations['query_start'],
            'query_end': part.annotations['query_end'],
            'subject_start': part.annotations['subject_start'],
            'subject_end': part.annotations['subject_end']
        }
        for part in parts
    ]
    
    with open(temp_file, 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fields, restval='NONE')
        writer.writeheader()
        writer.writerows(csv_list)

    solution_model.parts_file.save(file_name, File(open(temp_file, newline='')))
    solution_model.save()
    os.remove(temp_file)

def primers_csv(solution_model, parts):
    """
    Exports the assembly primers to the chosen file in csv format.


    Parameters
    ----------
    solution_model : AssemblySolution
    
    parts : list
        List of Amplicon objects of the assembly



    Returns
    -------        
    None
    """
    file_name = f'{solution_model.name.replace(" ", "-")}-{solution_model.pk}-primers.csv'
    temp_file = f'/home/hthoma/projects/smithy-app/smithy/media/csv/{file_name}'

    fields = ['id', 'primer_type', 'sequence', 'footprint', 'tail', 'tm_footprint', 'tm_total', 'gc', 
              'hairpin', 'hairpin_tm', 'hairpin_dg', 'hairpin_dh', 'hairpin_ds',
              'homodimer', 'homodimer_tm', 'homodimer_dg', 'homodimer_dh', 'homodimer_ds'
    ]
    csv_list = []

    for part in parts:
        csv_list.extend([{
                'id': f'{part.name}-fwd',
                'primer_type': 'fwd',
                'sequence': part.forward_primer.seq,
                'footprint': part.forward_primer.footprint,
                'tail': part.forward_primer.tail,
                'tm_footprint': part.annotations['forward_primer']['tm_footprint'], 
                'tm_total': part.annotations['forward_primer']['tm_total'], 
                'gc': part.annotations['forward_primer']['gc'], 
                'hairpin': part.annotations['forward_primer']['hairpin'], 
                'hairpin_tm': part.annotations['forward_primer']['hairpin_tm'], 
                'hairpin_dg': part.annotations['forward_primer']['hairpin_dg'], 
                'hairpin_dh': part.annotations['forward_primer']['hairpin_dh'], 
                'hairpin_ds': part.annotations['forward_primer']['hairpin_ds'],
                'homodimer': part.annotations['forward_primer']['homodimer'], 
                'homodimer_tm': part.annotations['forward_primer']['homodimer_tm'], 
                'homodimer_dg': part.annotations['forward_primer']['homodimer_dg'], 
                'homodimer_dh': part.annotations['forward_primer']['homodimer_dh'], 
                'homodimer_ds': part.annotations['forward_primer']['homodimer_ds'],
            },
            {
                'id': f'{part.name}-rvs',
                'primer_type': 'rvs',
                'sequence': part.reverse_primer.seq,
                'footprint': part.reverse_primer.footprint,
                'tail': part.reverse_primer.tail,
                'tm_footprint': part.annotations['reverse_primer']['tm_footprint'], 
                'tm_total': part.annotations['reverse_primer']['tm_total'], 
                'gc': part.annotations['reverse_primer']['gc'], 
                'hairpin': part.annotations['reverse_primer']['hairpin'], 
                'hairpin_tm': part.annotations['reverse_primer']['hairpin_tm'], 
                'hairpin_dg': part.annotations['reverse_primer']['hairpin_dg'], 
                'hairpin_dh': part.annotations['reverse_primer']['hairpin_dh'], 
                'hairpin_ds': part.annotations['reverse_primer']['hairpin_ds'],
                'homodimer': part.annotations['reverse_primer']['homodimer'], 
                'homodimer_tm': part.annotations['reverse_primer']['homodimer_tm'], 
                'homodimer_dg': part.annotations['reverse_primer']['homodimer_dg'], 
                'homodimer_dh': part.annotations['reverse_primer']['homodimer_dh'], 
                'homodimer_ds': part.annotations['reverse_primer']['homodimer_ds'],
            }])

    with open(temp_file, 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fields, restval='NONE')
        writer.writeheader()
        writer.writerows(csv_list)

    solution_model.primers_file.save(file_name, File(open(temp_file, newline='')))
    solution_model.save()
    os.remove(temp_file)

def solution_analysis(assembly, fragments, query_length):
    primer_lengths = []
    primer_tms = []
    part_lengths = [
        part.template.seq.length
        for part in assembly[:-1]
    ]
    for part in assembly: 
        primer_lengths.extend(
            [
                len(part.forward_primer.seq),
                len(part.reverse_primer.seq)
            ]
        )
        primer_tms.extend(
            [
                part.annotations['forward_primer']['tm_footprint'],
                part.annotations['reverse_primer']['tm_footprint']
            ]
        )
    score_sum = sum([fragment.score for fragment in fragments])
    match_p = round(
        score_sum / query_length,
        2
    ) * 100
    synth_p = 100.00 - match_p
    part_ave = int(mean(part_lengths))
    primer_ave = int(mean(primer_lengths))
    primer_tm_ave = round(mean(primer_tms), 2)
    part_max = max(part_lengths)
    part_min = min(part_lengths)
    pass
    return match_p, synth_p, part_ave, primer_ave, primer_tm_ave, part_max, part_min    

def write_uploaded_file(f, fname):
    with open(fname, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)

def db_list(addgene, igem, dnasu):
    """
    


    Parameters
    ----------



    Returns
    -------
    """
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
    """
    


    Parameters
    ----------



    Returns
    -------
    """
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
    plt.close(ax.get_figure())

    part_model.part_map.save(part_plot_name, File(open(temp_plot, 'rb')))
    part_model.save()
    os.remove(temp_plot)

def plasmid_map(solution_model, assembly, assembly_name, space, total_len):
    """
    


    Parameters
    ----------



    Returns
    -------
    """
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
    plt.close(ax.get_figure())

    solution_model.plasmid_map.save(plot_name, File(open(temp_plot, 'rb')))
    solution_model.save()
    os.remove(temp_plot)

def gibson_create_service(obj):
    """
    


    Parameters
    ----------



    Returns
    -------
    """
    assembler = GibsonAssembler(
        obj.mv_conc, 
        obj.dv_conc, 
        obj.dna_conc,
        obj.dntp_conc, 
        obj.tm, 
        obj.backbone_file.path, 
        obj.insert_file.path, 
        db_list(obj.addgene, obj.igem, obj.dnasu), 
        min_frag=obj.min_blast, 
        max_frag=obj.max_blast, 
        min_synth=obj.min_synth, 
        max_synth=obj.max_synth,
        overlap=obj.overlap,
        multi_query=obj.multi_query
    )

    if obj.multi_query:
        results, error = assembler.run_multi_query()
        assembler.multi_query_solution_building(results)
    else:
        results, error = assembler.query()
        assembler.solution_building(results)
    assembly, fragments = assembler.design(solution=0)

    gibson_solution_service(obj, assembler, assembly, fragments)

def goldengate_create_service(obj):
    """
    


    Parameters
    ----------



    Returns
    -------
    """
    assembler = GoldenGateAssembler(
        obj.mv_conc, 
        obj.dv_conc, 
        obj.dna_conc,
        obj.dntp_conc, 
        obj.tm, 
        obj.backbone_file.path, 
        obj.insert_file.path, 
        db_list(obj.addgene,obj.igem, obj.dnasu), 
        min_frag=obj.min_blast, 
        max_frag=obj.max_blast, 
        min_synth=obj.min_synth, 
        max_synth=obj.max_synth,
        ovhngs=obj.overhangs,
        multi_query=obj.multi_query,
        scarless=obj.scarless
    )

    if obj.multi_query:
        results, error = assembler.run_multi_query()
        assembler.multi_query_solution_building(results)
    else:
        results, error = assembler.query()
        assembler.solution_building(results)
    assembly, fragments = assembler.design(solution=0)

    goldengate_solution_service(obj, assembler, assembly, fragments)

def biobricks_create_service(obj):
    """
    


    Parameters
    ----------



    Returns
    -------
    """
    assembler = BioBrickAssembler(
        obj.mv_conc, 
        obj.dv_conc, 
        obj.dna_conc,
        obj.dntp_conc, 
        obj.tm, 
        obj.backbone_file.path, 
        obj.insert_file.path, 
        db_list(obj.addgene, obj.igem, obj.dnasu), 
        min_frag=obj.min_blast, 
        max_frag=obj.max_blast, 
        min_synth=obj.min_synth, 
        max_synth=obj.max_synth,
        multi_query=obj.multi_query  
    )

    if obj.multi_query:
        results, error = assembler.run_multi_query()
        assembler.multi_query_solution_building(results)
    else:
        results, error = assembler.query()
        assembler.solution_building(results)
    assembly, fragments = assembler.design(solution=0)

    biobricks_solution_service(obj, assembler, assembly, fragments)   

def pcr_create_service(obj):
    """
    


    Parameters
    ----------



    Returns
    -------
    """
    assembler = PCRAssembler(
        obj.mv_conc,
        obj.dv_conc,
        obj.dna_conc,
        obj.dntp_conc,
        obj.tm,
        obj.backbone_file.path,
        obj.insert_file.path,
        db_list(obj.addgene, obj.igem, obj.dnasu),
        min_frag=obj.min_blast,
        max_frag=obj.max_blast,
        min_synth=obj.min_synth,
        max_synth=obj.max_synth,
        overlap=obj.overlap,
        multi_query=obj.multi_query
    )

    if obj.multi_query:
        results, error = assembler.run_multi_query()
        assembler.multi_query_solution_building(results)
    else:
        results, error = assembler.query()
        assembler.solution_building(results)
    assembly, fragments = assembler.design(solution=0)

    pcr_solution_service(obj, assembler, assembly, fragments)

def slic_create_service(obj):
    """
    


    Parameters
    ----------



    Returns
    -------
    """
    assembler = SLICAssembler(
        obj.mv_conc,
        obj.dv_conc,
        obj.dna_conc,
        obj.dntp_conc,
        obj.tm,
        obj.backbone_file.path,
        obj.insert_file.path,
        db_list(obj.addgene, obj.igem, obj.dnasu),
        min_frag=obj.min_blast,
        max_frag=obj.max_blast,
        min_synth=obj.min_synth,
        max_synth=obj.max_synth,
        overlap=obj.overlap,
        multi_query=obj.multi_query
    )

    if obj.multi_query:
        results, error = assembler.run_multi_query()
        assembler.multi_query_solution_building(results)
    else:
        results, error = assembler.query()
        assembler.solution_building(results)
    assembly, fragments = assembler.design(solution=0)

    slic_solution_service(obj, assembler, assembly, fragments)    

def gibson_solution_service(obj, assembler, assembly, fragments):
    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length
    primer_lengths = []
    primer_tms = []
    part_lengths = [
        part.template.seq.length
        for part in assembly[:-1]
    ]
    for part in assembly: 
        primer_lengths.extend(
            [
                len(part.forward_primer.seq),
                len(part.reverse_primer.seq)
            ]
        )
        primer_tms.extend(
            [
                part.annotations['forward_primer']['tm_footprint'],
                part.annotations['reverse_primer']['tm_footprint']
            ]
        )
    score_sum = sum([fragment.score for fragment in fragments])
    match_percentage = round(
        score_sum / assembler.query_record.seq.length,
        2
    ) * 100
    synth_percentage = 100.00 - match_percentage
    part_ave = int(mean(part_lengths))
    primer_ave = int(mean(primer_lengths))
    primer_tm_ave = round(mean(primer_tms), 2)
    part_max = max(part_lengths)
    part_min = min(part_lengths)

    # save assembly parts with meta/annotations and their primers here
    # TODO update to have a match % and BLAST solution sequence
    # TODO add a foreach solution in for the assembly
    gibson_solution = GibsonSolution(
        name=f'{obj.title} Solution',
        backbone=assembler.backbone.seq,
        query=assembler.query_record.seq,
        solution='',
        parts_count=len(fragments),
        primers_count=len(fragments) * 2,
        match=match_percentage,
        synth_amount=synth_percentage,
        re_enzymes=False,
        part_length_average=part_ave,
        primer_length_average=primer_ave,
        longest_part=part_max,
        shortest_part=part_min,
        tm_average=primer_tm_ave,
        solution_length=assembler.query_record.seq.length,
        assembly=obj
    )
    gibson_solution.save()

    plasmid_map(gibson_solution, assembly, obj.title, 0, total_len)
    parts_csv(gibson_solution, assembly)
    primers_csv(gibson_solution, assembly)

    for i, part in enumerate(assembly):
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
            left_index = len(assembly) - 1
        else:
            left_index = i - 1

        if i == len(assembly) - 1:
            right_index = 0
        else:
            right_index = i + 1

        part_map(
            gibson_part_entry, 
            part, 
            assembly[left_index], 
            assembly[right_index], 
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

def goldengate_solution_service(obj, assembler, assembly, fragments):
    space = 0 if obj.scarless else 4
    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length + space

    # TODO update to have a match % and BLAST solution sequence
    # TODO add a foreach solution in for the assembly
    goldengate_solution = GoldenGateSolution(
        name=f'{obj.title} Solution',
        backbone=assembler.backbone.seq,
        query=assembler.query_record.seq,
        solution='',
        parts_count=len(fragments),
        primers_count=len(fragments) * 2,
        match=0.0,
        assembly=obj
    )
    goldengate_solution.save()

    plasmid_map(goldengate_solution, assembly, obj.title, space, total_len)
    parts_csv(goldengate_solution, assembly)
    primers_csv(goldengate_solution, assembly)

    for i, part in enumerate(assembly):
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
            subject_end = part.annotations['subject_end'],
            cuts=part.annotations['cuts'],
            cut_locations=json.dumps(part.annotations['cut_locations'])           
        )
        goldengate_part_entry.save()

        if i == 0:
            left_index = len(assembly) - 1
        else:
            left_index = i - 1

        if i == len(assembly) - 1:
            right_index = 0
        else:
            right_index = i + 1

        part_map(
            goldengate_part_entry, 
            part, 
            assembly[left_index], 
            assembly[right_index], 
            f'{part.name}-{i}', 
            space
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

def biobricks_solution_service(obj, assembler, assembly, fragments):
    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length

    biobricks_solution = BioBricksSolution(
        name=f'{obj.title} Solution',
        backbone=assembler.backbone.seq,
        query=assembler.query_record.seq,
        solution='',
        parts_count=len(fragments),
        primers_count=len(fragments) * 2,
        match=0.0,
        assembly=obj
    )
    biobricks_solution.save()

    plasmid_map(biobricks_solution, assembly, obj.title, 0, total_len)
    parts_csv(biobricks_solution, assembly)
    primers_csv(biobricks_solution, assembly)

    for i, part in enumerate(assembly):
        biobricks_part_entry = BioBricksPart(
            name=part.name,
            database=part.annotations['db'],
            length=part.template.seq.length, 
            length_extended=part.seq.length,
            seq=part.template.seq,
            seq_extended=part.seq,
            position=i,
            solution=biobricks_solution,
            query_start = part.annotations['query_start'],
            query_end = part.annotations['query_end'],
            subject_start = part.annotations['subject_start'],
            subject_end = part.annotations['subject_end'] ,
            cuts=0
        )
        biobricks_part_entry.save()

        if i == 0:
            left_index = len(assembly) - 1
        else:
            left_index = i - 1

        if i == len(assembly) - 1:
            right_index = 0
        else:
            right_index = i + 1

        part_map(
            biobricks_part_entry, 
            part, 
            assembly[left_index], 
            assembly[right_index], 
            f'{part.name}-{i}', 
            0
        )

        forward_primer = BioBricksPrimer(
            name= f'{biobricks_part_entry.name} forward primer',
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
            part=biobricks_part_entry
        )
        forward_primer.save()

        reverse_primer = BioBricksPrimer(
            name= f'{biobricks_part_entry.name} reverse primer ',
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
            part=biobricks_part_entry 
        )
        reverse_primer.save()  

def pcr_solution_service(obj, assembler, assembly, fragments):
    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length

    pcr_solution = PCRSolution(
        name=f'{obj.title} Solution',
        backbone=assembler.backbone.seq,
        query=assembler.query_record.seq,
        solution='',
        parts_count=len(fragments),
        primers_count=len(fragments) * 2, 
        match=0.0,
        assembly=obj
    )
    pcr_solution.save()

    plasmid_map(pcr_solution, assembly, obj.title, 0, total_len)
    parts_csv(pcr_solution, assembly)
    primers_csv(pcr_solution, assembly)

    for i, part in enumerate(assembly):
        pcr_part_entry = PCRPart(
            name=part.name,
            database=part.annotations['db'],
            length=part.template.seq.length, 
            length_extended=part.seq.length,
            seq=part.template.seq,
            seq_extended=part.seq,
            position=i,
            solution=pcr_solution,
            query_start = part.annotations['query_start'],
            query_end = part.annotations['query_end'],
            subject_start = part.annotations['subject_start'],
            subject_end = part.annotations['subject_end']
        )
        pcr_part_entry.save()

        if i == 0:
            left_index = len(assembly) - 1
        else:
            left_index = i - 1

        if i == len(assembly) - 1:
            right_index = 0
        else:
            right_index = i + 1

        part_map(
            pcr_part_entry, 
            part, 
            assembly[left_index], 
            assembly[right_index], 
            f'{part.name}-{i}', 
            0
        )

        forward_primer = PCRPrimer(
            name= f'{pcr_part_entry.name} forward primer',
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
            part=pcr_part_entry
        )
        forward_primer.save()

        reverse_primer = PCRPrimer(
            name= f'{pcr_part_entry.name} reverse primer ',
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
            part=pcr_part_entry 
        )
        reverse_primer.save()

def slic_solution_service(obj, assembler, assembly, fragments):
    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length

    slic_solution = SLICSolution(
        name=f'{obj.title} Solution',
        backbone=assembler.backbone.seq,
        query=assembler.query_record.seq,
        solution='',
        parts_count=len(fragments),
        primers_count=len(fragments) * 2, 
        match=0.0,
        assembly=obj
    )
    slic_solution.save()

    plasmid_map(slic_solution, assembly, obj.title, 0, total_len)
    parts_csv(slic_solution, assembly)
    primers_csv(slic_solution, assembly)

    for i, part in enumerate(assembly):
        slic_part_entry = SLICPart(
            name=part.name,
            database=part.annotations['db'],
            length=part.template.seq.length, 
            length_extended=part.seq.length,
            seq=part.template.seq,
            seq_extended=part.seq,
            position=i,
            solution=slic_solution,
            query_start = part.annotations['query_start'],
            query_end = part.annotations['query_end'],
            subject_start = part.annotations['subject_start'],
            subject_end = part.annotations['subject_end']
        )
        slic_part_entry.save()

        if i == 0:
            left_index = len(assembly) - 1
        else:
            left_index = i - 1

        if i == len(assembly) - 1:
            right_index = 0
        else:
            right_index = i + 1

        part_map(
            slic_part_entry, 
            part, 
            assembly[left_index], 
            assembly[right_index], 
            f'{part.name}-{i}', 
            0
        )

        forward_primer = SLICPrimer(
            name= f'{slic_part_entry.name} forward primer',
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
            part=slic_part_entry
        )
        forward_primer.save()

        reverse_primer = SLICPrimer(
            name= f'{slic_part_entry.name} reverse primer ',
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
            part=slic_part_entry 
        )
        reverse_primer.save()

def bundle_create_service(bundle_data):
    # save backbone and insert files here
    # smithy-app/smithy/media
    path = '/home/hthoma/projects/smithy-app/smithy/media/fasta'
    backbone_name, backbone_extension = os.path.splitext(bundle_data['backbone_file'].name) 
    insert_name, insert_extension = os.path.splitext(bundle_data['insert_file'].name)
    backbone_file_path = f'{path}/{backbone_name}_{datetime.now():%S%f}{backbone_extension}'
    insert_file_path = f'{path}/{insert_name}_{datetime.now():%S%f}{insert_extension}'
    
    write_uploaded_file(bundle_data['backbone_file'], backbone_file_path)
    write_uploaded_file(bundle_data['insert_file'], insert_file_path)

    # create base assembler object for running queries and creating
    # a solution tree
    query_assembler = Assembler(
        bundle_data['mv_conc'], 
        bundle_data['dv_conc'], 
        bundle_data['dna_conc'],
        bundle_data['dntp_conc'], 
        bundle_data['tm'], 
        backbone_file_path, 
        insert_file_path, 
        db_list(bundle_data['addgene'], bundle_data['igem'], bundle_data['dnasu']), 
        min_frag=bundle_data['min_blast'], 
        max_frag=bundle_data['max_blast'], 
        min_synth=bundle_data['min_synth'], 
        max_synth=bundle_data['max_synth'],
        multi_query=bundle_data['multi_query']
    )
    if bundle_data['multi_query']:
        results, error = query_assembler.run_multi_query()
        query_assembler.multi_query_solution_building(results)
    else:
        results, error = query_assembler.query()
        query_assembler.solution_building(results)

    bundle = AssemblyBundle(
        title=bundle_data['title'],
        description=bundle_data['description']
    )
    bundle.save()

    if bundle_data['gibson']:
        gibson_obj = GibsonAssembly(
            title=bundle_data['title'],
            multipart=bundle_data['multipart'],
            addgene=bundle_data['addgene'],
            igem=bundle_data['igem'],
            dnasu=bundle_data['dnasu'],
            min_blast=bundle_data['min_blast'],
            max_blast=bundle_data['max_blast'],
            min_synth=bundle_data['min_synth'],
            max_synth=bundle_data['max_synth'],
            mv_conc=bundle_data['mv_conc'],
            dv_conc=bundle_data['dv_conc'],
            dntp_conc=bundle_data['dntp_conc'],
            dna_conc=bundle_data['dna_conc'],
            tm=bundle_data['tm'],
            multi_query=bundle_data['multi_query'],
            overlap=bundle_data['overlap']
        )
        gibson_obj.backbone_file.save(bundle_data['backbone_file'].name, File(open(backbone_file_path, 'rb')))
        gibson_obj.insert_file.save(bundle_data['insert_file'].name, File(open(insert_file_path, 'rb')))
        gibson_obj.save()

        gibson_assembler = GibsonAssembler(
            bundle_data['mv_conc'], 
            bundle_data['dv_conc'], 
            bundle_data['dna_conc'],
            bundle_data['dntp_conc'], 
            bundle_data['tm'], 
            backbone_file_path, 
            insert_file_path, 
            db_list(bundle_data['addgene'], bundle_data['igem'], bundle_data['dnasu']), 
            min_frag=bundle_data['min_blast'], 
            max_frag=bundle_data['max_blast'], 
            min_synth=bundle_data['min_synth'], 
            max_synth=bundle_data['max_synth'],
            overlap=bundle_data['overlap'],
            multi_query=bundle_data['multi_query']
        )
        gibson_assembler.solution_tree = query_assembler.solution_tree
        gibson_assembly, gibson_fragments = gibson_assembler.design(solution=0)
        gibson_solution_service(gibson_obj, gibson_assembler, gibson_assembly, gibson_fragments)
        bundle.gibson.add(gibson_obj)
    if bundle_data['goldengate']:
        goldengate_obj = GoldenGateAssembly(
            title=bundle_data['title'],
            multipart=bundle_data['multipart'],
            addgene=bundle_data['addgene'],
            igem=bundle_data['igem'],
            dnasu=bundle_data['dnasu'],
            min_blast=bundle_data['min_blast'],
            max_blast=bundle_data['max_blast'],
            min_synth=bundle_data['min_synth'],
            max_synth=bundle_data['max_synth'],
            mv_conc=bundle_data['mv_conc'],
            dv_conc=bundle_data['dv_conc'],
            dntp_conc=bundle_data['dntp_conc'],
            dna_conc=bundle_data['dna_conc'],
            tm=bundle_data['tm'],
            multi_query=bundle_data['multi_query'],
            overhangs=bundle_data['overhangs'],
            scarless=bundle_data['scarless']
        )
        goldengate_obj.backbone_file.save(bundle_data['backbone_file'].name, File(open(backbone_file_path, 'rb')))
        goldengate_obj.insert_file.save(bundle_data['insert_file'].name, File(open(insert_file_path, 'rb')))
        goldengate_obj.save()

        goldengate_assembler = GoldenGateAssembler(
            bundle_data['mv_conc'], 
            bundle_data['dv_conc'], 
            bundle_data['dna_conc'],
            bundle_data['dntp_conc'], 
            bundle_data['tm'], 
            backbone_file_path, 
            insert_file_path, 
            db_list(bundle_data['addgene'], bundle_data['igem'], bundle_data['dnasu']), 
            min_frag=bundle_data['min_blast'], 
            max_frag=bundle_data['max_blast'], 
            min_synth=bundle_data['min_synth'], 
            max_synth=bundle_data['max_synth'],
            ovhngs=int(bundle_data['overhangs']),
            multi_query=bundle_data['multi_query'],
            scarless=bundle_data['scarless']
        )
        goldengate_assembler.solution_tree = query_assembler.solution_tree
        goldengate_assembly, goldengate_fragments = goldengate_assembler.design(solution=0)
        goldengate_solution_service(goldengate_obj, goldengate_assembler, goldengate_assembly, goldengate_fragments)
        bundle.goldengate.add(goldengate_obj)
    if bundle_data['biobricks']:
        biobricks_obj = BioBricksAssembly(
            title=bundle_data['title'],
            multipart=bundle_data['multipart'],
            addgene=bundle_data['addgene'],
            igem=bundle_data['igem'],
            dnasu=bundle_data['dnasu'],
            min_blast=bundle_data['min_blast'],
            max_blast=bundle_data['max_blast'],
            min_synth=bundle_data['min_synth'],
            max_synth=bundle_data['max_synth'],
            mv_conc=bundle_data['mv_conc'],
            dv_conc=bundle_data['dv_conc'],
            dntp_conc=bundle_data['dntp_conc'],
            dna_conc=bundle_data['dna_conc'],
            tm=bundle_data['tm'],
            multi_query=bundle_data['multi_query'],
        )
        biobricks_obj.backbone_file.save(bundle_data['backbone_file'].name, File(open(backbone_file_path, 'rb')))
        biobricks_obj.insert_file.save(bundle_data['insert_file'].name, File(open(insert_file_path, 'rb')))        
        biobricks_obj.save()

        biobricks_assembler = BioBrickAssembler(
            bundle_data['mv_conc'], 
            bundle_data['dv_conc'], 
            bundle_data['dna_conc'],
            bundle_data['dntp_conc'], 
            bundle_data['tm'], 
            backbone_file_path, 
            insert_file_path, 
            db_list(bundle_data['addgene'], bundle_data['igem'], bundle_data['dnasu']), 
            min_frag=bundle_data['min_blast'], 
            max_frag=bundle_data['max_blast'], 
            min_synth=bundle_data['min_synth'], 
            max_synth=bundle_data['max_synth'],
            multi_query=bundle_data['multi_query']
        )
        biobricks_assembler.solution_tree = query_assembler.solution_tree
        biobricks_assembly, biobricks_fragments = biobricks_assembler.design(solution=0)
        biobricks_solution_service(biobricks_obj, biobricks_assembler, biobricks_assembly, biobricks_fragments)
        bundle.biobricks.add(biobricks_obj)
    if bundle_data['pcr']:
        pcr_obj = PCRAssembly(
            title=bundle_data['title'],
            multipart=bundle_data['multipart'],
            addgene=bundle_data['addgene'],
            igem=bundle_data['igem'],
            dnasu=bundle_data['dnasu'],
            min_blast=bundle_data['min_blast'],
            max_blast=bundle_data['max_blast'],
            min_synth=bundle_data['min_synth'],
            max_synth=bundle_data['max_synth'],
            mv_conc=bundle_data['mv_conc'],
            dv_conc=bundle_data['dv_conc'],
            dntp_conc=bundle_data['dntp_conc'],
            dna_conc=bundle_data['dna_conc'],
            tm=bundle_data['tm'],
            multi_query=bundle_data['multi_query'],
            overlap=bundle_data['overlap']
        )
        pcr_obj.backbone_file.save(bundle_data['backbone_file'].name, File(open(backbone_file_path, 'rb')))
        pcr_obj.insert_file.save(bundle_data['insert_file'].name, File(open(insert_file_path, 'rb')))
        pcr_obj.save()

        pcr_assembler = PCRAssembler(
            bundle_data['mv_conc'], 
            bundle_data['dv_conc'], 
            bundle_data['dna_conc'],
            bundle_data['dntp_conc'], 
            bundle_data['tm'], 
            backbone_file_path, 
            insert_file_path, 
            db_list(bundle_data['addgene'], bundle_data['igem'], bundle_data['dnasu']), 
            min_frag=bundle_data['min_blast'], 
            max_frag=bundle_data['max_blast'], 
            min_synth=bundle_data['min_synth'], 
            max_synth=bundle_data['max_synth'],
            overlap=bundle_data['overlap'],
            multi_query=bundle_data['multi_query']
        )
        pcr_assembler.solution_tree = query_assembler.solution_tree
        pcr_assembly, pcr_fragments = pcr_assembler.design(solution=0)
        pcr_solution_service(pcr_obj, pcr_assembler, pcr_assembly, pcr_fragments)
        bundle.pcr.add(pcr_obj)
    if bundle_data['slic']:
        slic_obj = SLICAssembly(
            title=bundle_data['title'],
            multipart=bundle_data['multipart'],
            addgene=bundle_data['addgene'],
            igem=bundle_data['igem'],
            dnasu=bundle_data['dnasu'],
            min_blast=bundle_data['min_blast'],
            max_blast=bundle_data['max_blast'],
            min_synth=bundle_data['min_synth'],
            max_synth=bundle_data['max_synth'],
            mv_conc=bundle_data['mv_conc'],
            dv_conc=bundle_data['dv_conc'],
            dntp_conc=bundle_data['dntp_conc'],
            dna_conc=bundle_data['dna_conc'],
            tm=bundle_data['tm'],
            multi_query=bundle_data['multi_query'],
            overlap=bundle_data['overlap']
        )
        slic_obj.backbone_file.save(bundle_data['backbone_file'].name, File(open(backbone_file_path, 'rb')))
        slic_obj.insert_file.save(bundle_data['insert_file'].name, File(open(insert_file_path, 'rb')))
        slic_obj.save()

        slic_assembler = SLICAssembler(
            bundle_data['mv_conc'], 
            bundle_data['dv_conc'], 
            bundle_data['dna_conc'],
            bundle_data['dntp_conc'], 
            bundle_data['tm'], 
            backbone_file_path, 
            insert_file_path, 
            db_list(bundle_data['addgene'], bundle_data['igem'], bundle_data['dnasu']), 
            min_frag=bundle_data['min_blast'], 
            max_frag=bundle_data['max_blast'], 
            min_synth=bundle_data['min_synth'], 
            max_synth=bundle_data['max_synth'],
            overlap=bundle_data['overlap'],
            multi_query=bundle_data['multi_query']
        )
        slic_assembler.solution_tree = query_assembler.solution_tree
        slic_assembly, slic_fragments = slic_assembler.design(solution=0)
        slic_solution_service(slic_obj, slic_assembler, slic_assembly, slic_fragments)
        bundle.slic.add(slic_obj)

    os.remove(backbone_file_path)
    os.remove(insert_file_path)

    return bundle.pk




