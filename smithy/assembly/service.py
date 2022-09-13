from django.db.models.deletion import Collector
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
from datetime import datetime, time
import csv
import matplotlib.pyplot as plt
from statistics import mean
from math import log10
from collections import defaultdict
from math import inf

#with open('/home/hthoma/files/config.json') as config_file:
#    config = json.load(config_file)

smithy_path = '/home/dkoch/smithy-app/smithy/' 

def part_indexes(idx, length):
    left = length - 1 if idx == 0 else idx - 1
    right = 0 if idx == length - 1 else idx + 1 
    return left, right

cluster_bounds = (
    (1, 0, 1000),
    (2, 1000, 2000),
    (3, 2000, 3000),
    (4, 3000, 4000),
    (5, 4000, 5000),
    (6, 5000, 6000),
    (7, 6000, 7000),
    (8, 7000, 8000),
    (9, 8000, 9000),
    (10, 9000, inf),
)

def pcr_clusters(nt_lengths):
    clusters = defaultdict(list)

    for l in nt_lengths:
        for c, lower, upper in cluster_bounds:
            if l > lower and l <= upper:
                clusters[c].append(l)

    return clusters

def lengths_and_plasmids(assembly):
    primer_lengths = []
    part_lengths = []
    plasmid_count = 0

    for part in assembly[:-1]:
        primer_lengths.extend(
            [
                len(part.forward_primer.seq),
                len(part.reverse_primer.seq)
            ]
        )

        if part.annotations['db'] == 'NONE':
            part_lengths.append(part.template.seq.length)
        else:
            plasmid_count += 1

    primer_lengths.extend(
        [
            len(assembly[-1].forward_primer.seq),
            len(assembly[-1].reverse_primer.seq)
        ]
    )
    return primer_lengths, part_lengths, plasmid_count

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
    temp_file = f'{smithy_path}media/csv/{file_name}'

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
    temp_file = f'{smithy_path}media/csv/{file_name}'

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

def order_csv(solution_model, parts, enzymes):
    file_name = f'{solution_model.name.replace(" ", "-")}-{solution_model.pk}-order.csv'
    temp_file = f'{smithy_path}media/csv/{file_name}'

    fields = ['id', 'type', 'sequence']

    part_list = []
    primer_list = []
    enzyme_list = []

    for part in parts: 
        part_list.append(
            {
                'id': part.name,
                'type': 'part',
                'sequence': part.template.seq.watson
            }
        )
        primer_list.extend([
            {
                'id': f'{part.name}-fwd',
                'type': 'primer',
                'sequence': part.forward_primer.seq
            },
            {
                'id': f'{part.name}-rvs',
                'type': 'primer',
                'sequence': part.reverse_primer.seq
            }
        ])

    for enzyme in enzymes:
        enzyme_list.append(
            {
                'id': enzyme,
                'type': 'enzyme',
                'sequence': 'none'
            }
        )

    with open(temp_file, 'w', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fields, restval='NONE')
        writer.writeheader()
        writer.writerows(part_list)
        writer.writerows(primer_list)
        writer.writerows(enzyme_list)

    solution_model.order_file.save(file_name, File(open(temp_file, newline='')))
    solution_model.save()
    os.remove(temp_file)

def costs(nt_costs, nt_lengths, plasmid_count, enz_costs, enz_types):
    nt_types = ['oligos', 'blocks', 'genes', 'plasmids']
    nt_totals = [0.0, 0.0, 0.0, 0.0]
    costs = {}
    total_cost = 0

    # TODO: change categories to : oligo, block, gene only for non-db parts
    # TODO: add plasmid flat costs for db parts
    for length in nt_lengths: 
        if length <= 100:
            # cost of primer
            nt_totals[0] += length * nt_costs[0]
        elif length > 100 and length <= 1000:
            # cost of part
            nt_totals[1] += length * nt_costs[1]
        elif length > 1000:
            # cost of gene
            nt_totals[2] += length * nt_costs[2]
        
    # cost of plasmid
    nt_totals[3] = plasmid_count * 75.0
    
    total_cost = sum(nt_totals) + sum(enz_costs)

    cost_types = nt_types + enz_types
    cost_indv = nt_totals + enz_costs

    costs.update({'total': round(total_cost, 2)})
    costs.update({'types': cost_types})
    costs.update({'costs': cost_indv})

    return costs

def pcr_time(nt_lengths):
    pcr_total = 0.0

    clusters = pcr_clusters(nt_lengths)

    for cluster, lengths in clusters.items():
        max_length = max(lengths)
        step = 0.5 * cluster
        denature = step + 0.4
        anneal = step + 0.4
        polymerization = step + 0.4
        time = 1.0 + 30 * (denature + anneal + polymerization) + 10.0
        pcr_total += time

    return round(pcr_total / 60, 2)

def goldengate_times(pcr, insert_count):
    # NEB protocol
    times = {}
    time_types = ['pcr', 'digestion-ligation']
    time_vals = [pcr]

    if insert_count == 1:
        time_vals.append(0.1)
    elif insert_count >= 2 and insert_count < 5:
        time_vals.append(1.1)
    elif insert_count >= 5 and insert_count < 11:
        time_vals.append(1.1)
    elif insert_count >= 11:
        time_vals.append(5.01)

    times.update({'total': round(sum(time_vals), 2)})
    times.update({'types': time_types})
    times.update({'times': time_vals})

    return times

def gibson_times(pcr):
    # NEB protocol
    times = {}
    times_types = ['pcr', 'chewback, ligation, repair']
    time_vals = [pcr, 1.0]

    times.update({'total': round(sum(time_vals), 2)})
    times.update({'types': times_types})
    times.update({'times': time_vals})

    return times

def slic_times(pcr, overlap):
    # slic methods article
    times = {}
    times_types = ['pcr', 'chewback', 'ligation']
    chewback = 0
    cooling = 0.1
    ligation = 0.5 + cooling

    if overlap < 40: 
        chewback = 0.8
    else:
        chewback = 1.35

    time_vals = [pcr, chewback, ligation]

    times.update({'total': round(sum(time_vals), 2)})
    times.update({'types': times_types})
    times.update({'times': time_vals})

    return times

def biobricks_times(pcr, insert_count):
    # iGEM and Ginko Bioworks protocol 
    times = {}
    times_types = ['pcr', 'digestion', 'ligation']
    time_vals = [pcr, 0.8, 0.5]
    
    times.update({'total': round(sum(time_vals), 2)})
    times.update({'types': times_types})
    times.update({'times': time_vals})

    return times

def pcr_soe_times(nt_lengths):
    # TODO: needs more detail 
    times = {}
    time_types = ['pcr']
    time_vals = [pcr_time(nt_lengths)]

    times.update({'total': round(sum(time_vals), 2)})
    times.update({'types': time_types})
    times.update({'times': time_vals})

    return times

def goldengate_risk(insert_count):
    pcr_pots = insert_count
    assembly_pots = 1
    pass

def gibson_risk():
    pots = 1
    pass

def slic_risk(insert_count):
    pots = insert_count + 1
    pass

def biobricks_risk(insert_count):
    pots = (insert_count - 1) * 4
    pass

def pcr_risk(insert_count):
    pots = insert_count
    pass

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
    synth_p = round(100.00 - match_p, 2)
    part_ave = int(mean(part_lengths))
    primer_ave = int(mean(primer_lengths))
    primer_tm_ave = round(mean(primer_tms), 2)
    part_max = max(part_lengths)
    part_min = min(part_lengths)
    db_parts = sum(not f.synth for f in fragments)
    synth_parts = sum(f.synth for f in fragments) 
    return match_p, synth_p, part_ave, primer_ave, primer_tm_ave, part_max, part_min, db_parts, synth_parts    

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
    part_plot_name = f'{name.replace(" ", "-")}-map.png'
    temp_plot = f'{smithy_path}media/images/{part_plot_name}'

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
    plot_name = f'{assembly_name.replace(" ", "-")}_map.png'
    temp_plot = f'{smithy_path}media/images/{plot_name}'   
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
        assembler.multi_query_solution_building(
            results,
            assembly_type='gibson',
            part_costs=[
                obj.primer_cost,
                obj.part_cost,
                obj.gene_cost
            ],
            pcr_polymerase_cost=(obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts),
            polymerase_cost=(obj.polymerase_cost / obj.polymerase_n_reacts),
            exonuclease_cost=(obj.exonuclease_cost / obj.exonuclease_n_reacts),
            ligase_cost=(obj.ligase_cost / obj.ligase_n_reacts),
            parts_pref=obj.parts_pref,
            cost_pref=obj.cost_pref,
            pcr_ps=obj.pcr_ps,
            gibson_ps=obj.assembly_ps
        )
    else:
        results, error = assembler.query()
        assembler.solution_building(
            results,
            assembly_type='gibson',
            part_costs=[
                obj.primer_cost,
                obj.part_cost,
                obj.gene_cost
            ],
            pcr_polymerase_cost=(obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts),
            polymerase_cost=(obj.polymerase_cost / obj.polymerase_n_reacts),
            exonuclease_cost=(obj.exonuclease_cost / obj.exonuclease_n_reacts),
            ligase_cost=(obj.ligase_cost / obj.ligase_n_reacts),
            parts_pref=obj.parts_pref,
            cost_pref=obj.cost_pref,
            pcr_ps=obj.pcr_ps,
            gibson_ps=obj.assembly_ps
        )
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
        assembler.multi_query_solution_building(
            results,
            assembly_type='goldengate',
            part_costs=[
                obj.primer_cost,
                obj.part_cost,
                obj.gene_cost
            ],
            pcr_polymerase_cost=(obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts),
            ligase_cost=(obj.ligase_cost / obj.ligase_n_reacts),
            restenz_cost=(obj.re_cost / obj.re_n_reacts),
            parts_pref=obj.parts_pref,
            cost_pref=obj.cost_pref,
            pcr_ps=obj.pcr_ps,
            goldengate_ps=obj.assembly_ps
        )
    else:
        results, error = assembler.query()
        assembler.solution_building(
            results,
            assembly_type='goldengate',
            part_costs=[
                obj.primer_cost,
                obj.part_cost,
                obj.gene_cost
            ],
            pcr_polymerase_cost=(obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts),
            ligase_cost=(obj.ligase_cost / obj.ligase_n_reacts),
            restenz_cost=(obj.re_cost / obj.re_n_reacts),
            parts_pref=obj.parts_pref,
            cost_pref=obj.cost_pref,
            pcr_ps=obj.pcr_ps,
            goldengate_ps=obj.assembly_ps
        )
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

    digest_cost = 2 * ((obj.EcoRI_cost + obj.XbaI_cost + obj.SpeI_cost + obj.PstI_cost) / 4)
    digest_n = 2 * ((obj.EcoRI_n_reacts + obj.XbaI_n_reacts + obj.SpeI_n_reacts + obj.PstI_n_reacts) / 4)

    if obj.multi_query:
        results, error = assembler.run_multi_query()
        assembler.multi_query_solution_building(
            results,
            assembly_type='biobricks',
            part_costs=[
                obj.primer_cost,
                obj.part_cost,
                obj.gene_cost
            ],
            pcr_polymerase_cost=(obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts),
            ligase_cost=(obj.ligase_cost / obj.ligase_n_reacts),
            biobricks_digest_cost=(digest_cost / digest_n),
            parts_pref=obj.parts_pref,
            cost_pref=obj.cost_pref,
            pcr_ps=obj.pcr_ps,
            ligation_ps=obj.ligation_ps,
            biobricks_digest_ps=obj.digestion_ps
        )
    else:
        results, error = assembler.query()
        assembler.solution_building(
            results,
            assembly_type='biobricks',
            part_costs=[
                obj.primer_cost,
                obj.part_cost,
                obj.gene_cost
            ],
            pcr_polymerase_cost=(obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts),
            ligase_cost=(obj.ligase_cost / obj.ligase_n_reacts),
            biobricks_digest_cost=(digest_cost / digest_n),
            parts_pref=obj.parts_pref,
            cost_pref=obj.cost_pref,
            pcr_ps=obj.pcr_ps,
            ligation_ps=obj.ligation_ps,
            biobricks_digest_ps=obj.digestion_ps
        )
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
        assembler.multi_query_solution_building(
            results,
            assembly_type='pcrsoe',
            part_costs=[
                obj.primer_cost,
                obj.part_cost,
                obj.gene_cost
            ],
            pcr_polymerase_cost=(obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts),
            parts_pref=obj.parts_pref,
            cost_pref=obj.cost_pref,
            pcr_ps=obj.pcr_ps
        )
    else:
        results, error = assembler.query()
        assembler.solution_building(
            results,
            assembly_type='pcrsoe',
            part_costs=[
                obj.primer_cost,
                obj.part_cost,
                obj.gene_cost
            ],
            pcr_polymerase_cost=(obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts),
            parts_pref=obj.parts_pref,
            cost_pref=obj.cost_pref,
            pcr_ps=obj.pcr_ps
        )
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
        assembler.multi_query_solution_building(
            results,
            assembly_type='slic',
            part_costs=[
                obj.primer_cost,
                obj.part_cost,
                obj.gene_cost
            ],
            pcr_polymerase_cost=(obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts),
            exonuclease_cost=(obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts),
            parts_pref=obj.parts_pref,
            cost_pref=obj.cost_pref,
            pcr_ps=obj.pcr_ps,
            slic_exo_ps=obj.chewback_ps
        )
    else:
        results, error = assembler.query()
        assembler.solution_building(
            results,
            assembly_type='slic',
            part_costs=[
                obj.primer_cost,
                obj.part_cost,
                obj.gene_cost
            ],
            pcr_polymerase_cost=(obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts),
            exonuclease_cost=(obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts),
            parts_pref=obj.parts_pref,
            cost_pref=obj.cost_pref,
            pcr_ps=obj.pcr_ps,
            slic_exo_ps=obj.chewback_ps
        )
    assembly, fragments = assembler.design(solution=0)

    slic_solution_service(obj, assembler, assembly, fragments)    

def gibson_solution_service(obj, assembler, assembly, fragments):
    # Log based odds of success: risk = log((1 - P_s) / P_s)
    # pcr: P_s = 0.8
    # chewback, ligation, repair: P_s = 0.7

    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length
    # match_p, synth_p, part_ave, primer_ave, primer_tm_ave, part_max, part_min, db_parts, synth_parts
    analysis = solution_analysis(assembly, fragments, assembler.query_record.seq.length)
    primer_lengths, part_lengths, plasmid_count = lengths_and_plasmids(assembly)
    enzyme_orders = []
    
    if obj.exonuclease_cost > 0.0:
        enzyme_orders.append('T5 exonuclease')
    if obj.ligase_cost > 0.0:
        enzyme_orders.append('Taq ligase')
    if obj.polymerase_cost > 0.0:
        enzyme_orders.append('Phusion polymerase')

    pcr = pcr_time(part_lengths + primer_lengths)
    gibson_time = gibson_times(pcr)
    gibson_cost = costs(
        [obj.primer_cost, obj.part_cost, obj.gene_cost],
        part_lengths + primer_lengths,
        plasmid_count,
        [
            obj.exonuclease_cost / obj.exonuclease_n_reacts, 
            obj.ligase_cost / obj.ligase_n_reacts, 
            obj.polymerase_cost / obj.polymerase_n_reacts
        ],
        ['exonuclease', 'ligase', 'polymerase']
    )
    gibson_risk = {
        'total': 0.35,
        'types': ['pcr', 'chewback, ligation, repair'],
        'risks': [
            log10((0.2) / 0.8), 
            log10((0.3) / 0.7)
        ]
    }

    gibson_solution = GibsonSolution(
        name=f'{obj.title} Solution',
        backbone=assembler.backbone.seq,
        query=assembler.query_record.seq,
        solution='',
        parts_count=len(fragments),
        primers_count=len(fragments) * 2,
        match=analysis[0],
        synth_amount=analysis[1],
        re_enzymes=False,
        part_length_average=analysis[2],
        primer_length_average=analysis[3],
        tm_average=analysis[4],
        longest_part=analysis[5],
        shortest_part=analysis[6],
        db_parts=analysis[7],
        synth_parts=analysis[8],
        solution_length=assembler.query_record.seq.length,
        assembly=obj,
        time_summary=json.dumps(gibson_time),
        cost_summary=json.dumps(gibson_cost),
        risk_summary=json.dumps(gibson_risk)
    )
    gibson_solution.save()

    plasmid_map(gibson_solution, assembly, obj.title, 0, total_len)
    parts_csv(gibson_solution, assembly)
    primers_csv(gibson_solution, assembly)
    order_csv(gibson_solution, assembly, enzyme_orders)

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

        left_index, right_index = part_indexes(i, len(assembly))

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
    # Log based odds of success: risk = log((1 - P_s) / P_s)
    # pcr: P_s = 0.8
    # digestion, ligation: P_s = 0.9

    space = 0 if obj.scarless else 4
    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length + space
    # match_p, synth_p, part_ave, primer_ave, primer_tm_ave, part_max, part_min, db_parts, synth_parts
    analysis = solution_analysis(assembly, fragments, assembler.query_record.seq.length)
    primer_lengths, part_lengths, plasmid_count = lengths_and_plasmids(assembly)
    enzyme_orders = []
    
    if obj.re_cost > 0.0:
        enzyme_orders.append('Type2S')
    if obj.ligase_cost > 0.0:
        enzyme_orders.append('Ligase')

    pcr = pcr_time(part_lengths + primer_lengths)
    goldengate_time = goldengate_times(pcr, len(fragments))
    goldengate_cost = costs(
        [obj.primer_cost, obj.part_cost, obj.gene_cost],
        part_lengths + primer_lengths,
        plasmid_count,
        [
            obj.re_cost / obj.re_n_reacts, 
            obj.ligase_cost / obj.ligase_n_reacts
        ],
        ['type2s RE', 'ligase']
    )
    goldengate_risk = {
        'total': 0.35,
        'types': ['pcr', 'digestion, ligation'],
        'risks': [
            log10((0.2) / 0.8), 
            log10((0.1) / 0.9)
        ]
    }

    # TODO update to have a match % and BLAST solution sequence
    # TODO add a foreach solution in for the assembly
    goldengate_solution = GoldenGateSolution(
        name=f'{obj.title} Solution',
        backbone=assembler.backbone.seq,
        query=assembler.query_record.seq,
        solution='',
        parts_count=len(fragments),
        primers_count=len(fragments) * 2,
        match=analysis[0],
        synth_amount=analysis[1],
        re_enzymes=True,
        part_length_average=analysis[2],
        primer_length_average=analysis[3],
        tm_average=analysis[4],
        longest_part=analysis[5],
        shortest_part=analysis[6],
        db_parts=analysis[7],
        synth_parts=analysis[8],
        solution_length=assembler.query_record.seq.length,
        assembly=obj,
        time_summary=json.dumps(goldengate_time),
        cost_summary=json.dumps(goldengate_cost),
        risk_summary=json.dumps(goldengate_risk)
    )
    goldengate_solution.save()

    plasmid_map(goldengate_solution, assembly, obj.title, space, total_len)
    parts_csv(goldengate_solution, assembly)
    primers_csv(goldengate_solution, assembly)
    order_csv(goldengate_solution, assembly, enzyme_orders)

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

        left_index, right_index = part_indexes(i, len(assembly))

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
    # Log based odds of success: risk = log((1 - P_s) / P_s)
    # pcr: P_s = 0.8
    # digestion: P_s = 0.9 
    # ligation: P_s = 0.8
    cost_coefficient = 2 * len(fragments) - 1
    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length
    # match_p, synth_p, part_ave, primer_ave, primer_tm_ave, part_max, part_min, db_parts, synth_parts
    analysis = solution_analysis(assembly, fragments, assembler.query_record.seq.length)
    primer_lengths, part_lengths, plasmid_count = lengths_and_plasmids(assembly)
    enzyme_orders = []
    
    if obj.EcoRI_cost > 0.0:
        enzyme_orders.append('EcoRI')
    if obj.XbaI_cost > 0.0:
        enzyme_orders.append('XbaI')
    if obj.SpeI_cost > 0.0:
        enzyme_orders.append('SpeI')
    if obj.PstI_cost > 0.0:
        enzyme_orders.append('PstI')

    pcr = pcr_time(part_lengths + primer_lengths)
    biobricks_time = biobricks_times(pcr, len(fragments))
    biobricks_cost = costs(
        [obj.primer_cost, obj.part_cost, obj.gene_cost],
        part_lengths + primer_lengths,
        plasmid_count,
        [
            cost_coefficient * (obj.EcoRI_cost / obj.EcoRI_n_reacts), 
            cost_coefficient * (obj.XbaI_cost / obj.XbaI_n_reacts), 
            cost_coefficient * (obj.SpeI_cost / obj.SpeI_n_reacts), 
            cost_coefficient * (obj.PstI_cost / obj.PstI_n_reacts)
        ],
        ['EcoRI', 'XbaI', 'SpeI', 'PstI']
    )
    biobricks_risk = {
        'total': 0.35,
        'types': ['pcr', 'digestion', 'ligation'],
        'risks': [
            log10((0.2) / 0.8), 
            log10((0.1) / 0.9), 
            log10((0.2) / 0.8)
        ]
    }

    biobricks_solution = BioBricksSolution(
        name=f'{obj.title} Solution',
        backbone=assembler.backbone.seq,
        query=assembler.query_record.seq,
        solution='',
        parts_count=len(fragments),
        primers_count=len(fragments) * 2,
        match=analysis[0],
        synth_amount=analysis[1],
        re_enzymes=True,
        part_length_average=analysis[2],
        primer_length_average=analysis[3],
        tm_average=analysis[4],
        longest_part=analysis[5],
        shortest_part=analysis[6],
        db_parts=analysis[7],
        synth_parts=analysis[8],
        solution_length=assembler.query_record.seq.length,
        assembly=obj,
        time_summary=json.dumps(biobricks_time),
        cost_summary=json.dumps(biobricks_cost),
        risk_summary=json.dumps(biobricks_risk)
    )
    biobricks_solution.save()

    plasmid_map(biobricks_solution, assembly, obj.title, 0, total_len)
    parts_csv(biobricks_solution, assembly)
    primers_csv(biobricks_solution, assembly)
    order_csv(biobricks_solution, assembly, enzyme_orders)

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

        left_index, right_index = part_indexes(i, len(assembly))

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
    # Log based odds of success: risk = log((1 - P_s) / P_s)
    # pcr: P_s = 0.8

    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length
    # match_p, synth_p, part_ave, primer_ave, primer_tm_ave, part_max, part_min, db_parts, synth_parts
    analysis = solution_analysis(assembly, fragments, assembler.query_record.seq.length)
    primer_lengths, part_lengths, plasmid_count = lengths_and_plasmids(assembly)
    enzyme_orders = []
    
    if obj.pcr_polymerase_cost > 0.0:
        enzyme_orders.append('Phusion polymerase')

    pcr = pcr_time(part_lengths + primer_lengths)
    pcr_soe_time = pcr_soe_times(part_lengths + primer_lengths)
    pcr_cost = costs(
        [obj.primer_cost, obj.part_cost, obj.gene_cost],
        part_lengths + primer_lengths,
        plasmid_count,
        [obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts],
        ['polymerase']
    )
    pcr_risk = {
        'total': 0.35,
        'types': ['amplification-recombination'],
        'risks': [
            log10((0.2) / 0.8)
        ]
    }

    pcr_solution = PCRSolution(
        name=f'{obj.title} Solution',
        backbone=assembler.backbone.seq,
        query=assembler.query_record.seq,
        solution='',
        parts_count=len(fragments),
        primers_count=len(fragments) * 2, 
        match=analysis[0],
        synth_amount=analysis[1],
        re_enzymes=False,
        part_length_average=analysis[2],
        primer_length_average=analysis[3],
        tm_average=analysis[4],
        longest_part=analysis[5],
        shortest_part=analysis[6],
        db_parts=analysis[7],
        synth_parts=analysis[8],
        solution_length=assembler.query_record.seq.length,
        assembly=obj,
        time_summary=json.dumps(pcr_soe_time),
        cost_summary=json.dumps(pcr_cost),
        risk_summary=json.dumps(pcr_risk)
    )
    pcr_solution.save()

    plasmid_map(pcr_solution, assembly, obj.title, 0, total_len)
    parts_csv(pcr_solution, assembly)
    primers_csv(pcr_solution, assembly)
    order_csv(pcr_solution, assembly, enzyme_orders)

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

        left_index, right_index = part_indexes(i, len(assembly))

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
    # Log based odds of success: risk = log((1 - P_s) / P_s)
    # pcr: P_s = 0.8
    # chewback: P_s = 0.8 
    # ligation: P_s = 0.8

    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length
    # match_p, synth_p, part_ave, primer_ave, primer_tm_ave, part_max, part_min, db_parts, synth_parts
    analysis = solution_analysis(assembly, fragments, assembler.query_record.seq.length)
    primer_lengths, part_lengths, plasmid_count = lengths_and_plasmids(assembly)
    enzyme_orders = []
    
    if obj.exonuclease_cost > 0.0:
        enzyme_orders.append('T5 exonuclease')
    if obj.ligase_cost > 0.0:
        enzyme_orders.append('Taq ligase')

    pcr = pcr_time(part_lengths + primer_lengths)
    slic_time = slic_times(pcr, obj.overlap)
    slic_cost = costs(
        [obj.primer_cost, obj.part_cost, obj.gene_cost, obj.plasmid_cost],
        part_lengths + primer_lengths,
        plasmid_count,
        [
            obj.exonuclease_cost / obj.exonuclease_n_reacts, 
            obj.ligase_cost / obj.ligase_n_reacts
        ],
        ['exonuclease', 'ligase']
    )
    slic_risk = {
        'total': 0.35,
        'types': ['pcr', 'chewback', 'ligation'],
        'risks': [
            log10((0.2) / 0.8), 
            log10((0.2) / 0.8), 
            log10((0.2) / 0.8)
        ]
    }

    slic_solution = SLICSolution(
        name=f'{obj.title} Solution',
        backbone=assembler.backbone.seq,
        query=assembler.query_record.seq,
        solution='',
        parts_count=len(fragments),
        primers_count=len(fragments) * 2, 
        match=analysis[0],
        synth_amount=analysis[1],
        re_enzymes=False,
        part_length_average=analysis[2],
        primer_length_average=analysis[3],
        tm_average=analysis[4],
        longest_part=analysis[5],
        shortest_part=analysis[6],
        db_parts=analysis[7],
        synth_parts=analysis[8],
        solution_length=assembler.query_record.seq.length,
        assembly=obj,
        time_summary=json.dumps(slic_time),
        cost_summary=json.dumps(slic_cost),
        risk_summary=json.dumps(slic_risk)
    )
    slic_solution.save()

    plasmid_map(slic_solution, assembly, obj.title, 0, total_len)
    parts_csv(slic_solution, assembly)
    primers_csv(slic_solution, assembly)
    order_csv(slic_solution, assembly, enzyme_orders)

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

        left_index, right_index = part_indexes(i, len(assembly))

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
    path = f'{smithy_path}media/fasta'
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
        query_assembler.multi_query_solution_building(
            results,
            assembly_type='none',
            part_costs=[
                bundle_data['primer_cost'],
                bundle_data['part_cost'],
                bundle_data['gene_cost']
            ],
            pcr_polymerase_cost=(bundle_data['pcr_polymerase_cost'] / bundle_data['pcr_polymerase_n_reacts']),
            parts_pref=bundle_data['parts_pref'],
            cost_pref=bundle_data['cost_pref'],
            pcr_ps=bundle_data['pcr_ps']
        )
    else:
        results, error = query_assembler.query()
        query_assembler.solution_building(
            results,
            assembly_type='none',
            part_costs=[
                bundle_data['primer_cost'],
                bundle_data['part_cost'],
                bundle_data['gene_cost']
            ],
            pcr_polymerase_cost=(bundle_data['pcr_polymerase_cost'] / bundle_data['pcr_polymerase_n_reacts']),
            parts_pref=bundle_data['parts_pref'],
            cost_pref=bundle_data['cost_pref'],
            pcr_ps=bundle_data['pcr_ps']
        )

    bundle = AssemblyBundle(
        title=bundle_data['title'],
        description=bundle_data['description']
    )
    bundle.save()

    if bundle_data['gibson']:
        print('bundle - gibson')
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
            overlap=bundle_data['gib_overlap'],
            primer_cost=bundle_data['primer_cost'],
            part_cost=bundle_data['part_cost'],
            gene_cost=bundle_data['gene_cost'],
            exonuclease_cost=bundle_data['gib_exonuclease_cost'],
            ligase_cost=bundle_data['gib_ligase_cost'],
            polymerase_cost=bundle_data['gib_polymerase_cost'],
            exonuclease_n_reacts=bundle_data['gib_exonuclease_n_reacts'],
            ligase_n_reacts=bundle_data['gib_ligase_n_reacts'],
            polymerase_n_reacts=bundle_data['gib_polymerase_n_reacts'],
            assembly_ps=bundle_data['gib_assembly_ps']
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
            overlap=bundle_data['gib_overlap'],
            multi_query=bundle_data['multi_query']
        )
        gibson_assembler.solution_tree = query_assembler.solution_tree
        gibson_assembly, gibson_fragments = gibson_assembler.design(solution=0)
        gibson_solution_service(gibson_obj, gibson_assembler, gibson_assembly, gibson_fragments)
        bundle.gibson.add(gibson_obj)
    if bundle_data['goldengate']:
        print('bundle - goldengate')
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
            scarless=bundle_data['scarless'],
            primer_cost=bundle_data['primer_cost'],
            part_cost=bundle_data['part_cost'],
            gene_cost=bundle_data['gene_cost'],
            re_cost=bundle_data['gg_re_cost'],
            ligase_cost=bundle_data['gg_ligase_cost'],
            re_n_reacts=bundle_data['gg_re_n_reacts'],
            ligase_n_reacts=bundle_data['gg_ligase_n_reacts'],
            assembly_ps=bundle_data['gg_assembly_ps']
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
        print('bundle - biobricks')
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
            primer_cost=bundle_data['primer_cost'],
            part_cost=bundle_data['part_cost'],
            gene_cost=bundle_data['gene_cost'],
            EcoRI_cost=bundle_data['bb_EcoRI_cost'],
            XbaI_cost=bundle_data['bb_XbaI_cost'],
            SpeI_cost=bundle_data['bb_SpeI_cost'],
            PstI_cost=bundle_data['bb_PstI_cost'],
            EcoRI_n_reacts=bundle_data['bb_EcoRI_n_reacts'],
            XbaI_n_reacts=bundle_data['bb_XbaI_n_reacts'],
            SpeI_n_reacts=bundle_data['bb_SpeI_n_reacts'],
            PstI_n_reacts=bundle_data['bb_PstI_n_reacts'],
            ligase_cost=bundle_data['bb_ligase_cost'],
            ligase_n_reacts=bundle_data['bb_ligase_n_reacts'],
            digestion_ps=bundle_data['bb_digestion_ps'],
            ligation_ps=bundle_data['bb_ligation_ps']
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
        print('bundle - pcr')
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
            overlap=bundle_data['pcr_overlap'],
            primer_cost=bundle_data['primer_cost'],
            part_cost=bundle_data['part_cost'],
            gene_cost=bundle_data['gene_cost'],
            pcr_polymerase_cost=bundle_data['pcr_polymerase_cost'],
            pcr_polymerase_n_reacts=bundle_data['pcr_polymerase_n_reacts']
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
            overlap=bundle_data['pcr_overlap'],
            multi_query=bundle_data['multi_query']
        )
        pcr_assembler.solution_tree = query_assembler.solution_tree
        pcr_assembly, pcr_fragments = pcr_assembler.design(solution=0)
        pcr_solution_service(pcr_obj, pcr_assembler, pcr_assembly, pcr_fragments)
        bundle.pcr.add(pcr_obj)
    if bundle_data['slic']:
        print('bundle - slic')
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
            overlap=bundle_data['slic_overlap'],
            primer_cost=bundle_data['primer_cost'],
            part_cost=bundle_data['part_cost'],
            gene_cost=bundle_data['gene_cost'],
            exonuclease_cost=bundle_data['slic_exonuclease_cost'],
            exonuclease_n_reacts=bundle_data['slic_exonuclease_n_reacts'],
            chewback_ps=bundle_data['slic_chewback_ps']
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
            overlap=bundle_data['slic_overlap'],
            multi_query=bundle_data['multi_query']
        )
        slic_assembler.solution_tree = query_assembler.solution_tree
        slic_assembly, slic_fragments = slic_assembler.design(solution=0)
        slic_solution_service(slic_obj, slic_assembler, slic_assembly, slic_fragments)
        bundle.slic.add(slic_obj)

    os.remove(backbone_file_path)
    os.remove(insert_file_path)
    print('bundle - complete')
    return bundle.pk




