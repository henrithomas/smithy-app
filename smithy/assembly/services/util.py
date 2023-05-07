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
from django.conf import settings

img_path = settings.MEDIA_IMAGES_ROOT
csv_path = settings.MEDIA_CSV_ROOT
smithy_path = '/home/dkoch/smithy-app/smithy/' 

#with open('/home/hthoma/files/config.json') as config_file:
#    config = json.load(config_file)

def part_indexes(idx, length):
    left = length - 1 if idx == 0 else idx - 1
    right = 0 if idx == length - 1 else idx + 1 
    return left, right


def pcr_clusters(nt_lengths):
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
    clusters = defaultdict(list)

    for l in nt_lengths:
        for c, lower, upper in cluster_bounds:
            if l > lower and l <= upper:
                clusters[c].append(l)

    return clusters

def lengths_and_plasmids(assembly):
    primer_lengths = []
    part_lengths_pcr = []
    part_lengths_synth = []
    plasmid_count = 0

    for part in assembly[:-1]:
        primer_lengths.extend(
            [
                len(part.forward_primer.seq),
                len(part.reverse_primer.seq)
            ]
        )
        part_lengths_pcr.append(part.template.seq.length)

        if part.annotations['db'] == 'NONE':
            part_lengths_synth.append(part.template.seq.length)
        else:
            plasmid_count += 1

    primer_lengths.extend(
        [
            len(assembly[-1].forward_primer.seq),
            len(assembly[-1].reverse_primer.seq)
        ]
    )
    part_lengths_pcr.append(assembly[-1].template.seq.length)

    return primer_lengths, part_lengths_synth, part_lengths_pcr, plasmid_count

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
    temp_file = os.path.join(csv_path, file_name)
    # temp_file = f'{smithy_path}media/csv/{file_name}'

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
    temp_file = os.path.join(csv_path, file_name)
    # temp_file = f'{smithy_path}media/csv/{file_name}'

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
    temp_file = os.path.join(csv_path, file_name)
    # temp_file = f'{smithy_path}media/csv/{file_name}'

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

def assembly_costs(nt_costs, nt_lengths, plasmid_count, enz_costs, enz_types):
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

def solution_analysis(assembly, fragments, query_length):
    """
    
    Parameters
    ----------
    

    Returns
    -------
    tuple[match_p, synth_p, part_ave, primer_ave, primer_tm_ave, part_max, part_min, db_parts, synth_parts]
    """
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

def part_map(part_model, part, left, right, name, space):
    """
    

    Parameters
    ----------


    Returns
    -------
    """
    part_plot_name = f'{name.replace(" ", "-")}-map.png'
    temp_plot = os.path.join(img_path, part_plot_name)
    # temp_plot = f'{smithy_path}media/images/{part_plot_name}'

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
    with plt.ioff():
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
    temp_plot = os.path.join(img_path, plot_name)
    # temp_plot = f'{smithy_path}media/images/{plot_name}'   
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
    with plt.ioff():
        ax, _ = record.plot(figure_width=10)
        ax.figure.savefig(temp_plot)
        plt.close(ax.get_figure())

    solution_model.plasmid_map.save(plot_name, File(open(temp_plot, 'rb')))
    solution_model.save()
    os.remove(temp_plot)
