from .util import *
from ..models import GoldenGatePart, GoldenGatePrimer, GoldenGateSolution
from assemblies.goldengate import GoldenGateAssembler
import json
from math import log10

def calculate_time(len_fragments, part_lengths):
    pcr = pcr_time(part_lengths)
    time = assembly_times(pcr, len_fragments)
    return time

def assembly_times(pcr, parts_count):
    # NEB protocol
    times = {}
    time_types = ['pcr', 'assembly']
    time_vals = [pcr]

    if parts_count == 1:
        time_vals.append(0.1)
    elif parts_count >= 2 and parts_count < 5:
        time_vals.append(1.1)
    elif parts_count >= 5 and parts_count < 11:
        time_vals.append(1.1)
    elif parts_count >= 11:
        time_vals.append(5.01)

    times.update({'total': round(sum(time_vals), 2)})
    times.update({'types': time_types})
    times.update({'times': time_vals})

    return times

def enzymes_data(obj, len_fragments):
    orders = []
    costs = []
    types = []

    if obj.mastermix_cost > 0.0:
        orders.append('Master mix enzyme')
        costs = [ obj.mastermix_cost / obj.mastermix_n_reacts ]
        types = ['Master mix']
    else:
        if obj.re_cost > 0.0:
            orders.append('Type2S')

        if obj.ligase_cost > 0.0:
            orders.append('Ligase')

        costs.extend([
            obj.re_cost / obj.re_n_reacts, 
            obj.ligase_cost / obj.ligase_n_reacts
        ])

        types.extend(['Type2s RE', 'Ligase'])
        
    if obj.pcr_polymerase_cost > 0.0:
        pcr_poly_cost = (len_fragments + 1) * (obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts)
        orders.append('PCR polymerase')
        costs.append(pcr_poly_cost)
        types.append('PCR polymerase')

    return orders, costs, types

def assembly_risk(pcr_ps, assembly_ps):
    return {
        'total': 0.35,
        'types': ['PCR', 'Assembly'],
        'risks': [
            log10((1 - pcr_ps) / pcr_ps), 
            log10((1 - assembly_ps) / assembly_ps)
        ]
    }

def save_parts_and_primers(assembly, space, solution):
    for i, part in enumerate(assembly):
        db_part = GoldenGatePart(
            name=part.name,
            database=part.annotations['db'],
            length=part.template.seq.length, 
            length_extended=part.seq.length,
            seq=part.template.seq,
            seq_extended=part.seq,
            position=i,
            solution=solution,
            query_start = part.annotations['query_start'],
            query_end = part.annotations['query_end'],
            subject_start = part.annotations['subject_start'],
            subject_end = part.annotations['subject_end'],
            cuts=part.annotations['cuts'],
            cut_locations=json.dumps(part.annotations['cut_locations'])           
        )
        db_part.save()

        name = f'{part.name}-{i}'
        left_index, right_index = part_indexes(i, len(assembly))

        part_map(db_part, part, assembly[left_index], assembly[right_index], name, space)

        forward_primer = GoldenGatePrimer(
            name= f'{db_part.name} forward primer',
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
            part=db_part
        )
        forward_primer.save()

        reverse_primer = GoldenGatePrimer(
            name= f'{db_part.name} reverse primer ',
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
            part=db_part 
        )
        reverse_primer.save()

def make_goldengate_solution(obj, assembler, assembly, fragments):
    space = 0 if obj.scarless else 4
    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length + space
    analysis = solution_analysis(assembly, fragments, assembler.query_record.seq.length)
    primer_lengths, part_lengths, part_lengths_pcr, plasmid_count = lengths_and_plasmids(assembly)
    enzyme_orders, enzyme_costs, enzyme_types = enzymes_data(obj, len(fragments))
    nt_costs = [obj.primer_cost, obj.part_cost, obj.gene_cost]
    nt_lengths = part_lengths + primer_lengths
    enzyme_orders = []

    cost = assembly_costs(nt_costs, nt_lengths, plasmid_count, enzyme_costs, enzyme_types)
    time = calculate_time(len(fragments), part_lengths_pcr)
    risk = assembly_risk(obj.pcr_ps, obj.assembly_ps)

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
        time_summary=json.dumps(time),
        cost_summary=json.dumps(cost),
        risk_summary=json.dumps(risk)
    )
    goldengate_solution.save()

    plasmid_map(goldengate_solution, assembly, obj.title, space, total_len)
    parts_csv(goldengate_solution, assembly)
    primers_csv(goldengate_solution, assembly)
    order_csv(goldengate_solution, assembly, enzyme_orders)

    save_parts_and_primers(assembly, space, goldengate_solution)

def run_goldengate(obj):
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

    make_goldengate_solution(obj, assembler, assembly, fragments)