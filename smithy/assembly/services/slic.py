from .util import *
from ..models import SLICPart, SLICPrimer, SLICSolution
from assemblies.gibson import SLICAssembler
import json
from math import log10

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


def slic_solution_service(obj, assembler, assembly, fragments):
    # Log based odds of success: risk = log((1 - P_s) / P_s)
    # pcr: P_s = 0.8
    # chewback: P_s = 0.8 
    # ligation: P_s = 0.8

    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length
    # match_p, synth_p, part_ave, primer_ave, primer_tm_ave, part_max, part_min, db_parts, synth_parts
    analysis = solution_analysis(assembly, fragments, assembler.query_record.seq.length)
    primer_lengths, part_lengths, part_lengths_pcr, plasmid_count = lengths_and_plasmids(assembly)
    enzyme_orders = []
    
    if obj.pcr_polymerase_cost > 0.0:
        enzyme_orders.append('PCR polymerase')

    slic_enz_costs = [
        (len(fragments) + 1) * (obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts)
    ]
    slic_enz_types = ['PCR polymerase']

    pcr = pcr_time(part_lengths_pcr)
    slic_time = slic_times(pcr, obj.overlap)
    slic_cost = costs(
        [obj.primer_cost, obj.part_cost, obj.gene_cost, obj.plasmid_cost],
        part_lengths + primer_lengths,
        plasmid_count,
        slic_enz_costs,
        slic_enz_types
    )
    slic_risk = {
        'total': 0.35,
        'types': ['PCR', 'Chewback'],
        'risks': [
            log10((1 - obj.pcr_ps) / obj.pcr_ps), 
            log10((1 - obj.chewback_ps) / obj.chewback_ps)
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
