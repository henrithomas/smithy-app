from .utility import *
from ..models import GoldenGatePart, GoldenGatePrimer, GoldenGateSolution
from assemblies.goldengate import GoldenGateAssembler
import json
from math import log10

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

def goldengate_solution_service(obj, assembler, assembly, fragments):
    # Log based odds of success: risk = log((1 - P_s) / P_s)
    # pcr: P_s = 0.8
    # digestion, ligation: P_s = 0.9

    space = 0 if obj.scarless else 4
    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length + space
    # match_p, synth_p, part_ave, primer_ave, primer_tm_ave, part_max, part_min, db_parts, synth_parts
    analysis = solution_analysis(assembly, fragments, assembler.query_record.seq.length)
    primer_lengths, part_lengths, part_lengths_pcr, plasmid_count = lengths_and_plasmids(assembly)
    enzyme_orders = []

    if obj.mastermix_cost > 0.0:
        enzyme_orders.append('Master mix enzyme')
        gg_enz_costs = [ obj.mastermix_cost / obj.mastermix_n_reacts ]
        gg_enz_types = ['Master mix']
    else:
        if obj.re_cost > 0.0:
            enzyme_orders.append('Type2S')
        if obj.ligase_cost > 0.0:
            enzyme_orders.append('Ligase')
        gg_enz_costs = [
            obj.re_cost / obj.re_n_reacts, 
            obj.ligase_cost / obj.ligase_n_reacts
        ]
        gg_enz_types = ['Type2s RE', 'Ligase']
        
    if obj.pcr_polymerase_cost > 0.0:
        enzyme_orders.append('PCR polymerase')
        gg_enz_costs.append((len(fragments) + 1) * (obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts))
        gg_enz_types.append('PCR polymerase')

    pcr = pcr_time(part_lengths_pcr)
    goldengate_time = goldengate_times(pcr, len(fragments))
    goldengate_cost = assembly_costs(
        [obj.primer_cost, obj.part_cost, obj.gene_cost],
        part_lengths + primer_lengths,
        plasmid_count,
        gg_enz_costs,
        gg_enz_types
    )
    goldengate_risk = {
        'total': 0.35,
        'types': ['PCR', 'Assembly'],
        'risks': [
            log10((1 - obj.pcr_ps) / obj.pcr_ps), 
            log10((1 - obj.assembly_ps) / obj.assembly_ps)
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
