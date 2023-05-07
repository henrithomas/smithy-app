from .utility import *
from ..models import GibsonPart, GibsonSolution, GibsonPrimer
from assemblies.gibson import GibsonAssembler
import json
from math import log10

def calculate_time(fragments, part_lengths_pcr):
    pcr = pcr_time(part_lengths_pcr)
    time = assembly_times(pcr, len(fragments))
    return time

def assembly_times(pcr, parts_count):
    # NEB protocol
    assembly = 0.25 if parts_count <= 3 else 1.0
    times = {}
    times_types = ['pcr', 'assembly']
    time_vals = [pcr, assembly]

    times.update({'total': round(sum(time_vals), 2)})
    times.update({'types': times_types})
    times.update({'times': time_vals})

    return times

def mastermix_data(obj):
    """
    Formats enzyme mastermix data for the assembly.

    Parameters
    ----------
    mastermix_cost : float

    mastermix_n_reacts : int

    Returns
    -------
    Lists of enzyme orders, costs, and types data.
    """
    orders = ['Master mix enzyme']
    costs = [obj.mastermix_cost / obj.mastermix_n_reacts]
    types = ['Master mix'] 

    return orders, costs, types

def enzyme_data(obj):
    """
    Formats enzyme data for the assembly.
    
    Parameters
    ----------
    obj : GibsonAssembly

    Returns
    -------
    Lists of enzyme orders, costs, and types data.
    """
    orders = []
    costs = []
    types = ['Exonuclease', 'Ligase', 'Polymerase']
    
    if obj.exonuclease_cost > 0.0:
        orders.append('T5 exonuclease')
    if obj.ligase_cost > 0.0:
        orders.append('Taq ligase')
    if obj.polymerase_cost > 0.0:
        orders.append('Phusion polymerase')
    costs = [
        obj.exonuclease_cost / obj.exonuclease_n_reacts, 
        obj.ligase_cost / obj.ligase_n_reacts, 
        obj.polymerase_cost / obj.polymerase_n_reacts
    ]

    return orders, costs, types

def pcr_polymerase_data(enzyme_orders, enzyme_costs, enzyme_types, pcr_poly_cost):
    enzyme_orders.append('PCR polymerase')
    enzyme_costs.append(pcr_poly_cost)
    enzyme_types.append('PCR polymerase')

def assembly_enzymes(obj, fragments):
    enzyme_orders = []
    enzyme_costs = []
    enzyme_types = []
    
    if obj.mastermix_cost > 0.0:
        enzyme_orders, enzyme_costs, enzyme_types = mastermix_data(obj)
    else:
        enzyme_orders, enzyme_costs, enzyme_types = enzyme_data(obj)

    if obj.pcr_polymerase_cost > 0.0:
        pcr_poly_cost = (len(fragments) + 1) * (obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts)
        pcr_polymerase_data(enzyme_orders, enzyme_costs, enzyme_types, pcr_poly_cost)
        
    return enzyme_orders, enzyme_costs, enzyme_types

def assembly_risk(pcr_ps, assembly_ps):
    risk = {
        'total': 0.35,
        'types': ['PCR', 'Assembly'],
        'risks': [
            log10((1 - pcr_ps) / pcr_ps), 
            log10((1 - assembly_ps) / assembly_ps)
        ]
    }
    
    return risk

def save_parts_and_primers(assembly, solution):
    for i, part in enumerate(assembly):
        db_part = GibsonPart(
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
            subject_end = part.annotations['subject_end'] 
        )
        db_part.save()

        name = f'{part.name}-{i}'
        left_index, right_index = part_indexes(i, len(assembly)) 

        part_map(db_part, part, assembly[left_index], assembly[right_index], name, 0)

        forward_primer = GibsonPrimer(
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

        reverse_primer = GibsonPrimer(
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

def gibson_solution_service(obj, assembler, assembly, fragments):
    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length
    analysis = solution_analysis(assembly, fragments, assembler.query_record.seq.length)
    primer_lengths, part_lengths, part_lengths_pcr, plasmid_count = lengths_and_plasmids(assembly)
    enzyme_orders, enzyme_costs, enzyme_types = assembly_enzymes(obj, fragments)
    nt_costs = [obj.primer_cost, obj.part_cost, obj.gene_cost]
    nt_lengths = part_lengths + primer_lengths

    cost = assembly_costs(nt_costs, nt_lengths, plasmid_count, enzyme_costs, enzyme_types)
    time = calculate_time(fragments, part_lengths_pcr)
    risk = assembly_risk(obj.pcr_ps, obj.assembly_ps)

    solution = GibsonSolution(
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
        time_summary=json.dumps(time),
        cost_summary=json.dumps(cost),
        risk_summary=json.dumps(risk)
    )
    solution.save()

    plasmid_map(solution, assembly, obj.title, 0, total_len)
    parts_csv(solution, assembly)
    primers_csv(solution, assembly)
    order_csv(solution, assembly, enzyme_orders)

    save_parts_and_primers(assembly, solution)

