from .util import *
from ..models import BioBricksPart, BioBricksPrimer, BioBricksSolution
from assemblies.traditional_re import BioBrickAssembler
import json
from math import log10

def calculate_time(fragments, part_lengths_pcr):
    pcr = pcr_time(part_lengths_pcr)
    time = assembly_times(pcr, len(fragments))
    return time 

def assembly_times(pcr, insert_count):
    # iGEM and Ginko Bioworks protocol 
    times = {}
    times_types = ['pcr', 'digestion', 'ligation']
    time_vals = [pcr, 0.8, 0.5]
    
    times.update({'total': round(sum(time_vals), 2)})
    times.update({'types': times_types})
    times.update({'times': time_vals})

    return times

def enzymes_data(obj, len_fragments):
    orders = []
    costs = [
        obj.EcoRI_cost / obj.EcoRI_n_reacts, 
        obj.XbaI_cost / obj.XbaI_n_reacts, 
        obj.SpeI_cost / obj.SpeI_n_reacts, 
        obj.PstI_cost / obj.PstI_n_reacts,
        (len_fragments + 1) * (obj.pcr_polymerase_cost / obj.pcr_polymerase_n_reacts),
        obj.ligase_cost / obj.ligase_n_reacts
    ]
    types = ['EcoRI', 'XbaI', 'SpeI', 'PstI', 'PCR polymerase', 'Ligase']
    
    if obj.EcoRI_cost > 0.0:
        orders.append('EcoRI')
    if obj.XbaI_cost > 0.0:
        orders.append('XbaI')
    if obj.SpeI_cost > 0.0:
        orders.append('SpeI')
    if obj.PstI_cost > 0.0:
        orders.append('PstI')
    if obj.pcr_polymerase_cost > 0.0:
        orders.append('PCR polymerase')

    return orders, costs, types

def assembly_risk(pcr_ps, digestion_ps, ligation_ps):
    biobricks_risk = {
        'total': 0.35,
        'types': ['PCR', 'Digestion', 'Ligation'],
        'risks': [
            log10((1 - pcr_ps) / pcr_ps), 
            log10((1 - digestion_ps) / digestion_ps), 
            log10((1 - ligation_ps) / ligation_ps)
        ]
    }
    
    return biobricks_risk

def save_parts_and_primers(assembly, biobricks_solution):
    for i, part in enumerate(assembly):
        db_part = BioBricksPart(
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
        db_part.save()

        name = f'{part.name}-{i}'
        left_index, right_index = part_indexes(i, len(assembly))

        part_map(db_part, part, assembly[left_index], assembly[right_index], name, 0)

        forward_primer = BioBricksPrimer(
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

        reverse_primer = BioBricksPrimer(
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

def make_biobricks_solution(obj, assembler, assembly, fragments):
    total_len = assembler.backbone.seq.length + assembler.query_record.seq.length
    analysis = solution_analysis(assembly, fragments, assembler.query_record.seq.length)
    primer_lengths, part_lengths, part_lengths_pcr, plasmid_count = lengths_and_plasmids(assembly)
    enzyme_orders, enzyme_costs, enzyme_types = enzymes_data(obj, len(fragments))
    nt_costs = [obj.primer_cost, obj.part_cost, obj.gene_cost]
    nt_lengths = part_lengths + primer_lengths

    time = calculate_time(fragments, part_lengths_pcr)
    biobricks_cost = assembly_costs(nt_costs, nt_lengths, plasmid_count, enzyme_costs, enzyme_types)
    biobricks_risk = assembly_risk(obj.pcr_ps, obj.digestion_ps, obj.ligation_ps)

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
        time_summary=json.dumps(time),
        cost_summary=json.dumps(biobricks_cost),
        risk_summary=json.dumps(biobricks_risk)
    )
    biobricks_solution.save()

    plasmid_map(biobricks_solution, assembly, obj.title, 0, total_len)
    parts_csv(biobricks_solution, assembly)
    primers_csv(biobricks_solution, assembly)
    order_csv(biobricks_solution, assembly, enzyme_orders)

    save_parts_and_primers(assembly, biobricks_solution) 

def run_biobricks(obj):
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

    make_biobricks_solution(obj, assembler, assembly, fragments)   