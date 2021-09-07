from django.db.models.deletion import Collector
from .models import (
    GibsonAssembly, 
    GibsonPart,
    GibsonPrimer,
    GoldenGateAssembly,
    GoldenGatePart,
    GoldenGatePrimer,
    GibsonSolution,
    GoldenGateSolution,
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
from dna_features_viewer import (
    GraphicFeature, 
    GraphicRecord, 
    CircularGraphicRecord
)
from django.core.files import File
import os


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

    solution_model.plasmid_map.save(plot_name, File(open(temp_plot, 'rb')))
    solution_model.save()
    os.remove(temp_plot)

def gibson_create_service(gibson_obj):
    """
    


    Parameters
    ----------



    Returns
    -------
    """
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

def goldengate_create_service(goldengate_obj):
    """
    


    Parameters
    ----------



    Returns
    -------
    """
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
        multi_query=goldengate_obj.multi_query,
        scarless=goldengate_obj.scarless
    )
    
    space = 0 if goldengate_obj.scarless else 4

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

    plasmid_map(goldengate_solution, gg_assembly, goldengate_obj.title, space, total_len)

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

def biobricks_create_service(biobricks_obj):
    """
    


    Parameters
    ----------



    Returns
    -------
    """
    biobricks_assembler = BioBrickAssembler(
        biobricks_obj.mv_conc, 
        biobricks_obj.dv_conc, 
        biobricks_obj.dna_conc,
        biobricks_obj.dntp_conc, 
        biobricks_obj.tm, 
        biobricks_obj.backbone_file.path, 
        biobricks_obj.insert_file.path, 
        db_list(biobricks_obj.addgene, biobricks_obj.igem, biobricks_obj.dnasu), 
        min_frag=biobricks_obj.min_blast, 
        max_frag=biobricks_obj.max_blast, 
        min_synth=biobricks_obj.min_synth, 
        max_synth=biobricks_obj.max_synth,
        overlap=biobricks_obj.overlap,
        multi_query=biobricks_obj.multi_query  
    )

    if biobricks_obj.multi_query:
        results, error = biobricks_assembler.run_multi_query()
        biobricks_assembler.multi_query_solution_building(results)
    else:
        results, error = biobricks_assembler.query()
        biobricks_assembler.solution_building(results)
    biobricks_assembly, biobricks_fragments = biobricks_assembler.design(solution=0)

    total_len = biobricks_assembler.backbone.seq.length + biobricks_assembler.query_record.seq.length

    biobricks_solution = BioBricksSolution(
        name=f'Solution - {biobricks_obj.title}',
        backbone=biobricks_assembler.backbone.seq,
        query=biobricks_assembler.query_record.seq,
        solution='',
        parts_count=len(biobricks_fragments),
        primers_count=len(biobricks_fragments) * 2,
        match=0.0,
        assembly=biobricks_obj
    )
    biobricks_solution.save()

    plasmid_map(biobricks_solution, biobricks_assembly, biobricks_obj.title, 0, total_len)


    for i, part in enumerate(biobricks_assembly):
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
            subject_end = part.annotations['subject_end'] 
        )
        biobricks_part_entry.save()

        if i == 0:
            left_index = len(biobricks_assembly) - 1
        else:
            left_index = i - 1

        if i == len(biobricks_assembly) - 1:
            right_index = 0
        else:
            right_index = i + 1

        part_map(
            biobricks_part_entry, 
            part, 
            biobricks_assembly[left_index], 
            biobricks_assembly[right_index], 
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

def pcr_create_service(pcr_obj):
    """
    


    Parameters
    ----------



    Returns
    -------
    """
    pcr_assembler = PCRAssembler(
        pcr_obj.mv_conc,
        pcr_obj.dv_conc,
        pcr_obj.dna_conc,
        pcr_obj.dntp_conc,
        pcr_obj.tm,
        pcr_obj.backbone_file.path,
        pcr_obj.insert_file.path,
        db_list(pcr_obj.addgene, pcr_obj.igem, pcr_obj.dnasu),
        min_frag=pcr_obj.min_blast,
        max_frag=pcr_obj.max_blast,
        min_synth=pcr_obj.min_synth,
        max_synth=pcr_obj.max_synth,
        overlap=pcr_obj.overlap,
        multi_query=pcr_obj.multi_query
    )

    if pcr_obj.multi_query:
        results, error = pcr_assembler.run_multi_query()
        pcr_assembler.multi_query_solution_building(results)
    else:
        results, error = pcr_assembler.query()
        pcr_assembler.solution_building(results)
    pcr_assembly, pcr_fragments = pcr_assembler.design(solution=0)

    total_len = pcr_assembler.backbone.seq.length + pcr_assembler.query_record.seq.length

    pcr_solution = PCRSolution(
        name=f'Solution - {pcr_obj.title}',
        backbone=pcr_assembler.backbone.seq,
        query=pcr_assembler.query_record.seq,
        solution='',
        parts_count=len(pcr_fragments),
        primers_count=len(pcr_fragments) * 2, 
        match=0.0,
        assembly=pcr_obj
    )
    pcr_solution.save()

    plasmid_map(pcr_solution, pcr_assembly, pcr_obj.title, 0, total_len)


    for i, part in enumerate(pcr_assembly):
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
            left_index = len(pcr_assembly) - 1
        else:
            left_index = i - 1

        if i == len(pcr_assembly) - 1:
            right_index = 0
        else:
            right_index = i + 1

        part_map(
            pcr_part_entry, 
            part, 
            pcr_assembly[left_index], 
            pcr_assembly[right_index], 
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

def slic_create_service(slic_obj):
    """
    


    Parameters
    ----------



    Returns
    -------
    """
    slic_assembler = SLICAssembler(
        slic_obj.mv_conc,
        slic_obj.dv_conc,
        slic_obj.dna_conc,
        slic_obj.dntp_conc,
        slic_obj.tm,
        slic_obj.backbone_file.path,
        slic_obj.insert_file.path,
        db_list(slic_obj.addgene, slic_obj.igem, slic_obj.dnasu),
        min_frag=slic_obj.min_blast,
        max_frag=slic_obj.max_blast,
        min_synth=slic_obj.min_synth,
        max_synth=slic_obj.max_synth,
        overlap=slic_obj.overlap,
        multi_query=slic_obj.multi_query
    )

    if slic_obj.multi_query:
        results, error = slic_assembler.run_multi_query()
        slic_assembler.multi_query_solution_building(results)
    else:
        results, error = slic_assembler.query()
        slic_assembler.solution_building(results)
    slic_assembly, slic_fragments = slic_assembler.design(solution=0)

    total_len = slic_assembler.backbone.seq.length + slic_assembler.query_record.seq.length

    slic_solution = SLICSolution(
        name=f'Solution - {slic_obj.title}',
        backbone=slic_assembler.backbone.seq,
        query=slic_assembler.query_record.seq,
        solution='',
        parts_count=len(slic_fragments),
        primers_count=len(slic_fragments) * 2, 
        match=0.0,
        assembly=slic_obj
    )
    slic_solution.save()

    plasmid_map(slic_solution, slic_assembly, slic_obj.title, 0, total_len)


    for i, part in enumerate(slic_assembly):
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
            left_index = len(slic_assembly) - 1
        else:
            left_index = i - 1

        if i == len(slic_assembly) - 1:
            right_index = 0
        else:
            right_index = i + 1

        part_map(
            slic_part_entry, 
            part, 
            slic_assembly[left_index], 
            slic_assembly[right_index], 
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
