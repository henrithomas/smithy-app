from util import *
from gibson import make_gibson_solution
from goldengate import make_goldengate_solution
from pcr import make_pcr_solution
from slic import make_slic_solution
from biobricks import make_biobricks_solution
from ..models import (
    AssemblyBundle,
    GibsonAssembly,
    GoldenGateAssembly,
    PCRAssembly,
    SLICAssembly,
    BioBricksAssembly
)
from assemblies.gibson import GibsonAssembler,SLICAssembler, PCRAssembler
from assemblies.goldengate import GoldenGateAssembler
from assemblies.traditional_re import BioBrickAssembler
from assemblies.assembler import Assembler
from django.core.files import File
from django.conf import settings
import os

def save_gibson_assembly(bundle_data, backbone_file_path, insert_file_path):
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
            pcr_polymerase_cost=bundle_data['pcr_polymerase_cost'],
            pcr_polymerase_n_reacts=bundle_data['pcr_polymerase_n_reacts'],
            assembly_ps=bundle_data['gib_assembly_ps'],
            mastermix_cost=bundle_data['gib_mmix_cost'],
            mastermix_n_reacts=bundle_data['gib_mmix_n_reacts'],
            pcr_ps=bundle_data['pcr_ps']
        )
    gibson_obj.backbone_file.save(bundle_data['backbone_file'].name, File(open(backbone_file_path, 'rb')))
    gibson_obj.insert_file.save(bundle_data['insert_file'].name, File(open(insert_file_path, 'rb')))
    gibson_obj.save()
    return gibson_obj

def save_goldengate_assembly(bundle_data, backbone_file_path, insert_file_path):
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
            pcr_polymerase_cost=bundle_data['pcr_polymerase_cost'],
            pcr_polymerase_n_reacts=bundle_data['pcr_polymerase_n_reacts'],
            assembly_ps=bundle_data['gg_assembly_ps'],
            mastermix_cost=bundle_data['gg_mmix_cost'],
            mastermix_n_reacts=bundle_data['gg_mmix_n_reacts'],
            pcr_ps=bundle_data['pcr_ps']
        )
    goldengate_obj.backbone_file.save(bundle_data['backbone_file'].name, File(open(backbone_file_path, 'rb')))
    goldengate_obj.insert_file.save(bundle_data['insert_file'].name, File(open(insert_file_path, 'rb')))
    goldengate_obj.save()
    return goldengate_obj

def gibson_design(bundle_data, backbone_file_path, insert_file_path, solution_tree, gibson_obj):
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
    gibson_assembler.solution_tree = solution_tree
    gibson_assembly, gibson_fragments = gibson_assembler.design(solution=0)
    make_gibson_solution(gibson_obj, gibson_assembler, gibson_assembly, gibson_fragments)

def goldengate_design(bundle_data, backbone_file_path, insert_file_path, solution_tree, goldengate_obj):
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
    goldengate_assembler.solution_tree = solution_tree
    goldengate_assembly, goldengate_fragments = goldengate_assembler.design(solution=0)
    make_goldengate_solution(goldengate_obj, goldengate_assembler, goldengate_assembly, goldengate_fragments)

def pcr_design(bundle_data, backbone_file_path, insert_file_path, solution_tree, pcr_obj):
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
    pcr_assembler.solution_tree = solution_tree
    pcr_assembly, pcr_fragments = pcr_assembler.design(solution=0)
    make_pcr_solution(pcr_obj, pcr_assembler, pcr_assembly, pcr_fragments)

def save_pcr_assembly(bundle_data, backbone_file_path, insert_file_path):
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
            pcr_polymerase_n_reacts=bundle_data['pcr_polymerase_n_reacts'],
            mastermix_cost=bundle_data['pcr_mmix_cost'],
            mastermix_n_reacts=bundle_data['pcr_mmix_n_reacts'],
            pcr_ps=bundle_data['pcr_ps']
        )
    pcr_obj.backbone_file.save(bundle_data['backbone_file'].name, File(open(backbone_file_path, 'rb')))
    pcr_obj.insert_file.save(bundle_data['insert_file'].name, File(open(insert_file_path, 'rb')))
    pcr_obj.save()
    return pcr_obj

def slic_design(bundle_data, backbone_file_path, insert_file_path, solution_tree, slic_obj):
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
    slic_assembler.solution_tree = solution_tree
    slic_assembly, slic_fragments = slic_assembler.design(solution=0)
    make_slic_solution(slic_obj, slic_assembler, slic_assembly, slic_fragments)

def save_slic_assembly(bundle_data, backbone_file_path, insert_file_path):
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
            chewback_ps=bundle_data['slic_chewback_ps'],
            pcr_ps=bundle_data['pcr_ps']
        )
    slic_obj.backbone_file.save(bundle_data['backbone_file'].name, File(open(backbone_file_path, 'rb')))
    slic_obj.insert_file.save(bundle_data['insert_file'].name, File(open(insert_file_path, 'rb')))
    slic_obj.save()
    return slic_obj

def biobricks_design(bundle_data, backbone_file_path, insert_file_path, solution_tree, biobricks_obj):
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
    biobricks_assembler.solution_tree = solution_tree
    biobricks_assembly, biobricks_fragments = biobricks_assembler.design(solution=0)
    make_biobricks_solution(biobricks_obj, biobricks_assembler, biobricks_assembly, biobricks_fragments)

def save_biobricks_assembly(bundle_data, backbone_file_path, insert_file_path):
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
            pcr_polymerase_cost=bundle_data['pcr_polymerase_cost'],
            pcr_polymerase_n_reacts=bundle_data['pcr_polymerase_n_reacts'],
            digestion_ps=bundle_data['bb_digestion_ps'],
            ligation_ps=bundle_data['bb_ligation_ps'],
            pcr_ps=bundle_data['pcr_ps']
        )
    biobricks_obj.backbone_file.save(bundle_data['backbone_file'].name, File(open(backbone_file_path, 'rb')))
    biobricks_obj.insert_file.save(bundle_data['insert_file'].name, File(open(insert_file_path, 'rb')))        
    biobricks_obj.save()
    return biobricks_obj

def bundle_create_service(bundle_data):
    # save backbone and insert files here
    # smithy-app/smithy/media
    path = os.path.join(settings.MEDIA_ROOT, 'fasta')
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
        gibson_obj = save_gibson_assembly(bundle_data, backbone_file_path, insert_file_path)

        gibson_design(bundle_data, backbone_file_path, insert_file_path, query_assembler.solution_tree, gibson_obj)
        bundle.gibson.add(gibson_obj)
  
    if bundle_data['goldengate']:
        print('bundle - goldengate')
        goldengate_obj = save_goldengate_assembly(bundle_data, backbone_file_path, insert_file_path)

        goldengate_design(bundle_data, backbone_file_path, insert_file_path, query_assembler.solution_tree, goldengate_obj)
        bundle.goldengate.add(goldengate_obj)
  
    if bundle_data['pcr']:
        print('bundle - pcr')
        pcr_obj = save_pcr_assembly(bundle_data, backbone_file_path, insert_file_path)

        pcr_design(bundle_data, backbone_file_path, insert_file_path, query_assembler.solution_tree, pcr_obj)
        bundle.pcr.add(pcr_obj)

    if bundle_data['slic']:
        print('bundle - slic')
        slic_obj = save_slic_assembly(bundle_data, backbone_file_path, insert_file_path)

        slic_design(bundle_data, backbone_file_path, insert_file_path, query_assembler.solution_tree, slic_obj)
        bundle.slic.add(slic_obj)

    if bundle_data['biobricks']:
        print('bundle - biobricks')
        biobricks_obj = save_biobricks_assembly(bundle_data, backbone_file_path, insert_file_path)

        biobricks_design(bundle_data, backbone_file_path, insert_file_path, query_assembler.solution_tree, biobricks_obj)
        bundle.biobricks.add(biobricks_obj)
 
    os.remove(backbone_file_path)
    os.remove(insert_file_path)
    print('bundle - complete')
    return bundle.pk
