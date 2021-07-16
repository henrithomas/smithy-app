from assemblies.assembler import Assembler
from pydna.design import assembly_fragments
from Bio.Seq import Seq
from pydna.primer import Primer
from pydna.amplify import pcr
from importlib import import_module

# TODO citation
class OverlapExtensionAssembler(Assembler):
    cloning_type = 'Overlap Extension'

    def __init__(self, overlap, lin_enzyme, *args, **kwargs):
        super(OverlapExtensionAssembler, self).__init__(*args, **kwargs)
        self.overlap = overlap
        self.lin_enzyme = getattr(import_module('Bio.Restriction'), lin_enzyme)
    
    def primer_extension(self, fragments_pcr, backbone_pcr):
        # returns assembly, a list of amplicons with extensions
        fragments_pcr.insert(0, backbone_pcr)
        fragments_pcr.append(backbone_pcr)
        assembly = assembly_fragments(fragments_pcr, overlap=self.overlap)
        backbone_pcr = pcr(assembly[-1].forward_primer, assembly[0].reverse_primer, backbone_pcr)
        assembly = assembly[1:-1]
        assembly.insert(0, backbone_pcr)
        return assembly

    def design(self, solution=0):
        # return fragments
        fragments = self.get_solution(solution)
        nodes = self.solution_tree.solution_nodes(solution)
        # digest backbone
        # backbone_digest = self.digest_backbone(self.lin_enzyme)
        # create assembly primer complements for backbone and fragments
        fragments_pcr, backbone_pcr = self.primer_complement(fragments, self.backbone)
        # create primer extensions
        assembly = self.primer_extension(fragments_pcr, backbone_pcr)
        assembly = self.annotations(assembly, nodes)
        # create expected assembly construct/seq
        # run primer thermo analysis
        return assembly, nodes