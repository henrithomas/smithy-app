from assemblies.assembler import Assembler
from assemblies.traditional_re import TraditionalREAssembler
from Bio.Seq import Seq
from pydna.primer import Primer
from pydna.amplify import pcr
from importlib import import_module
import json

# TODO citation
class GoldenGateAssembler(TraditionalREAssembler):
    cloning_type = 'Golden Gate'
    ligase = 'T4-DNA'
    restriction_enzyme = getattr(import_module('Bio.Restriction'), 'BsaI')
    cutsite = restriction_enzyme.site
    # potapov et al 2018
    overhang_files = [
        '/home/hthoma/projects/smithy-app/assemblies/overhangs/goldengate/oh1.json',
        '/home/hthoma/projects/smithy-app/assemblies/overhangs/goldengate/oh2.json',
        '/home/hthoma/projects/smithy-app/assemblies/overhangs/goldengate/oh3.json',
        '/home/hthoma/projects/smithy-app/assemblies/overhangs/goldengate/oh4.json'
    ]

    def __init__(self, *args, ovhngs=0, re='BsaI', ligase='T4-DNA', **kwargs):
        super(GoldenGateAssembler, self).__init__(re, *args, ligase=ligase, **kwargs)
        with open(self.overhang_files[ovhngs]) as f:
            self.overhangs = json.load(f)

    def primer_extension(self, fragments_pcr, backbone_pcr):
        # returns assembly, a list of amplicons with extensions
        assembly = []
        fragments_pcr.append(backbone_pcr)
        set_size = len(fragments_pcr)
        for i, amp in enumerate(fragments_pcr):
            overhang_fwd = self.overhangs['seq'][i]
            if i == set_size - 1:
                overhang_rev = str(Seq(self.overhangs['seq'][0]).reverse_complement())
            else:
                overhang_rev = str(Seq(self.overhangs['seq'][i+1]).reverse_complement())

            re_site_fwd = self.re1.site + 'n' + overhang_fwd
            re_site_rev = self.re1.site + 'n' + overhang_rev
            amp_extended = self.add_cutsites(re_site_fwd, re_site_rev, amp)
            
            assembly.append(amp_extended)
            
        return assembly

    def design(self, solution=0):
        # return fragments
        fragments = self.get_solution(solution)
        nodes = self.solution_tree.solution_nodes(solution)
        # digest backbone
        # backbone_digest = self.digest_backbone(self.restriction_enzyme)
        # create assembly primer complements for backbone and fragments
        fragments_pcr, backbone_pcr = self.primer_complement(fragments, self.backbone)
        # create primer extensions
        assembly = self.primer_extension(fragments_pcr, backbone_pcr)
        assembly = self.annotations(assembly, nodes)
        # create expected assembly construct/seq
        # run primer thermo analysis
        return assembly
