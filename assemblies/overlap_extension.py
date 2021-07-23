from re import A, template
from assemblies.assembler import Assembler
from pydna.design import assembly_fragments
from Bio.Seq import Seq
from pydna.primer import Primer
from pydna.amplify import pcr, Anneal
from importlib import import_module
from primers.analysis import assembly_thermo


# TODO citation
class OverlapExtensionAssembler(Assembler):
    cloning_type = 'Overlap Extension'

    def __init__(self, overlap, lin_enzyme, *args, **kwargs):
        super(OverlapExtensionAssembler, self).__init__(*args, **kwargs)
        self.overlap = overlap
        self.lin_enzyme = getattr(import_module('Bio.Restriction'), lin_enzyme)
    
    def add_extensions(self, ext1, ext2, amp):
        ext_fwd = Seq(ext1) + amp.forward_primer
        ext_rvs = Seq(ext2) + amp.reverse_primer

        ext_fwd_p = Primer(ext_fwd, position=amp.forward_primer.position, footprint=amp.forward_primer._fp,id=amp.forward_primer.id)
        ext_rev_p = Primer(ext_rvs, position=amp.reverse_primer.position, footprint=amp.reverse_primer._fp,id=amp.reverse_primer.id)
        
        annealing = Anneal([ext_fwd_p, ext_rev_p], amp.template)  
        amp_ext = annealing.products[0]
        return amp_ext

    def primer_extension(self, fragments_pcr, backbone_pcr):
        assembly = []
        fragments_pcr.append(backbone_pcr)
        end_idx = len(fragments_pcr) - 1
        overlap = int(self.overlap / 2)

        for i, amp in enumerate(fragments_pcr):
            if i == 0:
                overlap_fwd = fragments_pcr[end_idx].seq.watson[-overlap:]
            else:
                overlap_fwd = fragments_pcr[i - 1].seq.watson[-overlap:]

            if i == end_idx:
                overlap_rvs = fragments_pcr[0].seq.crick[-overlap:]
            else:
                overlap_rvs = fragments_pcr[i + 1].seq.crick[-overlap:]

            amp_extended = self.add_extensions(overlap_fwd, overlap_rvs, amp)
            
            assembly.append(amp_extended)
        
        return assembly
        

    def design(self, solution=0):
        # return fragments
        fragments = self.get_solution(solution)
        nodes = self.solution_tree.solution_nodes(solution)
        
        # create assembly primer complements for backbone and fragments
        fragments_pcr, backbone_pcr = self.primer_complement(fragments, self.backbone)
        
        # create primer extensions
        assembly = self.primer_extension(fragments_pcr, backbone_pcr)
        assembly = self.annotations(assembly, nodes)
        # TODO create expected assembly construct/seq
        # run primer thermo analysis
        assembly = assembly_thermo(assembly, self.mv_conc, self.dv_conc, self.dna_conc, self.tm_custom)
  
        return assembly, nodes

    def design_no_query(self):
        # take in direct dseqrecord list of the parts and backbone from form/model input
        # run self.primer_complement(fragments, backbone)
        # assembly = =self.primer_extension(fragments_pcr, backbone_pcr)
        # make sure names are correct for the assembly fragments & backbone compared to the input frags & backbone
        # run thermo analysis

        pass