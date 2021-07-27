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
    """
    A class to implement overlap extension-based cloning methods. It is mostly used as a parent class for other 
    more specific classes like GibsonAssembler.


    Attributes
    ----------
    cloning_type : str
        A simple description of the cloning method

    overlap : int
        Amount of overlap for extension primer designs for an assembly

    lin_enzyme : Bio.Restriction enzyme
        An enzyme instance from Bio.Restriction used to linearize a backbone Dseqrecord

    *args
        See documentation for Assembler parent class

    **kwargs
        See documentation for Assembler parent class


    Returns
    -------
        An instance of OverlapExtensionAssembler
    

    Methods
    -------
    def add_extensions(self, ext1, ext2, amp):
        Adds primer extensions for a given assembly amplicon with unique forward and reverse primer extension sequences.

    def primer_extension(self, fragments_pcr, backbone_pcr):
        Determines the primer extensions for a given assembly amplicon in a parts set.
    
    def design(self, solution=0):
        Runs a full design procedure on the selected solution. Fetches the solution from the solution_tree, adds primer 
        complements, designs and add primer extensions, logs part annotations, and performs thermodynamic analysis.

    """

    cloning_type = 'Overlap Extension'

    def __init__(self, overlap, lin_enzyme, *args, **kwargs):
        """
        Constructs all necessary attributes for the OverlapExtensionAssembler object.


        Parameters
        ----------
        overlap : int
            Amount of overlap for extension primer designs for an assembly

        lin_enzyme : Bio.Restriction enzyme
            An enzyme instance from Bio.Restriction used to linearize a backbone Dseqrecord

        *args
            See documentation for Assembler parent class

        **kwargs
            See documentation for Assembler parent class

        """
        super(OverlapExtensionAssembler, self).__init__(*args, **kwargs)
        self.overlap = overlap
        self.lin_enzyme = getattr(import_module('Bio.Restriction'), lin_enzyme)
    
    def add_extensions(self, ext1, ext2, amp):
        """
        Adds primer extensions for a given assembly amplicon with unique forward and reverse primer extension sequences.


        Parameters
        ----------
        ext1 : str
            Extension sequence for the forward primer of the amplicon

        ext2 : str
            Extension sequence for the reverse primer of the amplicon 

        amp : pydna Amplicon
            A single pydna amplicon object from the assembly
        
        *args
            See documentation for Assembler parent class

        **kwargs
            See documentation for Assembler parent class


        Returns
        -------
        An extended amplicon object for assembly
        """
        ext_fwd = Seq(ext1) + amp.forward_primer
        ext_rvs = Seq(ext2) + amp.reverse_primer

        ext_fwd_p = Primer(
                        ext_fwd, 
                        position=amp.forward_primer.position, 
                        footprint=amp.forward_primer._fp,
                        id=amp.forward_primer.id)
        ext_rev_p = Primer(
                        ext_rvs, 
                        position=amp.reverse_primer.position, 
                        footprint=amp.reverse_primer._fp,
                        id=amp.reverse_primer.id)
        
        annealing = Anneal([ext_fwd_p, ext_rev_p], amp.template)  
        amp_ext = annealing.products[0]
        return amp_ext

    def primer_extension(self, fragments_pcr, backbone_pcr):
        """
        Determines the primer extensions for a given assembly amplicon in a parts set.


        Parameters
        ----------
        fragments_pcr : List of pydna Amplicons  
            A list of pydna Amplicons used for a given assembly solution

        backbone_pcr : A pydna Amplicon
            The backbone Amplicon for a given assembly


        Returns
        -------
        A list of new Amplicon objects, including fragments and backbone, that have gone through primer extension design and 
        extended accordingly
        """
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
        """
        Runs a full design procedure on the selected solution. Fetches the solution from the solution_tree, adds primer 
        complements, designs and add primer extensions, logs part annotations, and performs thermodynamic analysis.
        

        Parameters
        ----------
        solution : int
            The solution index in the solution_tree of the assembler


        Returns
        -------
        A fully design list of assembly parts for assembly with a list of the blast record data for each part (nodes) 
        """

        # return fragments and nodes
        fragments = self.get_solution(solution)
        nodes = self.solution_tree.solution_nodes(solution)
        
        # create assembly primer complements for backbone and fragments
        fragments_pcr, backbone_pcr = self.primer_complement(fragments, self.backbone)
        
        # create primer extensions
        assembly = self.primer_extension(fragments_pcr, backbone_pcr)

        # add simple annotations
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