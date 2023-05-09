from assemblies.traditional_re import TraditionalREAssembler
from importlib import import_module
from primers.analysis import assembly_thermo
from pydna.dseqrecord import Dseqrecord

class BioBrickAssembler(TraditionalREAssembler):
    """
    # citation: https://parts.igem.org/Help:BioBrick_Prefix_and_Suffix,
                https://parts.igem.org/Assembly:Standard_assembly
    """
    cloning_type = 'BioBrickAssembly'
    prefix = 'gaattcgcggccgcttctagag'
    prefix_cds = 'gaattcgcggccgcttctag'
    suffix = 'ctgcagcggccgctactagta'

    def __init__(self, *args, re1='EcoRI', re2='XbaI', re3='SpeI', re4='PstI', **kwargs):
        super(BioBrickAssembler, self).__init__(re1, *args, re2=re2, **kwargs)
        self.re3 = getattr(import_module('Bio.Restriction'), re3)
        self.re4 = getattr(import_module('Bio.Restriction'), re4)

    def primer_extension(self, fragments_pcr, backbone_pcr):
        """
        Creates the BioBrick primer extensions for a given assembly amplicon in a parts set.


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
        # returns assembly, a list of amplicons with extensions
        assembly = []
        fragments_pcr.append(backbone_pcr)
        
        for amp in fragments_pcr:
            amp_suffix = self.suffix
            if amp.template[:3].lower() == 'atg':
                amp_prefix = self.prefix_cds
            else:
                amp_prefix = self.prefix
            
            amp_extended = self.add_cutsites(amp_prefix, amp_suffix, amp)

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
        A fully designed list of assembly parts for assembly with a list of the blast record data for each part (nodes) 
        """
        # return fragments and nodes
        if self.multi_query:
            fragments = self.get_multi_query_solution(solution)
            nodes = self.solution_tree.multi_query_solution_nodes(solution)
            self.query_record = Dseqrecord(''.join([record.seq.watson for record in self.query_records]))
        else:
            fragments = self.get_solution(solution)
            nodes = self.solution_tree.solution_nodes(solution)

        # create assembly primer complements for backbone and fragments
        fragments_pcr, backbone_pcr = self.primer_complement(fragments, self.backbone)

        # create primer extensions
        assembly = self.primer_extension(fragments_pcr, backbone_pcr)

        # add simple annotations
        assembly = self.annotations(assembly, nodes)
        
        # run primer thermo analysis
        assembly = assembly_thermo(assembly, self.mv_conc, self.dv_conc, self.dna_conc, self.tm_custom)

        return assembly, nodes