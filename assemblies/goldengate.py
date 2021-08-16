from assemblies.assembler import Assembler
from assemblies.traditional_re import TraditionalREAssembler
from Bio.Seq import Seq
from pydna.primer import Primer
from pydna.amplify import pcr
from importlib import import_module
from primers.analysis import assembly_thermo
import json

# TODO citation
class GoldenGateAssembler(TraditionalREAssembler):
    """
    A class used to implement and execute Golden Gate assemblies. 

    citation: Golden Gate Cloning - Engler et al 2008


    Attributes
    ----------
    cloning_type : str
        A simple description of the cloning method

    overhang_files : list
        A list of file paths to overhang sets of various sizes 

    ovhngs : int
        An index in the overhang_files for the file path to use for the assembly

    re : str
        The name of the restriction enzyme to use for the assembly

    ligase : str
        The name of the ligase to use for the assembly

    *args
        See documentation for TraditionalREAssembler parent class

    **kwargs
        See documentation for TraditionalREAssembler parent class 


    Returns
    -------
    An instance of GoldenGateAssembler


    Methods
    -------
    def primer_extension(self, fragments_pcr, backbone_pcr):
        Determines the primer extensions for a given assembly amplicon in a parts set.

    def design(self, solution=0):
        Runs a full design procedure on the selected solution. Fetches the solution from the solution_tree, adds primer 
        complements, designs and add primer extensions, logs part annotations, and performs thermodynamic analysis.
    """
    
    cloning_type = 'Golden Gate'
    # potapov et al 2018
    overhang_files = [
        '/home/hthoma/projects/smithy-app/assemblies/overhangs/goldengate/oh1.json',
        '/home/hthoma/projects/smithy-app/assemblies/overhangs/goldengate/oh2.json',
        '/home/hthoma/projects/smithy-app/assemblies/overhangs/goldengate/oh3.json',
        '/home/hthoma/projects/smithy-app/assemblies/overhangs/goldengate/oh4.json'
    ]

    def __init__(self, *args, ovhngs=0, re='BsaI', ligase='T4-DNA', **kwargs):
        """
        Constructs all necessary attributes for the GoldenGateAssembler object


        Parameters
        ----------
        ovhngs : int
            An index in the overhang_files for the file path to use for the assembly

        re : str
            The name of the restriction enzyme to use for the assembly

        ligase : str
            The name of the ligase to use for the assembly

        *args
            See documentation for TraditionalREAssembler parent class

        **kwargs
            See documentation for TraditionalREAssembler parent class
        """
        super(GoldenGateAssembler, self).__init__(re, *args, ligase=ligase, **kwargs)
        with open(self.overhang_files[ovhngs]) as f:
            self.overhangs = json.load(f)

    def primer_extension(self, fragments_pcr, backbone_pcr):
        """
        Determines and creates the Golden Gate primer extensions for a given assembly amplicon in a parts set.


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
