from pydna.design import primer_design
from pydna.dseqrecord import Dseqrecord
from assemblies.traditional_re import TraditionalREAssembler
from Bio.Seq import Seq
from pydna.primer import Primer
from pydna.amplify import Anneal
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

    def __init__(self, *args, ovhngs=0, re='BsaI', ligase='T4-DNA', scarless=True, **kwargs):
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

        scarless : bool
            Boolean value that determines if scarless primer design is used

        *args
            See documentation for TraditionalREAssembler parent class

        **kwargs
            See documentation for TraditionalREAssembler parent class
        """
        super(GoldenGateAssembler, self).__init__(re, *args, ligase=ligase, **kwargs)
        with open(self.overhang_files[ovhngs]) as f:
            self.overhangs = json.load(f)

        self.scarless = scarless

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

    def primer_complement_scarless(self, fragments, backbone, interval=10):
        """
        Searches an interval on each side of a junction between two neighboring fragments for a 
        four-base sequence matching an overhang sequence in the high-fidelity overhang set (self.overhangs).
        
        When an overhang is found in either the current fragment or it's right neighbor, the indexes of the 
        overhang are recorded. These indexes are used for extracting sequences that extend either the current 
        or neighboring fragment with the overhang at the end. If no overhang sequence is found, the search 
        defaults to using the last four bases of the current fragment for the overhang sequence. Indexes of the 
        overhangs are also used for properly PCR'ing each fragment if/when it must be trimmed for the neighboring
        fragment to extend into the trimmed suences (saved to pcr_intervals). 

        The extension sequences are extended further with the BsaI cutsite ('GGTCTCN'). 
        
        Finally, this designs non-extension primers for all fragments and the backbone for the assembly using the 
        values in pcr_intervals.

        Example:

        cutsite = 'GGTCTCN'
        interval = 10
        test_seq = current[-interval:] + right[:interval]

                       j=0        9          19
                       |          |          |
                      |---[----]---|----------| <--test_seq
                           3    7   
        current-----------[----]---|------------------right
                           ^^^^
                           overhang

        Overhang found at j = 3 
        ovhng_start = -(interval - j) = -7 
        ovhng_end = -(interval - (j + 4)) = -3

        current = current[:ovhng_end]
        right += current[-ovhng_start:]
        current_extension = cutsite
        right_extension = cutsite + current[-ovhng_start:]

                        ovhng_start
                          |
              |~~cutsite~~[----]---|------------------right

        current-----------[----]~~cutsite~~|
                               |
                            ovhng_end 

        Parameters
        ----------
        fragments : list
            List of Dseqrecords objects for the assembly

        backbone : Dseqrecord
            Dseqrecord of the backbone sequence

        Returns
        -------
        A list of non-extension amplicons for the assembly (fragments_pcr) and the extensions to add to the primers
        of fragments_pcr that facilitate assembly (extensions).
        """
        fragments.append(backbone)
        pcr_intervals = [[0, record.seq.length] for record in fragments]
        overhang_set = set(self.overhangs['seq'])
        extensions = [['',''] for _ in range(len(fragments))]

        for i in range(len(fragments)):
            found = False
            i_next = 0 if i == len(fragments) - 1 else i + 1
            test_seq = fragments[i].seq.watson[-interval:] + fragments[i_next].seq.watson[:interval]

            # search the sequence range for an overhang sequence in the 
            # high fidelity set
            for j in range(len(test_seq)):
                test_ovhng = test_seq[j:j+4].upper()

                if test_ovhng in overhang_set:
                    found = True
                    overhang_set.remove(test_ovhng)
                    break

            # if an overhang from the set is not found, choose a default 
            # overhang of the last four bases of the left fragment 
            if not found:
                j = interval - 4

            # set the overhang bounds if it is found in the neighbor fragment (i_next)
            if j > interval - 1:
                ovhng_start = j - interval
                ovhng_end = (j + 4) - interval

                pcr_intervals[i_next][0] = ovhng_start

                extensions[i][1] = self.re1.site + 'N' + fragments[i_next].seq.crick[-ovhng_end:]
                extensions[i_next][0] = self.re1.site + 'N'
            # set the overhang bounds if it is found in the current fragment (i)
            else:
                ovhng_start = -(interval - j)
                ovhng_end = -(interval - (j + 4))

                if ovhng_end != 0:
                    pcr_intervals[i][1] = ovhng_end
                
                extensions[i][1] = self.re1.site + 'N' 
                extensions[i_next][0] = self.re1.site + 'N' + fragments[i].seq.watson[ovhng_start:]

        fragments_pcr = [
            primer_design(
                fragment[pcr_intervals[i][0]:pcr_intervals[i][1]], 
                target_tm=self.tm, 
                tm_func=self.tm_custom
            ) 
            for i, fragment in enumerate(fragments)
        ]

        return fragments_pcr, extensions

    def primer_extension_scarless(self, fragments_pcr, extensions):
        """
        Extends each amplicon in fragments_pcr using extensions which contains the proper extension sequences to use for the 
        forward and reverse primers of each amplicon that facilite proper scarless assembly.  


        Parameters
        ----------
        fragments_pcr : list of pydna Amplicons  
            A list of pydna Amplicons used for a given assembly solution

        extensions : list of lists 
            Extension sequences for forward and reverse primer sequence extensions, [[extension_fwd, extension_rvs],...] 

        Returns
        -------
        A list of new Amplicon objects, including fragments and backbone, that have gone through primer extension design and 
        extended accordingly
        """
        assembly = []

        # add proper extensions to each primer, either just the cutsite or additional
        # bases from a neighboring fragment 
        for extension, amp in zip(extensions, fragments_pcr):
            amp_extended = self.add_cutsites(extension[0], extension[1], amp)
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
        else:
            fragments = self.get_solution(solution)
            nodes = self.solution_tree.solution_nodes(solution)

        if self.scarless:
            # create assembly primer complements for backbone and fragments
            fragments_pcr, extensions = self.primer_complement_scarless(fragments, self.backbone, interval=20) 
            # create primer extensions that facilitate scarless assembly
            assembly = self.primer_extension_scarless(fragments_pcr, extensions)
        else:
            # create assembly primer complements for backbone and fragments
            fragments_pcr, backbone_pcr = self.primer_complement(fragments, self.backbone)
            # create primer extensions
            assembly = self.primer_extension(fragments_pcr, backbone_pcr)

        if self.multi_query:
            frag_seqs = [record.seq.watson[7:-11] for record in assembly[:-1]]
            self.query_record = Dseqrecord(''.join(frag_seqs))
            pass 

        # add simple annotations
        assembly = self.annotations(assembly, nodes, space=4)
        self.cutsite_annotations(assembly, self.re1)
        
        # run primer thermo analysis
        assembly = assembly_thermo(assembly, self.mv_conc, self.dv_conc, self.dna_conc, self.tm_custom)

        return assembly, nodes
