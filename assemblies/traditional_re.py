from assemblies.assembler import Assembler
from Bio.Seq import Seq
from pydna.primer import Primer
from pydna.amplify import pcr, Anneal
from importlib import import_module
import json
from primers.analysis import assembly_thermo
from pydna.dseqrecord import Dseqrecord


# TODO citation
class TraditionalREAssembler(Assembler):
    """
    A class to implement restriction enzyme-based cloning methods. It is mostly used as a parent class
    for other more specific classes like GoldenGateAssembler.
    
    citation: https://www.thermofisher.com/us/en/home/life-science/cloning/cloning-learning-center/invitrogen-school-of-molecular-biology/molecular-cloning/cloning/traditional-cloning-basics.html
    

    Attributes
    ----------
    cloning_type : str
        A simple description of the cloning method

    re1 : Bio.Restriction enzyme
        An enzyme used in a restriction enzyme assembly.

    re2 : optional, Bio.Restriction enzyme 
        A second enzyme used in a restriction enzyme assembly, if needed.

    ligase : optional, str  
        The name of the ligase used for the assembly experiment. Defaults to 'T4-DNA'.

    polymerase : str 
        The name of the polymerase used for the assembly experiment. Defaults to 'T4-DNA'.

    *args
        See documentation for Assembler parent class

    **kwargs
        See documentation for Assembler parent class


    Returns
    -------
    An instance of TraditionalREAssembler


    Methods
    -------
    def add_cutsite(self, site, seq):
        Adds primer extensions for a single RE cutsite to each ends of an insert sequence

    def add_cutsites(self, site1, site2, seq):
        Adds primer extensions for two unique RE cutsites, one on the forward primer one on the reverse primer

    def single_digest_1(self, insert_pcr):
        Performs a single digest using an insert and backbone that share the same RE cutsite and ligation sites

    def single_digest_2A(self, insert_pcr):
        Performs a single digest using an insert and backbone that use two REs, re1 and re2, with compatible ends

    def single_digest_2B(self, insert_pcr):
        Performs a single digest using an insert and backbone using two REs, re1 and re2, with non-compatible ends. 
        One RE for the insert, one RE for the backbone. 

    def double_digest(self, insert_pcr):
        Performs a double digest using an insert and backbone using two REs with non-compatible ends.
        Two REs, re1 and re2, used for the insert and just re2 used on the backbone.
    """
    cloning_type = 'Traditional Restriction Enzyme'

    def __init__(self, re1, *args, re2=None, ligase='T4-DNA', polymerase='T4-DNA', **kwargs):
        """
        Constructs all necessary attributes for the TraditionalREAssembler object.


        Parameters
        ----------
        re1 : Bio.Restriction enzyme
            An enzyme used in a restriction enzyme assembly.

        re2 : optional, Bio.Restriction enzyme 
            A second enzyme used in a restriction enzyme assembly, if needed.

        ligase : optional, str  
            The name of the ligase used for the assembly experiment. Defaults to 'T4-DNA'.

        polymerase : str 
            The name of the polymerase used for the assembly experiment. Defaults to 'T4-DNA'.

        *args
            See documentation for Assembler parent class

        **kwargs
            See documentation for Assembler parent class
        """
        super(TraditionalREAssembler, self).__init__(*args, **kwargs)
        self.ligase = ligase
        self.polymerase = polymerase
        self.re1 = getattr(import_module('Bio.Restriction'), re1)
        if re2:
            self.re2 = getattr(import_module('Bio.Restriction'), re2)

    def blunting(self):
        pass
    
    # TODO use parent class' self.query and self.backbone for design
    # TODO figure out if we can assume backbones already have their cutsites
    # TODO check end compatibility for single_digest_2A before it's called
    # one enzyme
    def add_cutsite(self, site, seq):
        """
        Adds primer extensions for a single RE cutsite to each ends of an insert sequence


        Parameters
        ----------
        site : str
            The RE cutsite for the insert

        seq : pydna Amplicon
            A pydna amplicon object for the insert


        Returns
        -------
        An extended amplicon object for the insert with RE cutsites added
        """
        ext_fwd = Seq(site) + seq.forward_primer
        ext_rev = Seq(site) + seq.reverse_primer

        ext_fwd_p = Primer(ext_fwd, position=seq.forward_primer.position, footprint=seq.forward_primer._fp,id=seq.forward_primer.id)
        ext_rev_p = Primer(ext_rev, position=seq.reverse_primer.position, footprint=seq.reverse_primer._fp,id=seq.reverse_primer.id)

        annealing = Anneal([ext_fwd_p, ext_rev_p], seq.template)
        amp_ext = annealing.products[0]

        return amp_ext

    def add_cutsites(self, site1, site2, seq):
        """
        Adds primer extensions for two RE cutsites, one for each ends of an insert sequence


        Parameters
        ----------
        site1 : str
            The first RE cutsite for the insert

        site2 : str
            The second RE cutsite for the insert

        seq : pydna Amplicon
            A pydna amplicon object for the insert


        Returns
        -------
        An extended amplicon object for the insert with RE cutsites added
        """
        ext_fwd = Seq(site1) + seq.forward_primer
        ext_rev = Seq(site2) + seq.reverse_primer

        ext_fwd_p = Primer(ext_fwd, position=seq.forward_primer.position, footprint=seq.forward_primer._fp,id=seq.forward_primer.id)
        ext_rev_p = Primer(ext_rev, position=seq.reverse_primer.position, footprint=seq.reverse_primer._fp,id=seq.reverse_primer.id)

        annealing = Anneal([ext_fwd_p, ext_rev_p], seq.template)
        amp_ext = annealing.products[0]

        return amp_ext
   
    def single_digest_1(self, insert_pcr):
        """
        Performs a single digest using an insert and backbone that share the same RE cutsite and ligation sites


        Parameters
        ----------
        insert_pcr : A pydna Amplicon 
            An insert Amplicon for the assembly


        Returns
        -------
        A list of the RE-based assembly prepared insert and backbone Amplicon objects
        """
        # forward and reverse primers of both the backbone and insert receive the same cutsite extension
        insert_prep = self.add_cutsite(self.re1.site, insert_pcr)
        backbone_digest = self.digest_backbone(self.re1)
        # TODO digest insert here?
        pass
        return [insert_prep, backbone_digest]

    # two enzymes, compatible ends
    def single_digest_2A(self, insert_pcr):
        """
        Performs a single digest using an insert and backbone that use two REs, re1 and re2, with compatible ends


        Parameters
        ----------
        insert_pcr : A pydna Amplicon 
            An insert Amplicon for the assembly


        Returns
        -------
        A list of the RE-based assembly prepared insert and backbone Amplicon objects
        """
        # insert gets same extensions on fwd and rev for re1
        insert_prep = self.add_cutsite(self.re1.site, insert_pcr)
        # backbone digested by re2
        backbone_digest = self.digest_backbone(self.re2)
        pass
        return [insert_prep, backbone_digest]

    # two enzymes, non-compatible ends
    def single_digest_2B(self, insert_pcr):
        """
        Performs a single digest using an insert and backbone using two REs, re1 and re2, with non-compatible ends. 
        One RE for the insert, one RE for the backbone.


        Parameters
        ----------
        insert_pcr : A pydna Amplicon 
            An insert Amplicon for the assembly


        Returns
        -------
        A list of the RE-based assembly prepared insert and backbone Amplicon objects
        """
        # insert gets same extensions on fwd and rev for re1
        insert_prep = self.add_cutsite(self.re1.site, insert_pcr)
        # backbone digested by re2
        backbone_digest = self.digest_backbone(self.re2)
        # TODO add blunting reactions here
        pass
        return [insert_prep, backbone_digest]

    # two enzymes, non-compatible ends
    def double_digest(self, insert_pcr):

        """
        Performs a double digest using an insert and backbone using two REs with non-compatible ends.
        Two REs, re1 and re2, used for the insert and just re2 used on the backbone.


        Parameters
        ----------
        insert_pcr : A pydna Amplicon 
            An insert Amplicon for the assembly


        Returns
        -------
        A list of the RE-based assembly prepared insert and backbone Amplicon objects
        """
        # insert gets same extensions on fwd and rev for re1
        insert_prep = self.add_cutsites(self.re1.site, self.re2.site, insert_pcr)
        # backbone digested by re2
        backbone_digest = self.digest_backbone(self.re2)
        # TODO add blunting reactions here
        pass
        return [insert_prep, backbone_digest]

    def cut_locations(self, enzyme, record, json=False):
        """
        


        Parameters
        ----------



        Returns
        -------
        """
        site_len = len(enzyme.site)
        fwd_tail_len = len(record.forward_primer.tail)
        enzyme_site_rc = str(Seq(enzyme.site).reverse_complement())

        # Check each 6nt subsequence for a match with the normal or reverse complement
        # of the restriction enzyme site 
        extended_locations = []
        for i in  range(len(record.seq.watson)):
            test_seq = record.seq.watson[i:i + site_len].upper()

            if test_seq == enzyme.site or test_seq == enzyme_site_rc:
                extended_locations.append((i, i + site_len))

        # Adjust the indexes to remove the 5' end extension amount
        locations = [
            (start - fwd_tail_len, end - fwd_tail_len)
            for start, end in extended_locations[1:-1]
        ]

        if json:
            json_locations = [
                {
                    'start': start,
                    'end': end
                }
            for start, end in locations]

            json_ext_locations = [
                {
                    'start': start,
                    'end': end
                }
            for start, end in extended_locations]

            return {'original': json_locations, 'extended': json_ext_locations}
        
        return [extended_locations, extended_locations]

    def cutsite_annotations(self, assembly, enzyme):
        for amplicon in assembly: 
            amplicon.annotations.update({'cuts': amplicon.number_of_cuts(enzyme)})
            amplicon.annotations.update({'cut_locations': self.cut_locations(enzyme, amplicon, json=True)})
            
    def design(self):
        pass 


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
