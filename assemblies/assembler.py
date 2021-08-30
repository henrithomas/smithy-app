from Bio import SeqIO
from pydna.dseqrecord import Dseqrecord
from pydna.amplify import pcr
from pydna.design import primer_design, assembly_fragments
from pydna.assembly import Assembly
from importlib import import_module
from datetime import datetime
from fragment.tree import FragmentTree
from query.blaster import Blaster
from Bio.SeqUtils import MeltingTemp as _mt
from primers.analysis import assembly_thermo
from Bio.Seq import Seq
from pydna.seqrecord import SeqRecord 
from pydna.primer import Primer
from pydna.amplify import pcr


class Assembler:
    """
    A base parent class for other assembly classes that will specify cloning methods.


    Attributes
    ----------
    cloning_type : str
        A simple description of the cloning method

    mv_conc : float
        Monovalent ion concentration
    
    dv_conc : float
        Divalent ion concentration
    
    dna_conc : float
        DNA concentration
    
    dntp_conc : float
        dNTP concentration
    
    tm : float
        Melting temperature
    
    backbone_file : string
        Path to the backbone fasta file
    
    query_file : string
        Path to the insert query file
   
    dbs : list
        List of strings for the BLAST databases to use
    
    min_frag : int, optional
        Minimum nucleotide size of alignments to use from BLAST databases
    
    max_frag : int, optional
        Maximum nucleotide size of alignments to use from BLAST databases
    
    min_synth : int, optional
        Minimum nucleotide size of synthetic fragments to design for assemblies
    
    max_synth : int, optional
        Maximum nucleotide size of synthetic fragments to design for assemblies
    
    max_seqs : int, optional
        Maximum count of alignments to return from BLAST queries
    
    query_record : Dseqrecord
        A pydna Dseqrecord of the insert query sequence built from query_file
    
    backbone : Dseqrecord
        A pydna Dseqrecord of the backbone sequence built from backbone_file
    
    solution_tree : FragmentTree
        An instance of FragmentTree that holds all information for assembly solutions 
    
    insert : str
        (unused) The sequence for a solution for an assembly 


    Returns
    -------
    An instance of Assembler

    
    Methods
    -------
    def tm_custom(
        self,
        seq,
        check=True,
        strict=True,
        c_seq=None,
        shift=0,
        nn_table=_mt.DNA_NN4,  # DNA_NN4: values from SantaLucia & Hicks (2004)
        tmm_table=None,
        imm_table=None,
        de_table=None,
        dnac1=50, 
        dnac2=50,  
        selfcomp=False,
        Na=50,
        K=0,
        Tris=75.0,  # We use the 10X Taq Buffer with (NH4)2SO4 (above)
        Mg=1.5,  # 1.5 mM Mg2+ is often seen in modern protocols
        dNTPs=0.8,  # I assume 200 µM of each dNTP
        saltcorr=7,  # Tm = 81.5 + 0.41(%GC) - 600/N + 16.6 x log[Na+]
        func=_mt.Tm_NN,  # Used by Primer3Plus to calculate the product Tm.
    ):
        Uses parameters from instantiation to build a custom pydna melting temperature function.

    def parts_csv(self, parts, file_name):
        Exports the assembly parts data to the chosen file in csv format.
        
    def primers_csv(self, parts, file_name):
        Exports the assembly primer data to the choses file in csv format.

    def assembly_construct_write(self, frags, file_name):
        Exports the full assembly construct sequence. 

    def assembly_output(self, frags, b_bone, assembly):
        A convenience function to bundle the parts_csv, primers_csv, and assembly_construct_write functions.

    def digest_backbone(self, enzyme):
        Uses a pydna enzyme to cut a pydna Dseqrecord backbone sequence.

    def get_solution(self, s):
        Collects a Dseqrecord list of assembly fragments from a solution tree.
    
    def primer_complement(self, fragments, backbone):
        Designs basic primers for all fragments and the backbone for the assembly.

    def primer_verification(self, assembly):
        Runs thermodynamic analysis on all sequences in an assembly using the assembly_thermo function from primers.analysis.

    def query(self):
        Creates a Blaster instance using attributes fron the current assembler object and runs BLAST queries using Blaster.

    def solution_building(self, query_results):
        Creates a Fragment tree using Blaster results and then initializes the assembler's solution_tree.

    def annotations(self, assembly, nodes):
        Adds annotations to the Dseqrecords in the assembly from the BLAST alignment data in nodes.
    """
    
    cloning_type = 'Abstract'

    def __init__(self, mv, dv, dna, dntp, tm, backbone_file, query_file, db_list, min_frag=0, max_frag=1000, min_synth=0, max_synth=1000, multi_query=False):
        """
        Constructs all necessary attributes for the Assembler object.


        Parameters
        ----------
            mv_conc : float
                Monovalent ion concentration
            
            dv_conc : float
                Divalent ion concentration
            
            dna_conc : float
                DNA concentration
            
            dntp_conc : float
                dNTP concentration
            
            tm : float
                Melting temperature
            
            backbone_file : string
                Path to the backbone fasta file
            
            query_file : string
                Path to the insert query file
        
            dbs : list
                List of strings for the BLAST databases to use
            
            min_frag : int, optional
                Minimum nucleotide size of alignments to use from BLAST databases
            
            max_frag : int, optional
                Maximum nucleotide size of alignments to use from BLAST databases
            
            min_synth : int, optional
                Minimum nucleotide size of synthetic fragments to design for assemblies
            
            max_synth : int, optional
                Maximum nucleotide size of synthetic fragments to design for assemblies
            
            max_seqs : int, optional
                Maximum count of alignments to return from BLAST queries
            
            query_record : Dseqrecord
                A pydna Dseqrecord of the insert query sequence built from query_file
            
            backbone : Dseqrecord
                A pydna Dseqrecord of the backbone sequence built from backbone_file
            
            solution_tree : FragmentTree
                An instance of FragmentTree that holds all information for assembly solutions 
            
            insert : str
                (unused) The sequence for a solution for an assembly 
        """

        self.mv_conc = mv
        self.dv_conc = dv
        self.dna_conc = dna
        self.dntp_conc = dntp
        self.tm = tm
        self.backbone_file = backbone_file   
        self.query_file = query_file    
        self.dbs = db_list  
        self.min_frag = min_frag
        self.max_frag = max_frag
        self.min_synth = min_synth
        self.max_synth = max_synth
        self.max_seqs = 1000
        self.backbone = Dseqrecord(SeqIO.read(backbone_file, 'fasta'), circular=False)
        self.solution_tree = None
        self.insert = None
        self.multi_query = multi_query
        self.max_solutions = 100
        self.query_record = None 
        # self.query_records = []

        if multi_query:
            with open(query_file) as query_file:
                self.query_records = [Dseqrecord(record) for record in SeqIO.parse(query_file, 'fasta')] 
        else:
            self.query_record = Dseqrecord(SeqIO.read(query_file, 'fasta'), circular=False)
        
    def tm_custom(
        self,
        seq,
        check=True,
        strict=True,
        c_seq=None,
        shift=0,
        nn_table=_mt.DNA_NN4,  # DNA_NN4: values from SantaLucia & Hicks (2004)
        tmm_table=None,
        imm_table=None,
        de_table=None,
        dnac1=50, 
        dnac2=50,  
        selfcomp=False,
        Na=50,
        K=0,
        Tris=75.0,  # We use the 10X Taq Buffer with (NH4)2SO4 (above)
        Mg=1.5,  # 1.5 mM Mg2+ is often seen in modern protocols
        dNTPs=0.8,  # I assume 200 µM of each dNTP
        saltcorr=7,  # Tm = 81.5 + 0.41(%GC) - 600/N + 16.6 x log[Na+]
        func=_mt.Tm_NN,  # Used by Primer3Plus to calculate the product Tm.
    ):
        """
        Constructs a custom pydna melting temperature for the assembler object.


        Parameters
        ----------
        See pydna docs.


        Returns
        -------
        A custom pydna melting temperature function to use throughout the assembler object.
        """
        dnac1 = self.dna_conc
        dnac2 = self.dna_conc
        Na = self.mv_conc
        Mg = self.mv_conc
        dNTPs = self.dntp_conc
        return func(
            seq,
            check=check,
            strict=strict,
            c_seq=c_seq,
            shift=shift,
            nn_table=nn_table,
            tmm_table=tmm_table,
            imm_table=imm_table,
            de_table=de_table,
            dnac1=dnac1,
            dnac2=dnac2,
            selfcomp=selfcomp,
            Na=Na,
            K=K,
            Tris=Tris,
            Mg=Mg,
            dNTPs=dNTPs,
            saltcorr=saltcorr,
        )

    def parts_csv(self, parts, file_name):
        """
        Exports the assembly parts data to the chosen file in csv format.


        Parameters
        ----------
        parts : list
            List of Dseqrecords objects for the assembly
            
        file_name : str
            Name of the csv file to export


        Returns
        -------
        None
        """
        import csv
        fields = ['id', 'db','length', 'length_ext', 'seq', 'seq_ext']
        csv_list = []

        for part in parts:
            csv_list.append({
                'id': part.name, 
                'db': part.annotations['db'],
                'length': part.template.seq.length, 
                'length_ext': part.seq.length, 
                'seq': part.template.seq.watson, 
                'seq_ext': part.seq.watson
            })
        
        with open(file_name, 'w', newline='') as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=fields, restval='NONE')
            writer.writeheader()
            writer.writerows(csv_list)

    def primers_csv(self, parts, file_name):
        """
        Exports the assembly primers to the chosen file in csv format.


        Parameters
        ----------
        parts : list
            List of Dseqrecords objects

        file_name : str
            Name of the csv file to export


        Returns
        -------        
        None
        """
        # this will become more complicated when JSON structure is added
        # TODO add JSON structure for parts
        import csv
        fields = ['id', 'primer_type', 'sequence', 'footprint', 'tail']
        csv_list = []

        for part in parts:
            csv_list.extend([{
                    'id': f'{part.name}-fp',
                    'primer_type': 'f',
                    'sequence': part.forward_primer.seq,
                    'footprint': part.forward_primer.footprint,
                    'tail': part.forward_primer.tail,
                },
                {
                    'id': f'{part.name}-rp',
                    'primer_type': 'r',
                    'sequence': part.reverse_primer.seq,
                    'footprint': part.reverse_primer.footprint,
                    'tail': part.reverse_primer.tail,
                }])

        with open(file_name, 'w', newline='') as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=fields, restval='NONE')
            writer.writeheader()
            writer.writerows(csv_list)

    def assembly_construct_write(self, frags, file_name):
        """
        Exports the full assembly sequence, backbone and insert, to the chosen file in csv format.


        Parameters
        ----------
        frags : list 
            List of strings of all assembly sequences

        file_name : str
            Name of the file to export


        Returns
        -------        
        None
        """
        # this will become more complicated when JSON structure is added
        # TODO add JSON structure for construct
        construct = ''.join([str(frag.seq.watson) for frag in frags])
        with open(file_name, 'w') as file:
            file.write(construct)

    def assembly_output(self, frags, b_bone, assembly):
        """
        A convenience function to bundle the parts_csv, primers_csv, and assembly_construct_write functions.


        Parameters
        ----------
        frags : list
            List of strings of all assembly sequences

        b_bone : Dseqrecord
            Dseqrecord of the backbone sequence

        assembly : list
            List of Dseqrecords objects for the assembly


        Returns
        -------        
        None
        """
        # format the fragments, extension primers, and digest['backbone'] sequences here
        # maybe write out two files here, one for PCR without extension and another with extension
        # combine the digested backbone with the original fragments into a single sequence
        # create a dated folder to put all these files in so that they are not just in the cwd
        frags.insert(0, b_bone)
        self.assembly_csv_write(frags, '.\\results\\primers_ne.csv')
        self.assembly_csv_write(assembly, '.\\results\\primers_e.csv')
        self.assembly_construct_write(frags, '.\\results\\construct.txt')

    # TODO change method to use pydna linearize function for exceptions
    # TODO add more thorough checks on linearizing the backbone with given enzyme
    def digest_backbone(self, enzyme):
        """
        Uses a pydna enzyme to cut a pydna Dseqrecord backbone sequence.        


        Parameters
        ----------
        enzyme : Bio.Restriction enzyme
            An enzyme instance from Bio.Restriction 


        Returns
        -------        
        A dict with the undigested or digested backbone and the number of cuts.
        """
        results = {
            'cuts': self.backbone.number_of_cuts(enzyme),
            'backbone': self.backbone,
        }
        if results['cuts'] == 1:
            results['backbone'] = self.backbone.cut(enzyme)[0]
        return results

    def get_solution(self, s):
        """
        Collects a Dseqrecord list of assembly fragments from a solution tree.


        Parameters
        ----------
        s : int
            The solution index to retrieve from the solution tree 


        Returns
        -------
        A Dseqrecord fragment list for the sequences of a solution        
        """
        tree_solution = self.solution_tree.solution_seqs(s)
        fragments = [Dseqrecord(record) for record in tree_solution]
        return fragments

    def primer_complement(self, fragments, backbone):
        """
        Designs basic non-extension primers for all fragments and the backbone for the assembly.


        Parameters
        ----------
        fragments : list
            List of Dseqrecords objects for the assembly

        backbone : Dseqrecord
            Dseqrecord of the backbone sequence


        Returns
        -------
        A tuple of a list of fragment Amplicons and a backbone Amplicon after non-extension primer design 
        """
        fragments_pcr = [primer_design(fragment, target_tm=self.tm, tm_func=self.tm_custom) for fragment in fragments]
        backbone_pcr = primer_design(backbone, target_tm=self.tm, tm_func=self.tm_custom)
        
        return fragments_pcr, backbone_pcr

    def query(self):
        """
        Creates a Blaster instance using attributes fron the current assembler object and runs BLAST queries using Blaster.


        Parameters
        ----------
        None


        Returns
        -------        
        Blaster query results from the BLAST query and stderr messages from run_blastn()
        """
        blaster = Blaster(self.max_seqs, self.dbs, self.min_frag, self.max_frag)
        queries = blaster.queries(self.query_file)
        query_results, stderr = blaster.run_blastn(queries)
        return query_results, stderr

    # TODO separate into two methods to first build the frag tree, then select a solution
    def solution_building(self, query_results):
        """
        Creates a Fragment tree using Blaster results and then initializes the assembler's solution_tree.


        Parameters
        ----------
        query_results : list
            A long list of BLAST alignments


        Returns
        -------
        None        
        """
        sol_tree = FragmentTree(self.query_record.seq.watson, nodes=[], 
                        query_len=self.query_record.seq.length, 
                        min_synth=self.min_synth, 
                        max_synth=self.max_synth)
        sol_tree.build(query_results)
        self.solution_tree = sol_tree

    def annotations(self, assembly, nodes, space=0):
        """
        Adds annotations to the Dseqrecords in the assembly from the BLAST alignment data in nodes.


        Parameters
        ----------
        assembly : list 
            A list of Amplicons from the assembly design

        nodes : list
            A list of FragmentTree assembly FragmentNodes


        Returns
        -------        
        A new Amplicon list for the assembly with annotations for names and database sources
        """
        new_assembly = assembly.copy()
        backbone_start = self.query_record.seq.length + space
        backbone_end = backbone_start + self.backbone.seq.length

        if self.multi_query:
            for i, node in enumerate(nodes):
                seq_lengths = sum([record.template.seq.length for record in assembly[:i]])
                node.data.hsps[0].query_start = space * (i + 1) + seq_lengths
                node.data.hsps[0].query_end = node.data.hsps[0].query_start + assembly[i].template.seq.length

        data = [[
                node.node_id, 
                node.db,
                node.data.hsps[0].query_start,
                node.data.hsps[0].query_end,
                node.data.hsps[0].sbjct_start,
                node.data.hsps[0].sbjct_end
                ] for node in nodes]
        data.append([self.backbone.name, 'NONE', backbone_start, backbone_end, 0, 0])

        for i, amplicon in enumerate(new_assembly):
            amplicon.name = data[i][0]
            amplicon.annotations.update({'db': data[i][1]})
            amplicon.annotations.update({'query_start': data[i][2]})
            amplicon.annotations.update({'query_end': data[i][3]})
            amplicon.annotations.update({'subject_start': data[i][4]})
            amplicon.annotations.update({'subject_end': data[i][5]})
            
        return new_assembly

    def run_multi_query(self):
        lengths = [record.seq.length for record in self.query_records]
        num_fragments = len(self.query_records)
        blaster = Blaster(self.max_seqs, self.dbs, self.min_frag, self.max_frag)
        queries = blaster.queries(self.query_file)
        query_results, stderr = blaster.run_multi_blastn(queries, lengths, num_fragments)
        return query_results, stderr

    def multi_query_solution_building(self, query_results):
        sequences = [record.seq.watson for record in self.query_records]
        sol_tree = FragmentTree('')
        sol_tree.multi_query_blast_input(query_results, sequences)
        sol_tree.build_multi_query_solutions(self.max_solutions)
        self.solution_tree = sol_tree

    def get_multi_query_solution(self, s):
        tree_solution = self.solution_tree.multi_query_solution_seqs(s)
        fragments = [Dseqrecord(record) for record in tree_solution]
        return fragments 