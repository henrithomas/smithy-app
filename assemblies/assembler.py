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
    cloning_type = 'Abstract'

    def __init__(self, mv, dv, dna, dntp, tm, bb_file, q_file, db_list, min_frag=0, max_frag=0, min_synth=0, max_synth=0):
        self.mv_conc = mv
        self.dv_conc = dv
        self.dna_conc = dna
        self.dntp_conc = dntp
        self.tm = tm
        self.backbone_file = bb_file    # 'C:\\Users\\mt200\\Desktop\\projects\\thesis\\input\\sequences\\backbone.fasta'
        self.query_file = q_file    # 'C:\\Users\\mt200\\Desktop\\projects\\thesis\\input\\queries\\rpr1-dcas9mxi1-gfp-adh1.fasta'
        self.dbs = db_list  # ['addgene', 'igem', 'dnasu']
        self.min_frag = min_frag
        self.max_frag = max_frag
        self.min_synth = min_synth
        self.max_synth = max_synth
        self.max_seqs = 1000
        self.query_record = Dseqrecord(SeqIO.read(q_file, 'fasta'), circular=False)
        self.backbone = Dseqrecord(SeqIO.read(bb_file, 'fasta'), circular=False)
        self.solution_tree = None
        self.insert = None
        

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
        # change these params to the assemler values
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
        # this will become more complicated when JSON structure is added
        # TODO add JSON structure for construct
        construct = ''.join([str(frag.seq.watson) for frag in frags])
        with open(file_name, 'w') as file:
            file.write(construct)


    def assembly_output(self, frags, b_bone, assembly):
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
        results = {
            'cuts': self.backbone.number_of_cuts(enzyme),
            'backbone': self.backbone,
        }
        if results['cuts'] == 1:
            results['backbone'] = self.backbone.cut(enzyme)[0]
        return results

    def get_solution(self, s):
        tree_solution = self.solution_tree.solution_seqs(s)
        fragments = [Dseqrecord(record) for record in tree_solution]
        return fragments

    def primer_complement(self, fragments, backbone):
        fragments_pcr = [primer_design(fragment, target_tm=self.tm) for fragment in fragments]
        backbone_pcr = primer_design(backbone, target_tm=self.tm)
        # assembly_set = self.prepare_assembly(fragments_pcr, backbone_pcr, overlap)
        return fragments_pcr, backbone_pcr

    def primer_verification(self, assembly):
        thermo_set = assembly_thermo(assembly, self.mv_conc, self.dv_conc, self.dna_concdna, self.tm_custom)
        return thermo_set

    def query(self):
        blaster = Blaster(self.max_seqs, self.dbs, self.min_frag, self.max_frag)
        queries = blaster.queries(self.query_file)
        query_results, stderr = blaster.run_blastn(queries)
        return query_results, stderr

    # TODO separate into two methods to first build the frag tree, then select a solution
    def solution_building(self, query_results):
        sol_tree = FragmentTree(self.query_record.seq.watson, nodes=[], 
                        query_len=self.query_record.seq.length, 
                        min_synth=self.min_synth, 
                        max_synth=self.max_synth)
        sol_tree.build(query_results)
        self.solution_tree = sol_tree

    def annotations(self, assembly, nodes):
        new_assembly = assembly.copy()
        data = [[node.node_id, node.db] for node in nodes]
        data.insert(0, [self.backbone.name, 'NONE'])
        for i, amplicon in enumerate(new_assembly):
            amplicon.name = data[i][0]
            amplicon.annotations.update({'db': data[i][1]})
        return new_assembly
