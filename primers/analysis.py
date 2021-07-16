import primer3 as p3
from pydna.design import primer_design, assembly_fragments
from pydna.dseqrecord import Dseqrecord
from Bio.Seq import Seq
from pydna.tm import tm_default 
from Bio.SeqUtils import MeltingTemp as _mt
from datetime import datetime


class PrimerAnalyzer:
    def __init__(self, mv, dv, dna, dntp, tm_func=None):
        self.mv_conc = mv
        self.dv_conc = dv
        self.dna_conc = dna
        self.dntp_conc = dntp
        self.tm = tm_func
        pass
    
    def tm_old(
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
            dnac1 = self.dna_conc
            dnac2 = self.dna_conc
            Na = self.mv_conc
            Mg = self.dv_conc
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

    def check_gc(seq):
        g_count = seq.count('G') + seq.count('g')
        c_count = seq.count('C') + seq.count('c')
        return (g_count + c_count) / len(seq)

    def get_primers(self, assembly):
        primers = []
        for amplicon in assembly:
            fwd_id = amplicon.name + '-fp'
            rvs_id = amplicon.name + '-rp'
            
            primers.append(
                {
                    'id': fwd_id,
                    'seq': amplicon.forward_primer.seq._data,
                    'footprint': amplicon.forward_primer.footprint._data
                })
                
            primers.append(
                {
                    'id': rvs_id,
                    'seq': amplicon.reverse_primer.seq._data,
                    'footprint': amplicon.reverse_primer.footprint._data
                })
        return primers

    def thermo_analysis(self, primers, mv, dv, dna):
        # a primer = {'id': X, 'footprint': Y, 'seq': Z}
        thermo = []
        for primer in primers:
            thermo.append(
                {
                    'id': primer['id'],
                    'tm_total': self.tm(primer['seq']),
                    'tm_footprint': self.tm(primer['footprint']),
                    'gc': self.check_gc(primer['seq']),
                    'hairpin': p3.calcHairpin(
                                    primer['seq'], 
                                    mv_conc=self.mv_conc,
                                    dv_conc=self.dv_conc, 
                                    dna_conc=self.dna_conc),
                    'homodimer': p3.calcHomodimer(
                                    primer['seq'], 
                                    mv_conc=self.mv_conc,
                                    dv_conc=self.dv_conc, 
                                    dna_conc=self.dna_conc)
                }
            )
        return thermo

# Biosoft Gibbs free energy design guidlines for secondary structures, cross-, and hetero-dimers
def check_hairpin():
    pass

def check_homodimer():
    pass

def check_heterodimer():
    pass

def check_tm():
    pass

def check_tm_mismatch():
    pass

def check_gc(seq):
    g_count = seq.count('G') + seq.count('g')
    c_count = seq.count('C') + seq.count('c')
    return (g_count + c_count) / len(seq)

def primer_thermo(primer, pid, mv, dv, dna, tm_func):
    # takes a pydna primer/amplicon object
    # TODO change gc check to use pydna gc() primer class method
    # str(overhang.seq)
    hp = p3.calcHairpin(primer.seq._data, mv_conc=mv, dv_conc=dv, dna_conc=dna)
    hd = p3.calcHomodimer(primer.seq._data, mv_conc=mv, dv_conc=dv, dna_conc=dna)
    return {
                'id': pid,
                'tm_total': tm_func(primer.seq._data),
                'tm_footprint': tm_func(primer.footprint._data),
                'gc': primer.gc(),
                'hairpin': hp.structure_found,
                'hairpin_tm': hp.tm,
                'hairpin_dg': hp.dg,
                'hairpin_dh': hp.dh,
                'hairpin_ds': hp.ds,
                'homodimer': hd.structure_found,
                'homodimer_tm': hd.tm,
                'homodimer_dg': hd.dg,
                'homodimer_dh': hd.dh,
                'homodimer_ds': hd.ds
            }

def assembly_thermo(assembly, mv, dv, dna, tm_func):
    # takes a pydna assembly object 
    thermo = []
    for amplicon in assembly: 
        # forward primer
        thermo.append(primer_thermo(amplicon.forward_primer, (amplicon.name + '-fp'), mv, dv, dna, tm_func))
        # reverse primer
        thermo.append(primer_thermo(amplicon.reverse_primer, (amplicon.name + '-rp'), mv, dv, dna, tm_func))
    return thermo