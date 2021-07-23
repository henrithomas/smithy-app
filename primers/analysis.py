from contextlib import AsyncExitStack
import primer3 as p3
from pydna.design import primer_design, assembly_fragments
from pydna.dseqrecord import Dseqrecord
from Bio.Seq import Seq
from pydna.tm import tm_default 
from Bio.SeqUtils import MeltingTemp as _mt
from datetime import datetime


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
    assembly_copy = assembly.copy()
    for amplicon in assembly: 
        # forward primer
        thermo_fwd = {'forward_primer': primer_thermo(amplicon.forward_primer, (amplicon.name + '-fwd'), mv, dv, dna, tm_func)} 
        amplicon.annotations.update(thermo_fwd)
        # reverse primer
        thermo_rvs = {'reverse_primer': primer_thermo(amplicon.reverse_primer, (amplicon.name + '-rvs'), mv, dv, dna, tm_func)} 
        amplicon.annotations.update(thermo_rvs)

    return assembly_copy