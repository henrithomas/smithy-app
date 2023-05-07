from contextlib import AsyncExitStack
import primer3 as p3
from pydna.design import primer_design, assembly_fragments
from pydna.dseqrecord import Dseqrecord
from Bio.Seq import Seq
from pydna.tm import tm_default 
from Bio.SeqUtils import MeltingTemp as _mt
from datetime import datetime

def gc_content(primer):
    count = primer.count('c')
    count = count + primer.count('g')
    return round((count / len(primer)) * 100, 2) 

def primer_thermo(primer, name, mv_conc, dv_conc, dna_conc, tm_func):
    """
    Calculates thermodynamic analysis on a single primer within an assembly solution


    Parameters
    ----------
    primer : pydna Primer 
        The primer to perform thermodynamic analysis on 

    name : str
        An identifying label for the primer

    mv_conc : float
        Monovalent ion concentration

    dv_conc : float
        Divalent ion concentration

    dna_conc : float
        DNA concentration

    tm_func : pydna Tm function
        The melting temperature function to use for the analysis


    Returns
    -------
    A dict containing thermodynamic for a primer
    """
    # takes a pydna primer/amplicon object
    hp = p3.calcHairpin(primer.seq._data, mv_conc=mv_conc, dv_conc=dv_conc, dna_conc=dna_conc)
    hd = p3.calcHomodimer(primer.seq._data, mv_conc=mv_conc, dv_conc=dv_conc, dna_conc=dna_conc)
    return {
                'id': name,
                'tm_total': tm_func(primer.seq._data),
                'tm_footprint': tm_func(primer.footprint._data),
                'gc': gc_content(primer),
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

def assembly_thermo(assembly, mv_conc, dv_conc, dna_conc, tm_func):
    """
    Runs thermodynamic analysis over forward and reverse primers for each part in an assembly solution


    Parameters
    ----------
    assembly : list
        An assembly solution set of pydna Amplicon objects

    mv_conc : float
        Monovalent ion concentration

    dv_conc : float
        Divalent ion concentration

    dna_conc : float
        DNA concentration

    tm_func : pydna Tm function
        The melting temperature function to use for the analysis


    Returns
    -------
    A new assembly set with annotations for each part's forward and reverse primer's thermodynamic analysis
    """
    # takes a pydna assembly object 
    assembly_copy = assembly.copy()
    for amplicon in assembly: 
        # forward primer
        thermo_fwd = {'forward_primer': primer_thermo(amplicon.forward_primer, (amplicon.name + '-fwd'), mv_conc, dv_conc, dna_conc, tm_func)} 
        amplicon.annotations.update(thermo_fwd)
        # reverse primer
        thermo_rvs = {'reverse_primer': primer_thermo(amplicon.reverse_primer, (amplicon.name + '-rvs'), mv_conc, dv_conc, dna_conc, tm_func)} 
        amplicon.annotations.update(thermo_rvs)

    return assembly_copy