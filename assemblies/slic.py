from assemblies.overlap_extension import OverlapExtensionAssembler

class SLICAssembler(OverlapExtensionAssembler):
    cloning_type = 'SLIC'
    exonuclease = 'T4-DNA'
    recombination_enz = 'RecA'
    polymerase = 'T4-DNA'

    def __init__(self, *args, overlap=40, lin_enzyme='EcoRV', **kwargs):
        super(SLICAssembler, self).__init__(overlap, lin_enzyme, *args, **kwargs)