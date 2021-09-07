from assemblies.overlap_extension import OverlapExtensionAssembler

# TODO citation
class GibsonAssembler(OverlapExtensionAssembler):
    cloning_type = 'Gibson'
    exonuclease = 'T5'
    ligase = 'Taq'
    polymerase = 'Phusion'

    def __init__(self, *args, overlap=30, lin_enzyme='EcoRV', **kwargs):
        super(GibsonAssembler, self).__init__(overlap, lin_enzyme, *args, **kwargs)


class SLICAssembler(OverlapExtensionAssembler):
    cloning_type = 'SLIC'
    exonuclease = 'T4-DNA'
    recombination_enz = 'RecA'
    polymerase = 'T4-DNA'

    def __init__(self, *args, overlap=40, lin_enzyme='EcoRV', **kwargs):
        super(SLICAssembler, self).__init__(overlap, lin_enzyme, *args, **kwargs)


class PCRAssembler(OverlapExtensionAssembler):
    cloning_type = 'PCR-SOE'
    polymerase = 'Taq'

    def __init__(self, *args, overlap=50, lin_enzyme='EcoRV', **kwargs):
        super(PCRAssembler, self).__init__(overlap, lin_enzyme, *args, **kwargs)
