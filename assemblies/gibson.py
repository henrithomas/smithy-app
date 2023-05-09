from assemblies.overlap_extension import OverlapExtensionAssembler

# TODO citation
class GibsonAssembler(OverlapExtensionAssembler):
    cloning_type = 'Gibson'
    exonuclease = 'T5'
    ligase = 'Taq'
    polymerase = 'Phusion'

    def __init__(self, *args, overlap=30, lin_enzyme='EcoRV', **kwargs):
        super(GibsonAssembler, self).__init__(overlap, lin_enzyme, *args, **kwargs)