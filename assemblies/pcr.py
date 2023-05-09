from assemblies.overlap_extension import OverlapExtensionAssembler

class PCRAssembler(OverlapExtensionAssembler):
    cloning_type = 'PCR-SOE'
    polymerase = 'Taq'

    def __init__(self, *args, overlap=50, lin_enzyme='EcoRV', **kwargs):
        super(PCRAssembler, self).__init__(overlap, lin_enzyme, *args, **kwargs)