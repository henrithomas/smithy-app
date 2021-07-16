from assemblies.traditional_re import TraditionalREAssembler

# TODO citation
class GoldenBraidAssembler(TraditionalREAssembler):
    cloning_type = 'GoldenBraid'
    def __init__(self, *args, re1='BsaI', re2='BsmBI'):
        super(GoldenBraidAssembler, self).__init__(re1, *args, re2=re2)