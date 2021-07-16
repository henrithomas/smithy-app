from assemblies.traditional_re import TraditionalREAssembler
from Bio.Seq import Seq
from pydna.primer import Primer
from pydna.amplify import pcr
from importlib import import_module
import json

# TODO citation
class MoCloAssembler(TraditionalREAssembler):
    cloning_type = 'MoClo'
    
    # import level 0-2 overhangs and ele linkers
    # import other end linkers when needed in large level 2 constructs
    with open('assemblies\overhangs\moclo\level0.json') as l0, \
        open('assemblies\overhangs\moclo\level1.json') as l1, \
        open('assemblies\overhangs\moclo\endlinkers\ele.json') as ele:
        l0_ovhngs = json.load(l0)
        l1_ovhngs = json.load(l1)
        ele_ovhngs = json.load(ele)
    
    def __init__(self, *args, re1='BsaI', re2='BpiI', re3='Esp3I'):
        super(MoCloAssembler, self).__init__(re1, *args, re2=re2)
        self.re3 = getattr(import_module('Bio.Restriction'), re3)
    
    def module_prep(self, seq, part='P'):
        # TODO add BpiI cutsite with proper overhangs for the part type
        pass

    def module(self, seq, dest, part='P'):
        # TODO either cut out a part from a plasmid with BsaI or add BsaI sites 
        # and the proper overhangs for the part type
        pass

    def tu(self, modules, dest):
        # TODO take modules and assemble them together based on overhangs,
        # add proper BpiI sites and overhangs or ignore if already in an
        # l1 vecre with the BpiI sites and overhangs
        pass

    def multigene(self, t_units, dest, size=2, sublevel=1):
        # TODO iterate over tu's to assemble them into an l2 vector, 
        # use the proper and necessary end linkers for level 2i-X 
        # assemblies
        # TODO user the size parameter to determine which and how many end
        # linkers to use  
        pass