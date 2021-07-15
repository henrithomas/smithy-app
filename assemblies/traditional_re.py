from assemblies.assembler import Assembler
from Bio.Seq import Seq
from pydna.primer import Primer
from pydna.amplify import pcr
from importlib import import_module
import json

# TODO citation
class TraditionalREAssembler(Assembler):
    cloning_type = 'Trad. Restriction Enzyme'

    def __init__(self, re1, *args, re2=None, ligase='T4-DNA', polymerase='T4-DNA', **kwargs):
        super(TraditionalREAssembler, self).__init__(*args, **kwargs)
        self.ligase = ligase
        self.polymerase = polymerase
        self.re1 = getattr(import_module('Bio.Restriction'), re1)
        if re2:
            self.re2 = getattr(import_module('Bio.Restriction'), re2)

    def blunting(self):
        pass
    
    # citation: https://www.thermofisher.com/us/en/home/life-science/cloning/cloning-learning-center/invitrogen-school-of-molecular-biology/molecular-cloning/cloning/traditional-cloning-basics.html#vector
    # TODO use parent class' self.query and self.backbone for design
    # TODO figure out if we can assume backbones already have their cutsites
    # TODO check end compatibility for single_digest_2A before it's called
    # one enzyme
    def add_cutsite(self, site, seq):
        ext_fwd = Seq(site) + seq.forward_primer
        ext_rev = Seq(site) + seq.reverse_primer

        ext_fwd_p = Primer(ext_fwd, position=seq.forward_primer.position, footprint=seq.forward_primer._fp,id=seq.forward_primer.id)
        ext_rev_p = Primer(ext_rev, position=seq.reverse_primer.position, footprint=seq.reverse_primer._fp,id=seq.reverse_primer.id)

        return pcr(ext_fwd_p, ext_rev_p, seq.template)

    def add_cutsites(self, site1, site2, seq):
        ext_fwd = Seq(site1) + seq.forward_primer
        ext_rev = Seq(site2) + seq.reverse_primer

        ext_fwd_p = Primer(ext_fwd, position=seq.forward_primer.position, footprint=seq.forward_primer._fp,id=seq.forward_primer.id)
        ext_rev_p = Primer(ext_rev, position=seq.reverse_primer.position, footprint=seq.reverse_primer._fp,id=seq.reverse_primer.id)

        return pcr(ext_fwd_p, ext_rev_p, seq.template)

    
    def single_digest_1(self, insert_pcr):
        # forward and reverse primers of both the backbone and insert receive the same cutsite extension
        insert_prep = self.add_cutsite(self.re1.site, insert_pcr)
        backbone_digest = self.digest_backbone(self.re1)
        # TODO digest insert here?
        pass
        return [insert_prep, backbone_digest]

    # two enzymes, compatible ends
    def single_digest_2A(self, insert_pcr):
        # insert gets same extensions on fwd and rev for re1
        insert_prep = self.add_cutsite(self.re1.site, insert_pcr)
        # backbone digested by re2
        backbone_digest = self.digest_backbone(self.re2)
        pass
        return [insert_prep, backbone_digest]

    # two enzymes, non-compatible ends
    def single_digest_2B(self, insert_pcr):
        # insert gets same extensions on fwd and rev for re1
        insert_prep = self.add_cutsite(self.re1.site, insert_pcr)
        # backbone digested by re2
        backbone_digest = self.digest_backbone(self.re2)
        # TODO add blunting reactions here
        pass
        return [insert_prep, backbone_digest]

    # two enzymes, non-compatible ends
    def double_digest(self, insert_pcr):
        # insert gets same extensions on fwd and rev for re1
        insert_prep = self.add_cutsites(self.re1.site, self.re2.site, insert_pcr)
        # backbone digested by re2
        backbone_digest = self.digest_backbone(self.re2)
        # TODO add blunting reactions here
        pass
        return [insert_prep, backbone_digest]