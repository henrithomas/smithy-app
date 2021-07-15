# TODO find a good min_synth fragment value
class FragmentNode:
    def __init__(self, data=None, start=0, end=0, score=0, i='NONE', db='NONE'):
        self.data = data
        self.node_id = i
        self.db = db
        self.start = start
        self.end = end
        self.score = score
        self.level = 0
        self.cost = 0 

    @property
    def info(self):
        return f'{self.node_id}-{self.start}-{self.end}'

    @property
    def query_seq(self):
        return self.data.hsps[0].query

    @property
    def subject_seq(self):
        return self.data.hsps[0].sbjct

    @property
    def coordinates(self):
        return self.start, self.end

    def __add__(self, other):
        if isinstance(other, FragmentNode):
            node_data = self.data
            node_data.hsps.append(other.data.hsps[0])
            return FragmentNode(data=node_data,
                                start=self.start,
                                end=other.end,
                                score=self.score + other.score)
        return NotImplemented
    
    def __str__(self):
        pass

    def __repr__(self):
        pass