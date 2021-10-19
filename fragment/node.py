# TODO find a good min_synth fragment value
class FragmentNode:
    """
    The class for FragmentNodes that is used in the FragmentTree class


    Attributes
    ----------
    data : Biopython BLAST alignment
        A Biopython BLAST Alignment object for a single query result

    start : optional, int
        The start index of the alignment on the query sequence

    end : optional, int
        The end index of the alignment on the query sequence

    score : optional, int 
        The BLAST alignment score for the alignment saved in data

    i : optional, str
        An identifying string for the FragmentNode        

    db : optional, str
        The name of the database for the FragmentNode


    Returns
    -------
    An instance of FragmentNode


    Methods
    -------
    def info(self): 
        @property for a preformatted label for the FragmentNode

    def query_seq(self):
        @property for a pre-fetched query sequence for the FragmentNode

    def subject_seq(self):
        @property for a pre-fetched subject sequence for the FragmentNode

    def coordinates(self):
        @property for a tuple of the FragmentNode's start and end values

    def __add__(self, other):
        Allows two FragmentNodes to be added using the addition operator
    """ 
    def __init__(self, data=None, start=0, end=0, score=0, node_id='NONE', db='NONE', synth=True):
        """
        Constructs all necessary attributes for a FragmentNode


        Parameters
        ----------
        data : Biopython BLAST alignment
            A Biopython BLAST Alignment object for a single query result

        start : optional, int
            The start index of the alignment on the query sequence

        end : optional, int
            The end index of the alignment on the query sequence

        score : optional, int 
            The BLAST alignment score for the alignment saved in data

        i : optional, str
            An identifying string for the FragmentNode        

        db : optional, str
            The name of the database for the FragmentNode
        """
        self.data = data
        self.node_id = node_id
        self.db = db
        self.start = start
        self.end = end
        self.score = score
        self.level = 0
        self.cost = 0 
        self.synth = synth

    @property
    def info(self):
        """
        @property for a preformatted label for the FragmentNode


        Returns
        -------
        A formatted string label for the FragmentNode 
        """
        return f'{self.node_id}: ({self.start}, {self.end})'

    @property
    def query_seq(self):
        """
        @property for a pre-fetched query sequence for the FragmentNode

        
        Returns
        -------
        The query sequence from the BLAST alignement data of the FragmentNode
        """
        return self.data.hsps[0].query

    @property
    def subject_seq(self):
        """
        @property for a pre-fetched subject sequence for the FragmentNode


        Returns
        -------
        The subject sequence from the BLAST alignement data of the FragmentNode
        """
        return self.data.hsps[0].sbjct

    @property
    def coordinates(self):
        """
        @property for a tuple of the FragmentNode's start and end values


        Returns
        -------
        A tuple of the start and end indexes for the FragmentNode
        """
        return self.start, self.end

    def __add__(self, other):
        """
        Allows two FragmentNodes to be added using the addition operator

        
        Parameters
        ----------
        other : FragmentNode
            A FragmentNode instance


        Returns
        -------
        The summation of two FragmentNodes
        """
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