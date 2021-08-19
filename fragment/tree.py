from collections import defaultdict
from datetime import datetime
from itertools import combinations
from fragment.node import FragmentNode
from Bio.Blast.Record import HSP, Alignment
import random

class FragmentTree:
    """
    A class for building cloning assembly solutions for insert queries using BLAST query results.


    Attributes
    ----------
    node_list : optional, list
        A list of FragmentNodes for building a solution tree

    min_synth : optional, int
        Minimum nucleotide size of synthetic fragments to design for assemblies

    max_synth : optional, int
        Maximum nucleotide size of synthetic fragments to design for assemblies

    adjacency_set : defaultdict(set)
        A set of adjacencies between FragmentNodes for building insert solutions

    node_end : set
        A set of FragmentNode ends used for 

    visited : list
        A list of visited nodes to use when building solution

    solutions : list
        A list of lists, where each list is a collection of indexes of FragmentNodes in the nodes list that make
        an ordered solution for the assembly.

    scores : list
        A list of scores for each respective solution. Indexes of the scores list correspond to those of the solutions list. 

    query_len : optional, int
        The length of the insert query sequence

    query : str
        The sequence of the insert sequence

    max : float
        The maximum total score of all solutions  

    Returns
    -------
    An instance of FragmentTree


    Methods
    -------
    def build_edge_set(self):
        Constructs the fragment solution tree/network based on FragmentNode adjacencies

    def root_node(self):   
        Constructs and appends the "root"/starter node for the solution tree to the nodes list

    def filler_node(self, start, end):
        Constructs and appends a filler/synthetic node for the solution tree to the nodes list

    def blast_input(self, fragments):
        Takes a BLAST alignment list and builds the node_list and node_ends collections

    def build(self, input):
        Performs the full FragmentTree building and constructs the node_list, node_ends, adjacency_set, 
        solutions, and finds the max score of all solutions

    def dfs(self, v, path=None):
        Performs a depth-first-search on the FragmentTree network and finds solutions

    def print_solution(self, s):
        Prints a solution, s, in a readable format

    def formatted_solution(self, s):
        Creates a writable line of a solution, s, in comma-separated format

    def solution_seqs(self, s, text=False):
        Creates a list of all sequences for a solution to an assembly insert

    def solution_nodes(self, s):
        Creates a list of all FragmentNodes for a solution to an assembly insert

    def solution_score(self, path):
        Calculates the score of a given solution

    def complete_solutions(self):
        Checks the end of each solution and adds a necessary filler FragmentNode if necessary for 
        a complete insert solution 

    def max_score(self, v, path=None):
        Performs a depth-first search on the FragmentTree to find the max solution score
    """
    def __init__(self, query, nodes=None, query_len=0, min_synth=0, max_synth=1000):
        """
        Constructs all attributes for a FragmentTree


        Parameters
        ----------
        query : str
            The sequence of the insert sequence

        nodes : optional, list
            A list of FragmentNodes for building a solution tree

        query_len : optional, int
            The length of the insert query sequence

        min_synth : optional, int
            Minimum nucleotide size of synthetic fragments to design for assemblies

        max_synth : optional, int
            Maximum nucleotide size of synthetic fragments to design for assemblies
        """
        self.node_list = nodes
        self.min_synth = min_synth
        self.max_synth = max_synth
        self.adjacency_set = defaultdict(set) 
        self.node_ends = set()
        self.visited = None
        self.solutions = []
        self.scores = []
        self.query_len = query_len
        self.query = query
        self.max = 0
        self.multi_query_node_list = []
        self.multi_query_solutions = []

    def build_edge_set(self):
        """
        Constructs the fragment solution tree/network based on FragmentNode adjacencies


        Parameters
        ----------
        None


        Returns
        -------
        None
        """
        for i, j in combinations(range(len(self.node_list)), 2):
            node_a_start, node_a_end = self.node_list[i].coordinates
            node_b_start, node_b_end = self.node_list[j].coordinates
            
            if node_b_start > (node_a_end + 1 + self.min_synth) \
                and node_b_start < (node_a_end + 1 + self.max_synth) \
                    and (node_b_start - 1) not in self.node_ends:
                self.filler_node(node_a_end + 1, node_b_start - 1)
                self.node_ends.add(node_b_start - 1)
            if node_a_start > (node_b_end + 1 + self.min_synth) \
                and node_a_start < (node_b_end + 1 + self.max_synth) \
                    and (node_a_start - 1) not in self.node_ends:
                self.filler_node(node_b_end + 1, node_a_start - 1)
                self.node_ends.add(node_a_start - 1)

        for i, j in combinations(range(len(self.node_list)), 2):
            node_a_start, node_a_end = self.node_list[i].coordinates
            node_b_start, node_b_end = self.node_list[j].coordinates

            if (node_b_start == node_a_end) or (node_b_start == (node_a_end + 1)) and j not in self.adjacency_set[i]:
                self.adjacency_set[i].add(j)
            if (node_a_start == node_b_end) or (node_a_start == (node_b_end + 1)) and i not in self.adjacency_set[j]:
                self.adjacency_set[j].add(i)

    def root_node(self):
        """
        Constructs and appends the "root"/starter node for the solution tree to the nodes list


        Parameters
        ----------
        None


        Returns
        -------
        None
        """
        root_data = Alignment()
        root_data.title = 'ROOT'
        root_data.hit_id = 'ROOT'
        root_data.hit_def = 'ROOT'
        temp_hsp = HSP()
        temp_hsp.score = 0
        temp_hsp.query = '_____'
        temp_hsp.sbjct = '_____'
        temp_hsp.query_start = 0
        temp_hsp.query_end = 0
        temp_hsp.sbjct_start = 0
        temp_hsp.sbjct_end = 0
        root_data.hsps.append(temp_hsp)
        self.node_list.append(FragmentNode(data=root_data, start=0, end=0, score=0, i='ROOT'))

    def filler_node(self, start, end): 
        """
        Constructs and appends a filler/synthetic node for the solution tree to the nodes list


        Parameters
        ----------
        None


        Returns
        -------
        None
        """
        temp = Alignment()
        temp.title = f'FILLER-{start}-{end}'
        temp.hit_id = f'FILLER-{start}-{end}'
        temp.hit_def = f'FILLER-{start}-{end}'
        temp.length = end - start
        temp_hsp = HSP()
        temp_hsp.score = 0
        temp_hsp.query = self.query[start - 1:end]
        temp_hsp.sbjct = self.query[start - 1:end]
        temp_hsp.query_start = start
        temp_hsp.query_end = end
        temp_hsp.sbjct_start = 0
        temp_hsp.sbjct_end = 0
        temp.hsps.append(temp_hsp)
        self.node_list.append(FragmentNode(data=temp, start=start, end=end, score=0, i=f'FILLER-{start}-{end}'))

    def blast_input(self, fragments):
        """
        Takes a BLAST alignment list and builds the node_list and node_ends collections


        Parameters
        ----------
        fragments : list
            A list of Biopython BLAST alignment objects from a BLAST query 


        Returns
        -------
        None
        """
        self.root_node()
        for fragment in fragments:
            self.node_list.append(FragmentNode(data=fragment,
                                               start=fragment.hsps[0].query_start,
                                               end=fragment.hsps[0].query_end,
                                               score=fragment.hsps[0].score,
                                               i=fragment.hit_id,
                                               db=fragment.hit_id.split('-')[1]))
            self.node_ends.add(fragment.hsps[0].query_end)

    def build(self, input):
        """
        Performs the full FragmentTree building and constructs the node_list, node_ends, adjacency_set, 
        solutions, and finds the max score of all solutions


        Parameters
        ----------
        input : list
            A list of Biopython BLAST alignment objects from a BLAST query 


        Returns
        -------
        None
        """
        self.blast_input(input)
        print(f'number of fragments: {len(self.node_list)}')

        print('building fragment tree')
        startime = datetime.now()
        self.build_edge_set()
        print(f'fragment network build time: {datetime.now() - startime}')
        
        self.visited = [False] * len(self.node_list)
        self.dfs(0, [])
        print(f'number of solutions: {len(self.solutions)}')

        self.max = max(self.scores)
        print(f'highest score: {self.max}')

        self.solutions.sort(reverse=True)
        self.complete_solutions()

    def dfs(self, v, path=None):
        """
        Performs a depth-first-search on the FragmentTree network and finds solutions


        Parameters
        ----------
        v : int
            The node list index

        path : optional, list
            A list of node list indexes for a path to a solution in the FragmentTree


        Returns
        -------
        None
        """
        if not path:
            path = []
        self.visited[v] = True
        path.append(v)

        if not self.adjacency_set[v]:
            score = self.solution_score(path)
            self.scores.append(score)

            p = path.copy()
            p.insert(0, score)
            self.solutions.append(p)
        else:
            for i in self.adjacency_set[v]:
                if not self.visited[i]:
                    self.dfs(i, path)
        path.pop()
        self.visited[v] = False

    def print_solution(self, s):
        """
        Prints a solution, s, in a readable format


        Parameters
        ----------
        s : int
            The index for a solution in the solutions list

        Returns
        -------
        A readable string for a solution
        """
        if self.solutions:
            for i in self.solutions[s]:
                fragment = self.node_list[i]
                print(f'{fragment.node_id}-{fragment.start}-{fragment.end}, seq={fragment.data.hsps[0].sbjct}')
            test = self.node_list[self.solutions[s][-1]].end
            if test != self.query_len:
                print(f'FILLER-{test}-{self.query_len}, seq=_____')

    def formatted_solution(self, s):
        """
        Creates a writable line of a solution, s, in comma-separated format


        Parameters
        ----------
        s : int
            The index for a solution in the solutions list

        Returns
        -------
        A writable string for a solution
        """
        line = ''
        if self.solutions:
            nodes = []
            for i in self.solutions[s][1:]:
                nodes.append(self.node_list[i])
            line = ', '.join([f'{node.node_id}-{node.start}-{node.end}' for node in nodes])
            test = self.node_list[self.solutions[s][-1]].end
            if test != self.query_len:
                line += f', FILLER-{test}-{self.query_len}'
        return line

    def solution_seqs(self, s, text=False):
        """
        Creates a list of all sequences for a solution to an assembly insert


        Parameters
        ----------
        s : int
            The index for a solution in the solutions list

        text : optional, str
            An option to return the list as a comma-separated string for file writing or a list of strings


        Returns
        -------
        Either a list of solution sequences or a string in comma-separated format of the sequences
        """
        if self.solutions:
            seqs = [node.subject_seq for node in [self.node_list[i] for i in self.solutions[s][2:]]]
            if text:
                return ', '.join(seqs)
            else:
                return seqs

    def solution_nodes(self, s):
        """
        Creates a list of all FragmentNodes for a solution to an assembly insert



        Parameters
        ----------
        s : int
            The index for a solution in the solutions list


        Returns
        -------
        A list of FragmentNodes for a solution
        """
        nodes = [node for node in [self.node_list[i] for i in self.solutions[s][2:]]]
        return nodes
    
    def solution_names(self, s):
        pass

    def solution_score(self, path):
        """
        Calculates the score of a given solution


        Parameters
        ----------
        path : int
            A list of node list indexes for a solution


        Returns
        -------
        The calculated total of a solution's score
        """
        nodes = []
        for i in path:
            nodes.append(self.node_list[i])
        # nodes = [self.node_list[i] for i in path]
        # return sum([node.score for node in [self.node_list[i] for i in path]])
        return sum([node.score for node in nodes])

    def solution_score_v2(self, path):
        return sum([node.score for node in [self.node_list[i] for i in path]])
    
    def complete_solutions(self):
        """
        Checks the end of each solution and adds a necessary filler FragmentNode if necessary for 
        a complete insert solution 


        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        del_list = []
        
        for i, solution in enumerate(self.solutions):
            last = solution[-1]
            test = self.node_list[last].end
            remainder = self.query_len - test
            if remainder <= self.max_synth and remainder >= self.min_synth:
                self.filler_node(test, self.query_len)
                self.adjacency_set[last].add(len(self.node_list) - 1)
                solution.append(len(self.node_list) - 1)
            else: 
                del_list.append(i)

        self.solutions = [solution for i, solution in enumerate(self.solutions) if i not in del_list]        
        pass

    def max_score(self, v, path=None):

        """
        Performs a depth-first search on the FragmentTree to find the max solution score


        Parameters
        ----------
        v : int
            The node list index

        path : optional, list
            A list of node list indexes for a path to a solution in the FragmentTree


        Returns
        -------
        None
        """
        if not path:
            path = []

        self.visited[v] = True
        path.append(v)

        if not self.adjacency_set[v]:
            score = self.solution_score(path)
            if score > self.max:
                self.max = score
        else:
            for i in self.adjacency_set[v]:
                if not self.visited[i]:
                    self.max_score(i, path)

        path.pop()
        self.visited[v] = False

    def multi_query_blast_input(self, fragment_lists, sequences):
        for i, frag_list in enumerate(fragment_lists):
            temp = []
            if not frag_list:
                # If there were no results in any of the databases make a filler node
                seq_len = len(sequences[i])

                temp_align = Alignment()
                temp_align.title = f'FILLER-0-{seq_len}'
                temp_align.hit_id = f'FILLER-0-{seq_len}'
                temp_align.hit_def = f'FILLER-0-{seq_len}'
                temp_align.length = seq_len
                temp_hsp = HSP()
                temp_hsp.score = 0
                temp_hsp.query = sequences[i]
                temp_hsp.sbjct = sequences[i]
                temp_hsp.query_start = 0
                temp_hsp.query_end = seq_len
                temp_hsp.sbjct_start = 0
                temp_hsp.sbjct_end = 0
                temp_align.hsps.append(temp_hsp)
                
                temp.append(temp_align)
            else:
                for fragment in frag_list:
                    temp.append(FragmentNode(data=fragment,
                                            start=fragment.hsps[0].query_start,
                                            end=fragment.hsps[0].query_end,
                                            score=fragment.hsps[0].score,
                                            i=fragment.hit_id,
                                            db=fragment.hit_id.split('-')[1]))
            self.multi_query_node_list.append(temp)

    def solution_count(self, s):
        minimum  = min([len(nodes) for nodes in self.multi_query_node_list])
        if minimum < s:
            s = minimum
        return s 

    def build_multi_query_solutions(self, maximum):
        for i in range(maximum):
            solution = []
            for nodes in self.multi_query_node_list:
                solution.append(random.randint(0, len(nodes) - 1))
            self.multi_query_solutions.append(solution)

    def multi_query_solution_seqs(self, s):
        solution_indexes = self.multi_query_solutions[s]
        seqs = []
        for i, node_list in enumerate(self.multi_query_node_list):
            seqs.append(node_list[solution_indexes[i]].subject_seq)
        return seqs
    
    def multi_query_solution_nodes(self, s):
        solution_indexes = self.multi_query_solutions[s]
        nodes = []
        for i, node_list in enumerate(self.multi_query_node_list):
            nodes.append(node_list[solution_indexes[i]])
        return nodes
