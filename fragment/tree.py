from collections import defaultdict
from datetime import datetime
from itertools import combinations
from re import split
from fragment.node import FragmentNode
from Bio.Blast.Record import HSP, Alignment

class FragmentTree:
    def __init__(self, query, nodes=None, query_len=0, min_synth=0, max_synth=0):
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

    def build_edge_set(self):
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
        temp = Alignment()
        temp.title = 'FILLER'
        temp.hit_id = 'FILLER'
        temp.hit_def = 'FILLER'
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
        self.node_list.append(FragmentNode(data=temp, start=start, end=end, score=0, i=f'FILLER'))

    def blast_input(self, fragments):
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
        if self.solutions:
            for i in self.solutions[s]:
                fragment = self.node_list[i]
                print(f'{fragment.node_id}-{fragment.start}-{fragment.end}, seq={fragment.data.hsps[0].sbjct}')
            test = self.node_list[self.solutions[s][-1]].end
            if test != self.query_len:
                print(f'FILLER-{test}-{self.query_len}, seq=_____')

    def formatted_solution(self, s):
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
        if self.solutions:
            seqs = [node.subject_seq for node in [self.node_list[i] for i in self.solutions[s][2:]]]
            if text:
                return ', '.join(seqs)
            else:
                return seqs

    def solution_nodes(self, s):
        nodes = [node for node in [self.node_list[i] for i in self.solutions[s][2:]]]
        return nodes
    
    def solution_names(self, s):
        pass

    def solution_score(self, path):
        nodes = []
        for i in path:
            nodes.append(self.node_list[i])
        # nodes = [self.node_list[i] for i in path]
        # return sum([node.score for node in [self.node_list[i] for i in path]])
        return sum([node.score for node in nodes])

    def solution_score_v2(self, path):
        return sum([node.score for node in [self.node_list[i] for i in path]])
    
    def complete_solutions(self):
        for solution in self.solutions:
            last = solution[-1]
            test = self.node_list[last].end
            if test < self.query_len:
                self.filler_node(test, self.query_len)
                self.adjacency_set[last].add(len(self.node_list) - 1)
                solution.append(len(self.node_list) - 1)

    
    def max_score(self, v, path=None):
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