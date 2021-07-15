import os
import pprint
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime
from io import StringIO
from itertools import combinations

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import (NcbiblastnCommandline,
                                    NcbimakeblastdbCommandline)
from Bio.Blast.Record import HSP, Alignment
from pydna.dseqrecord import Dseqrecord
from fragment.tree import FragmentTree

# TODO either make the query files attribute kwarg or make it a param on the query func
class Blaster:
    def __init__(self, max_seqs, databases, min_frag, max_frag):
        self.max_seqs = max_seqs
        self.databases = databases
        self.solution_file = None
        self.sequence_file = None
        self.blastn_format = 5
        self.fragment_min = min_frag
        self.fragment_max = max_frag
        self.synthetic_min = 0
        self.query_length = 0
        self.query_lengths = []

    def index_func(self, align):
        return align.hsps[0].query_start


    def score_func(self, align):
        return align.hsps[0].score


    def solution_func(self, solution):
        return float(solution[0])


    def get_entry(self, entry, database):
        result = subprocess.run(['blastdbcmd', '-db', database, '-entry', entry], capture_output=True, text=True)
        return result.stdout, result.stderr


    def get_entry_batch(self, parts, database, file):
        with open(file, 'w') as f:
            for part in parts:
                f.write(f'{part}\n')
        entries = f'{os.getcwd()}\\{file}'
        result = subprocess.run(['blastdbcmd', '-db', database, '-entry_batch', entries], capture_output=True, text=True)
        return result.stdout, result.stderr

    def queries(self, f_query, f_out=None):
        query_params = []
        for db in self.databases:
            params = {
                'db': db,
                'query': f_query,
                'max_target_seqs': self.max_seqs,
                'outfmt': 5,
                'gapextend': 6,
                'gapopen': 6,
                'reward': 1,
                'penalty': -5,
            }
            if f_out:
                params.update({'out': f_out})
            query_params.append(params)
        return query_params


    def merge_blastn(self, blast_list):
        flattened = []
        # TODO make this alignment/hit flattening not so messy
        for blast in blast_list:
            db_name = blast.database
            for alignment in blast.alignments:
                for hit in alignment.hsps:
                    temp = Alignment()
                    temp.title = alignment.title
                    temp.hit_id = f'{alignment.hit_id}-{db_name}'
                    temp.hit_def = alignment.hit_def
                    temp.length = alignment.length
                    temp.hsps.append(hit)
                    flattened.append(temp)

        return flattened


    def filter_blastn(self, alignments):
        print('filtering blastn results')
        verify_set = {(-1, -1)}
        blast_hits_filtered = []

        alignments.sort(key=self.score_func, reverse=True)
        # TODO what alignment filtering needs to happen based on common practice??
        # TODO filtering tasks: check for internal cut sites --> mutate out internal cut sites, proper fragment length for primer compatibility
        for alignment in alignments:
            start = alignment.hsps[0].query_start
            end = alignment.hsps[0].query_end
            if (start, end) not in verify_set \
                and alignment.hsps[0].align_length >= self.fragment_min \
                    and alignment.hsps[0].align_length <= self.fragment_max:
                blast_hits_filtered.append(alignment)
                verify_set.add((start, end))
            # TODO add redundancy list here

        blast_hits_filtered.sort(key=self.index_func)

        return blast_hits_filtered



    def query_blastn(self, cmd_params, batch=False):
        print(f'running {cmd_params["db"]} blastn query')
        blastn_cmd = NcbiblastnCommandline(**cmd_params)
        stdout, stderr = blastn_cmd()
        if stdout:
            # TODO what other out formats should we consider, or is always assuming XML for convenience alright?
            if batch:
                stdout = list(NCBIXML.parse(StringIO(stdout)))
            else:
                stdout = NCBIXML.read(StringIO(stdout))
        return stdout, stderr


    def multi_query_blastn(self, cmd_params):
        print(f'running {cmd_params["db"]} blastn query')
        blastn_cmd = NcbiblastnCommandline(**cmd_params)
        stdout, stderr = blastn_cmd()
        if stdout:
            stdout = list(NCBIXML.parse(StringIO(stdout)))
        return stdout, stderr


    def run_blastn(self, params_list):
        print('building blastn query')
        query_results = []
        query_errors = []
        query_stderr = ''

        for param in params_list:
            stdout, stderr = self.query_blastn(param)
            query_results.append(stdout)
            query_errors.append(stderr)

        self.query_length = query_results[0].query_length
        alignment_list = self.merge_blastn(query_results)
        filtered = self.filter_blastn(alignment_list)

        return filtered, query_stderr

