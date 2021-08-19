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
    """
    A class for performing local BLAST queries


    Attributes
    ----------
    max_seqs : int
        A maximum count of BLAST alignment results to use from each database query

    databases : list
        A list of databases to query in the local BLAST search

    blastn_format : int
        The output format from the local BLAST queries. See Biopython documentation. 

    fragment_min : int
        The minimum nucleotide size of alignments to keep from BLAST queries

    fragment_max : int
        The maximum nucleotide size of alignments to keep from BLAST queries

    synthetic_min : int
        The minimum nucleotide size of 

    query_length : int
        The length of the BLAST query sequence


    Returns
    -------
    An instance of Blaster


    Methods
    -------
    def index_func(self, align):
        Returns the query_start value from a BLAST alignment

    def score_func(self, align):
        Returns the score value from a BLAST alignment

    def get_entry(self, entry, database):
        Fetches a single fasta entry in a local BLAST database

    def get_entry_batch(self, parts, database, file):
        Fetches a batch of fasta entries in a local BLAST database

    def queries(self, f_query, f_out=None):
        Builds dicts parameters for each BLAST query

    def merge_blastn(self, blast_list):
        Merges results from each BLAST database query

    def filter_blastn(self, alignments):
        Filters results from the merged BLAST results according to parameters 

    def query_blastn(self, cmd_params, batch=False):
        Executes a local BLASTN search for the insert sequence

    def multi_query_blastn(self, cmd_params):
        Executes a local BLASTN search for multiple fasta sequences

    def run_blastn(self, params_list):
        Performs the full BLASTN procedure using the parameters dicts created in queries()
    """

    def __init__(self, max_seqs, databases, min_frag, max_frag):
        """
        Constructs all necessary attributes for Blaster


        Parameters
        ----------
        max_seqs : int
            A maximum count of BLAST alignment results to use from each database query

        databases : list
            A list of databases to query in the local BLAST search

        blastn_format : int
            The output format from the local BLAST queries. See Biopython documentation. 

        fragment_min : int
            The minimum nucleotide size of alignments to keep from BLAST queries

        fragment_max : int
            The maximum nucleotide size of alignments to keep from BLAST queries

        synthetic_min : int
            The minimum nucleotide size of 

        query_length : int
            The length of the BLAST query sequence
        """

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
        """
        Returns the query_start value from a BLAST alignment


        Parameters
        ----------
        align : Biopython BLAST Alignment object
            A BLAST query alignment result


        Returns
        -------
        The query_start value from align
        """
        return align.hsps[0].query_start


    def score_func(self, align):
        """
        Returns the score value from a BLAST alignment


        Parameters
        ----------
        align : Bio.Alignment BLAST Alignment object
            A BLAST query alignment result


        Returns
        -------
        The BLAST score value from align
        """
        return align.hsps[0].score


    def solution_func(self, solution):
        """
        


        Parameters
        ----------



        Returns
        -------
        """
        return float(solution[0])


    def get_entry(self, entry, database):
        """
        Fetches a single fasta entry in a local BLAST database


        Parameters
        ----------
        entry : str
            The local BLAST database entry's name

        database : str
            The local BLAST database to query 


        Returns
        -------
        The results of the database query in stdout and any error messages in stderr
        """
        result = subprocess.run(['blastdbcmd', '-db', database, '-entry', entry], capture_output=True, text=True)
        return result.stdout, result.stderr


    def get_entry_batch(self, entries, database, file):
        """
        Fetches a batch of fasta entries in a local BLAST database


        Parameters
        ----------
        entries : list
            A list of entry names to query from the local BLAST database

        database : str
            The local BLAST database to query

        file : str
            The file path to use to write entry names into for a batch query


        Returns
        -------
        The results of the database query in stdout and any error messages in stderr
        """
        with open(file, 'w') as f:
            for entry in entries:
                f.write(f'{entry}\n')
        entries = f'{os.getcwd()}\\{file}'
        result = subprocess.run(['blastdbcmd', '-db', database, '-entry_batch', entries], capture_output=True, text=True)
        return result.stdout, result.stderr


    def queries(self, f_query, f_out=None):
        """
        Builds dicts parameters for each BLAST query


        Parameters
        ----------
        f_query : str
            The file path for the fasta sequence of the insert to query over the BLAST database

        f_out : optional, str
            An optional file path for the BLAST query results 
            ***Not intended for use in the main Smithy pipeline***

        Returns
        -------
        A list of dicts for each database query with all necessary BLASTN parameters
        """
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
        """
        Merges results from each BLAST database query


        Parameters
        ----------
        blast_list : list
            A list of lists, each list containing individual BLAST query alignment results


        Returns
        -------
        A flattened, combined list of all alignments from the BLAST queries
        """
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
        """
        Filters results from the merged BLAST results according to parameters


        Parameters
        ----------
        alignments : list
            A list of Biopython alignment objects from all local BLAST queries


        Returns
        -------
        An ordered and filtered list of BLAST alignments according to user parameters
        """
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


    def query_blastn(self, cmd_params):
        """
        Executes a local BLASTN search for the insert sequence


        Parameters
        ----------
        cmd_params : dict
            A dict of BLASTN parameters for the query

        batch : bool
            Declares the query to be using a multi- or single-entry fasta file


        Returns
        -------
        An XML formatted output from the BLASTN query in stdout, any error messages in stderr
        """
        print(f'running {cmd_params["db"]} blastn query')
        blastn_cmd = NcbiblastnCommandline(**cmd_params)
        stdout, stderr = blastn_cmd()
        if stdout:
            stdout = NCBIXML.read(StringIO(stdout))
        return stdout, stderr


    def multi_query_blastn(self, cmd_params):
        """
        Executes a local BLASTN search for multiple fasta sequences


        Parameters
        ----------
        cmd_params : dict
            A dict of BLASTN parameters for the query


        Returns
        -------
        An XML formatted list output from the BLASTN query in stdout, any error messages in stderr
        """
        print(f'running {cmd_params["db"]} blastn query')
        blastn_cmd = NcbiblastnCommandline(**cmd_params)
        stdout, stderr = blastn_cmd()
        if stdout:
            stdout = list(NCBIXML.parse(StringIO(stdout)))
        return stdout, stderr


    def run_blastn(self, params_list):
        """
        Performs the full BLASTN procedure using the parameters dicts created in queries()


        Parameters
        ----------
        params_list : list
            A list of dicts for each BLAST query to run


        Returns
        -------
        An ordered and filtered BLASTN alignment results list in filtered and any error messages 
        in query_stderr   
        """
        print('building blastn query')
        query_results = []
        query_errors = []
        query_stderr = ''

        for params in params_list:
            stdout, stderr = self.query_blastn(params)
            query_results.append(stdout)
            query_errors.append(stderr)

        self.query_length = query_results[0].query_length
        alignment_list = self.merge_blastn(query_results)
        filtered = self.filter_blastn(alignment_list)

        return filtered, query_stderr

    def multi_query_filter_blastn(self, alignments, length):
        blast_hits_filtered = []

        for alignment in alignments:
            if alignment.hsps[0].align_length == length:
                blast_hits_filtered.append(alignment)
        
        return blast_hits_filtered

    def run_multi_blastn(self, params_list, lengths, num_fragments):
        query_results = []
        query_errors = []
        merged_query_results = []
        filtered_query_results = []
        query_stderr = ''

        # run the multi-sequence query
        for params in params_list:
            stdout, stderr = self.multi_query_blastn(params)
            query_results.append(stdout)
            query_errors.append(stderr)

        # for each database alignment set for a sequence, merge these results into a flattend list
        for i in range(num_fragments):
            merge_list = []
            for results in query_results:
                merge_list.append(results[i])
            merged_query_results.append(self.merge_blastn(merge_list))

        # filter out incomplete alignments for each sequence queried
        for i, alignments in enumerate(merged_query_results):
            temp = self.multi_query_filter_blastn(alignments, lengths[i])
            filtered_query_results.append(temp)
            

        return filtered_query_results, query_stderr 
