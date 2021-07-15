from utils.file import read_large_file
import pandas as pd
import numpy as np
import csv
from collections import defaultdict

def make_seq_map(file):
    mapping = []
    with open(file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            mapping.append([int(row['insertid']), int(row['sequenceid'])])
    return mapping

def make_dnainsert_map(file):
    mapping = []
    with open(file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            mapping.append([int(row['cloneid']), int(row['insertid'])])
    return mapping

def make_seq_map_v2(file):
    mapping = defaultdict(int)
    with open(file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            mapping[int(row['insertid'])] = int(row['sequenceid'])
    return mapping

def make_dnainsert_map_v2(file):
    mapping = defaultdict(int)
    with open(file, newline='', encoding='utf-8-sig') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            mapping[int(row['cloneid'])] = int(row['insertid'])
    return mapping    

def get_id(entry_id, mapping):
    for row in mapping:
        if row[0] == entry_id:
            return row[1]
    return 0

def get_dnainsert(dnaid, reader):
    for row in reader: 
        if int(row['insertid']) == dnaid:
            return row 
    return None

def get_seqtext(seqid, reader):
    for row in reader: 
        if int(row['sequenceid']) == seqid:
            return row 
    return None
 

if __name__ == '__main__':
    clones = 'C:\\Users\\mt200\\Desktop\\projects\\data\\clones_info.csv'
    clone_dna_map = 'C:\\Users\\mt200\\Desktop\\projects\\data\\cloneid_dnainsertid_map.csv'
    dna_path = 'C:\\Users\\mt200\\Desktop\\projects\\data\\dnainsert.csv'
    insert_sequence_map = 'C:\\Users\\mt200\\Desktop\\projects\\data\\insert_sequence_map.csv'
    seqtext = 'C:\\Users\\mt200\\Desktop\\projects\\data\\seqtext.csv'

    flattened_csv = 'E:\\Thesis\\Databases\\DNASU\\flattened\\dnasu_flat.csv'
    fasta = 'C:\\Users\\mt200\\Desktop\\projects\\data\\dnasu_new.fasta'

    fields = ['cloneid', 'clonename', 'clonetype', 'status', 'verified', 'vermethod', 'domain', 
              'subdomain', 'restriction', 'comments', 'clonemapfilename', 'vectorid', 'vectorname', 
              'specialtreatment', 'source', 'description', 'nextgen', 'insertid', 'insertid', 
              'insertorder', 'sizeinbp', 'species', 'hasdiscrepancy', 'hasmutation', 'format', 
              'source', 'geneid', 'name', 'description', 'targetseqid', 'targetgenbank', 'region', 
              'refseqid', 'hasupdate', 'annotation', 'sequenceid', 'seqorder', 'seqtext']

    dnas = []
    dna_lookup = []
    seqs = []
    seq_lookup = []


    print('making clone to dna map')
    clone_to_dnainsert = make_dnainsert_map_v2(clone_dna_map)
    print('making insert to sequence map')
    insert_to_sequence = make_seq_map_v2(insert_sequence_map)

        
    # with open(dna_path, newline='', encoding='utf-8-sig') as dna_file:
    #     dna_reader = csv.DictReader(dna_file)
    #     print('making dna lookup')
    #     dnas = [dna_row for dna_row in dna_reader]

    with open(seqtext, newline='', encoding='utf-8-sig') as seq_file: 
        seq_reader = csv.DictReader(seq_file)
        print('making seq lookup')
        seqs = [seq_row for seq_row in seq_reader]

    # print(f'dna lookup length: {len(dnas)}')
    print(f'seq lookup length: {len(seqs)}')

    # dna_lookup = [(i, int(dna['insertid'])) for i, dna in enumerate(dnas)]
    seq_lookup = [[i, int(seq['sequenceid'])] for i, seq in enumerate(seqs)]

    print('writing data')
    line_count = 0
    with open(clones, newline='', encoding='utf-8-sig') as clones_file, \
        open(fasta, 'w', newline='') as fa_file:
        clones_reader = csv.DictReader(clones_file)

        for clone in clones_reader:
            if line_count % 10000 == 0:
                print(f'lines processed: {line_count}')
            
            insert_id = clone_to_dnainsert[int(clone['cloneid'])]
            seq_id = insert_to_sequence[insert_id]
            seq_text = None

            for lookup in seq_lookup:   
                if lookup[1] == seq_id:
                    seq_text = seqs[lookup[0]]
                    break

            if seq_text:
                fa_file.write(f'>{clone["vectorname"]}-{clone["clonename"]}\n{seq_text["seqtext"]}\n')
            
            line_count += 1