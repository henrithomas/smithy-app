import os
import json
from Bio.Blast.Applications import NcbimakeblastdbCommandline


def create_database(cmd_params, db_path='', project_path=''):
    makeblastdb = NcbimakeblastdbCommandline(**cmd_params)
    os.chdir(db_path)
    stdout, stderr = makeblastdb()
    os.chdir(project_path)
    print(stdout, stderr)

if __name__ == "__main__":
    with open('utils\\db.json') as json_file:
        data = json.load(json_file)
    # TODO add a seqkit rmdup cleanup here
    create_database(data['db_params'], data['db_path'], data['project_path'])