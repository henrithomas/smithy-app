import re
import csv
from datetime import datetime
from util import read_large_file


# def read_large_file(file_object):
#     while True:
#         data = file_object.readline()
#         if not data:
#             break
#         yield data


def generator_parse(parts_in, parts_out):
    start = datetime.now()
    fieldnames = ['part_id', 'part_name', 'sequence', 'sequence_length']
    part = []
    line_count = 0
    re_expression = fr'<field name=\"({fieldnames[0]}|{fieldnames[1]}|{fieldnames[2]}|{fieldnames[3]})\">(.*?)<\/field>'
    try:
        with open(parts_in, 'r', encoding="utf8") as f_in, open(parts_out, 'w', newline='') as f_out:
            writer = csv.writer(f_out)
            for line in read_large_file(f_in):
                if line_count % 1000000 == 0:
                    print(f'lines processed: {line_count}')
                search = re.findall(re_expression, line)
                line_count += 1
                if search:
                    (field_name, value) = search[0]
                    if field_name == fieldnames[0]:
                        part = [value]
                    else:
                        part.append(value)
                        if field_name == fieldnames[-1]:
                            writer.writerow(part)
    except (IOError, OSError):
        print("error opening or processing file")
    print(f'generator_parse done\ntime {datetime.now() - start}')


def iterator_parse(parts_in):
    start = datetime.now()
    try:
        with open(parts_in, 'r', encoding="utf8") as f_in:
            while True:
                print(next(f_in)[:100], end='')
    except (IOError, OSError):
        print("error opening or processing file")
    except StopIteration:
        print(f'iterator_parse done\ntime {datetime.now() - start}')


def open_igem(parts_in):
    with open(parts_in, "rb") as part_file:
        parsed = part_file.read().decode("utf-8", errors="ignore")
        for row in parsed.split("<row>"):
            print(row)
            name_match = re.search(r"<field name=\"(part_name)\">(.*?)</field>", row)
            if name_match:
                grouptwo = name_match.groups()

def igem_fasta():
    file_in = 'E:\\Thesis\\Databases\\iGEM\\my_igem_parts.csv'
    file_out = 'E:\\Thesis\\Databases\\iGEM\\igem.fasta'
    try:
        with open(file_in, 'r', newline='') as f_in, open(file_out, 'w') as f_out:
            partreader = csv.reader(f_in)
            for line in partreader:
                header = f'>{line[1]}\n'
                f_out.write(header)
                seq = f'{line[2]}\n'
                f_out.write(seq)
    except (IOError, OSError):
        print("error opening or processing file")


if __name__ == "__main__":
    file_in = 'E:\\Thesis\\Databases\\iGEM\\xml_parts.xml'
    file_out = 'E:\\Thesis\\Databases\\iGEM\\my_igem_parts.csv'
    # open_igem(file_in)
    # generator_parse(file_in, file_out)
    # igem_fasta()
