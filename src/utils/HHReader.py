"""
Parser for hhr result files created with hhblits|hhsearch|hhalign -o <hhr_file>
This script was copied from:
https://github.com/soedinglab/hh-suite/blob/master/scripts/hh_reader.py
and then modified slighty.
"""
import sys
import os
import argparse
from collections import namedtuple

hhr_alignment = namedtuple('hhr_alignment', ['query_id', 'query_length', 'query_neff',
                                             'template_id', 'template_length', 'template_info',
                                             'template_neff', 'query_ali', 'template_ali',
                                             'start', 'end', 'probability', 'evalue', 'score',
                                             'aligned_cols', 'identity', 'similarity', 'sum_probs'])


class HHRFormatError(Exception):
    def __init__(self, value):
        self.value = "ERROR: "+value

    def __str__(self):
        return repr(self.value)


def get_sequence_name(header):
    name = header.replace(">", "").split()[0]
    return name


def parse_result(lines):
    results = []

    query_id = None
    query_length = None
    query_neff = None
    query_seq = []
    template_id = None
    template_length = None
    template_seq = []
    template_info = None
    query_start = None
    query_end = None
    template_start = None
    template_end = None
    probability = None
    evalue = None
    score = None
    identity = None
    similarity = None
    template_neff = None
    sum_probs = None
    aligned_cols = None

    skipped_ali_tags = ["ss_dssp", "ss_pred", "Consensus"]

    is_alignment_section = False

    for line in lines:
        if(line.startswith("Query")):
            query_id = line.split()[1]
        elif(line.startswith("Match_columns")):
            query_length = int(line.split()[1])
        elif(line.startswith("Neff")):
            query_neff = float(line.split()[1])
        elif(is_alignment_section and (line.startswith("No") or line.startswith("Done!"))):
            if query_start is not None:
                result = hhr_alignment(query_id, query_length, query_neff,
                                       template_id, template_length, template_info, template_neff,
                                       "".join(query_seq), "".join(template_seq), (query_start, template_start),
                                       (query_end, template_end), probability, evalue, score,
                                       aligned_cols, identity, similarity, sum_probs)
                results.append(result)
            template_id = None
            template_info = None
            query_seq = []
            template_seq = []

            query_start = None
            query_end = None
            template_start = None
            template_end = None
        elif(line.startswith("Probab")):
            tokens = line.split()
            probability = float(tokens[0].split("=")[1])
            evalue = float(tokens[1].split("=")[1])
            score = float(tokens[2].split("=")[1])
            aligned_cols = int(tokens[3].split("=")[1])
            identity = float(tokens[4].split("=")[1].replace("%", ""))
            similarity = float(tokens[5].split("=")[1])
            sum_probs = float(tokens[6].split("=")[1])
            if(len(tokens) > 7):
                template_neff = float(tokens[7].split("=")[1])
            continue
        elif(line.startswith(">")):
            is_alignment_section = True
            template_id = line[1:].split()[0]
            template_info = line[1:].strip('\n')
        elif(line.startswith("Q")):
            tokens = line.split()
            if(tokens[1] in skipped_ali_tags):
                continue

            try:
                token_2 = tokens[2].replace("(", "").replace(")", "")
                token_2 = int(token_2)
            except:
                raise HHRFormatError(("Converting failure of start index ({}) "
                                      "of query alignment").format(tokens[2]))

            if query_start is None:
                query_start = token_2
            query_start = min(query_start, token_2)

            try:
                token_4 = tokens[4].replace("(", "").replace(")", "")
                token_4 = int(token_4)
            except:
                raise HHRFormatError(("Converting failure of end index ({}) "
                                      "of query alignment").format(tokens[4]))

            if query_end is None:
                query_end = token_4
            query_end = max(query_end, token_4)
            query_seq.append(tokens[3])
        elif(line.startswith("T")):
            tokens = line.split()
            if(tokens[1] in skipped_ali_tags):
                continue
            template_seq.append(tokens[3])

            try:
                token_2 = tokens[2].replace("(", "").replace(")", "")
                token_2 = int(token_2)
            except:
                raise HHRFormatError(("Converting failure of start index ({}) "
                                      "of template alignment").format(tokens[2]))

            if template_start is None:
                template_start = token_2
            template_start = min(template_start, token_2)

            try:
                token_4 = tokens[4].replace("(", "").replace(")", "")
                token_4 = int(token_4)
            except:
                raise HHRFormatError(("Converting failure of end index ({}) "
                                      "of template alignment").format(tokens[4]))

            if template_end is None:
                template_end = token_4
            template_end = max(template_end, token_4)

            try:
                token_5 = tokens[4].replace("(", "").replace(")", "")
                token_5 = int(token_5)
            except:
                raise HHRFormatError(("Converting failure of template length ({}) "
                                      "in template alignment").format(tokens[5]))
            template_length = token_5


    if(template_id is not None and query_start is not None):
        result = hhr_alignment(query_id, query_length, query_neff,
                               template_id, template_length, template_info, template_neff,
                               "".join(query_seq), "".join(template_seq), (query_start, template_start),
                               (query_end, template_end), probability, evalue, score,
                               aligned_cols, identity, similarity, sum_probs)
        results.append(result)

    return results


def read_result(input_file):
    with open(input_file, encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()
        return parse_result(lines)


def main():
    parser = argparse.ArgumentParser(description="HMM output analysis", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-ifile", "--input_hhr_file", required=True, help="input hhr file")
    parser.add_argument("-ofile", "--output_tsv_file", required=True, help="output tsv file")
    parser.add_argument("-maxeval", "--max_evalue", required=True, type=float, help="maximum e-value for output")

    args = parser.parse_args()
    ifile = args.input_hhr_file
    if not os.path.isfile(ifile) or not ifile.endswith('.hhr'):
        print ("Input file " + ifile + " does not exist or is not hhr. Exiting ...")
        sys.exit()

    ofile = args.output_tsv_file
    maxeval = float(args.max_evalue)
    fh = open(ofile, "w")

    fh.write("alignment\tevalue\tprobability\tannotation\tidentity\tcoverage\tQUERYstart\tQUERYend\tTEMPLATEstart\tTEMPLATEend\n")
    counter = 1
    for result in read_result(ifile):
        if result.evalue > maxeval:
            break
        coverage = (result.end[0] - result.start[0] + 1) * 100.0 / result.query_length
        # NOTE(ophira) this fix is due to annotations which contain tabs
        annotation = str(result.template_info).replace("\t", "")
        line = str(counter) + "\t" + str(result.evalue) + "\t" + str(result.probability) + "\t" + \
            annotation + "\t" + str(result.identity) + "\t" + str(coverage) + "\t" + \
            str(result.start[0]) + "\t" + str(result.end[0]) + "\t" + str(result.start[1]) + "\t" + str(result.end[1])
        fh.write(line + "\n")

        counter += 1

    fh.close()

def convert_file(ifile, ofile, maxeval):
    fh = open(ofile, "w")

    fh.write("alignment\tevalue\tprobability\tannotation\tidentity\tcoverage\tQUERYstart\tQUERYend\tTEMPLATEstart\tTEMPLATEend\tScore\n")
    counter = 1
    for result in read_result(ifile):
        if result.evalue > maxeval:
            break
        coverage = (result.end[0] - result.start[0] + 1) * 100.0 / result.query_length
        # NOTE(ophira) this fix is due to annotations which contain tabs
        annotation = str(result.template_info).replace("\t", "")
        line = str(counter) + "\t" + str(result.evalue) + "\t" + str(result.probability) + "\t" + \
            annotation + "\t" + str(result.identity) + "\t" + str(coverage) + "\t" + \
            str(result.start[0]) + "\t" + str(result.end[0]) + "\t" + str(result.start[1]) + "\t" + str(result.end[1]) + "\t" + str(result.score)
        fh.write(line + "\n")

        counter += 1

    fh.close()



if __name__ == "__main__":
    main()
