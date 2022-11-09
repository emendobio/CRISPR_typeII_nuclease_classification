import os
import subprocess
import argparse
import math
import re
import json
import ete3
import tempfile
import shutil
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import sys
import os
src_code_path = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1])
sys.path.append(src_code_path)
print(src_code_path)
from utils import FileUtils as eu

def main():

    learn_phylogenetic_tree(args.wd, args.i, args.c)

def learn_phylogenetic_tree(work_directory, input_fasta, clades_file):

    eu.BasicUtils.make_path(work_directory)

    ## first pass to select subtypes
    subtypes = ()
    with open(input_fasta, "r") as ofh:
        subtypes = {id2subtype(record.id) for record in SeqIO.parse(input_fasta, "fasta")}

    ## Get sequence names for all sequences in all clades

    clade_ids, seq2clade = get_clade_ids(clades_file)

    all_sequences = get_fasta_sequences(input_fasta)

    create_alignments(work_directory, clade_ids, input_fasta)

    create_hhsuite_profiles(work_directory, clade_ids)

    create_hhm_db(work_directory, clade_ids)

    clades_scores = compare_clades(work_directory, clade_ids)

    if None in subtypes:
        subtypes.remove(None)

    # for each clade alignment scores of clademembers and outmembers will be collected. Output json with all scores

    get_scores(work_directory, clade_ids, all_sequences)


def create_ffindex_from_msa_files(aligned_clades_files, wd):
    with tempfile.TemporaryDirectory() as temp_dir:
        for muscle_file in aligned_clades_files:
            temp_file = "{}/{}".format(temp_dir, muscle_file.split("/")[-1])
            shutil.copyfile(muscle_file, temp_file)
        runCommand(f'ffindex_build -s {wd}/hhdb_msa.ffdata {wd}/hhdb_msa.ffindex {temp_dir}')


def create_hhm_db(wd, clade_ids):

    aligned_clades_files = [wd + '/' + c + '_muscle.fa' for c in clade_ids]

    # create msa db
    create_ffindex_from_msa_files(aligned_clades_files, wd)

    # create a3m db

    runCommand(f'ffindex_apply {wd}/hhdb_msa.ffdata {wd}/hhdb_msa.ffindex \
                    -i {wd}/hhdb_a3m.ffindex -d {wd}/hhdb_a3m.ffdata -- hhconsensus \
                    -M 50 -maxres 65535 -i stdin -oa3m stdout -v 0')

    # create hhm db

    runCommand(f'ffindex_apply {wd}/hhdb_a3m.ffdata {wd}/hhdb_a3m.ffindex \
                 -i {wd}/hhdb_hhm.ffindex -d {wd}/hhdb_hhm.ffdata -- hhmake \
                 -i stdin -o stdout -v 0')

    # Computing context states for pre-filtering

    runCommand(f'cstranslate -f -x 0.3 -c 4 -I a3m -i {wd}/hhdb_a3m -o {wd}/hhdb_cs219')


def id2subtype(id):
    subtype_re = re.search("\/([cC][aA][sS][\w\-]+)$" , id)
    if subtype_re is not None:
        subtype = subtype_re.groups()[0]
        return subtype
    else:
        return None

def remove_sequence(align, id, check_empty_columns=True):
    align_neto = MultipleSeqAlignment([])
    for seq in align:
        if seq.id != id:
            align_neto.append(seq)
    if check_empty_columns:
        align_neto = remove_empty_columns(align_neto)

    return(align_neto)

def remove_empty_columns(align):
    align_cut = MultipleSeqAlignment(align)
    good_columns = []
    for i in range(0, align.get_alignment_length()):
        if set(align[:,i]) != {'-'}:
            good_columns.append(i)

    for c in align_cut:
        orig_seq = c.seq
        c.seq = ''
        for i in good_columns:
            c.seq += orig_seq[i]
    return(align_cut)

def get_clade_ids(clades_list):
    clade_seqs = {}
    seq2clade = {}
    with open(clades_list) as ifh:
        for treefile in ifh:
            #print(treefile)
            treefile = treefile.rstrip()
            cladename = os.path.splitext(os.path.basename(treefile))[0]
            # get the tree leaves and insert to the dict
            t = ete3.Tree(treefile, format=1)
            clade_seqs[cladename] = {n.name for n in t.get_leaves()}
            seq2clade = {n.name: cladename for n in t.get_leaves()}

    return clade_seqs, seq2clade

def get_fasta_sequences(fasta_file):
    sequences_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        #print(record.id, record.seq)
        sequences_dict[record.id] = record.seq
    return sequences_dict

def create_alignments(wd, clade_ids, fastafile):
    # create alignment for each clade separately using muscle
    for c in clade_ids:
        retained_sequences = []
        for record in SeqIO.parse(fastafile, "fasta"):
            if record.id in clade_ids[c]:
                retained_sequences.append(record)
        retained_sequence_iterator = iter(retained_sequences)
        SeqIO.write(retained_sequence_iterator, f'{wd}/{c}.fa', "fasta")
        muscle_in = f'{wd}/{c}.fa'
        muscle_out = f'{wd}/{c}_muscle.fa'
        assert os.path.exists(muscle_in), muscle_in
        runCommand(f'muscle -in {muscle_in} -out  {muscle_out}')
        assert os.path.exists(muscle_out), muscle_out

def create_hhsuite_profiles(wd, clade_ids):
    for c in clade_ids:
        # generate generic profile (a3m file) with hhmake
        runCommand(f'hhmake -i {wd}/{c}_muscle.fa -o {wd}/{c}_muscle.a3m -M first')

def get_scores(wd, clade_ids, all_sequences):
    # take each sequence - remove it from the clades alignment - create a3m without it
    # make alignment to this a3m and to all other clades and keep the scores
    train_scores = {}

    for c in clade_ids:
        train_scores[c] = {'internal':[], 'external':[]}
        si = 0
        for s in clade_ids[c]:
            self_clade_results = extract_hhalign_score(f'{wd}/{c}.{c}.hhalign')
            si += 1
            # create fasta file with the sequence
            ofh = open(f'{wd}/seqtmp.fa', 'w')
            print(f'>{s}\n{all_sequences[s]}\n', file=ofh)
            ofh.close()
            clade_align = AlignIO.read(f'{wd}/{c}_muscle.fa', "fasta")
            align_neto = remove_sequence(clade_align, s)
            write_fasta(align_neto, f'{wd}/cladetmp.fa')

            # calculate seq self score here

            runCommand(f'hhalign -i {wd}/seqtmp.fa -t {wd}/seqtmp.fa', stdout=f'{wd}/selftmp.hhalign')
            self_results = extract_hhalign_score(f'{wd}/selftmp.hhalign')


            # generate generic profile (a3m file) with hhmake
            runCommand(f'hhmake -i {wd}/cladetmp.fa -o {wd}/cladetmp.a3m -M first')
            runCommand(f'hhalign -i {wd}/seqtmp.fa -t {wd}/cladetmp.a3m', stdout=f'{wd}/inttmp.hhalign')
            results = extract_hhalign_score(f'{wd}/inttmp.hhalign')
            results['SelfScore'] = self_results['Score']
            results['SelfCladeScore'] = self_clade_results['Score']
            results['d'] = -math.log(float(results['Score'])/min(float(self_results['Score']), float(self_clade_results['Score'])))
            train_scores[c]['internal'].append(results)

            for c2 in clade_ids:
                if c == c2:
                    continue
                self_clade_results = extract_hhalign_score(f'{wd}/{c2}.{c2}.hhalign')
                runCommand(f'hhalign -i {wd}/seqtmp.fa -t {wd}/{c2}_muscle.a3m', stdout=f'{wd}/exttmp.hhalign')
                results = extract_hhalign_score(f'{wd}/exttmp.hhalign')
                results['SelfScore'] = self_results['Score']
                results['SelfCladeScore'] = self_clade_results['Score']
                if 'Score' not in results:
                    print(f'hhalign -i {wd}/seqtmp.fa -t {wd}/{c2}_muscle.a3m')
                    print("Score zero: " + ' ' + c + ' ' + s + ' ' + c2)
                    print("Score zero: " + ' ' + c + ' ' + s + ' ' + c2, file=ofh)
                    results['Score'] = 0.0001
                    exit(0)
                results['d'] = -math.log(float(results['Score'])/min(float(self_results['Score']), float(self_clade_results['Score'])))
                train_scores[c]['external'].append(results)

    with open(f'{wd}/data.json', 'w') as jfh:
        json.dump(train_scores, jfh)


def runCommand(command, stdout=None, stderr=None):
    # process = subprocess.Popen(command, shell=True)
    if stdout is not None:
        ofh = open(stdout, 'w')
    else:
        ofh = None
    if stderr is not None:
        efh = open(stderr, 'w')
    else:
        efh = None

    process = subprocess.Popen(command, shell=True, stdout=ofh, stderr=efh)
    process.wait()

def compare_clades(wd, clade_ids, keep=True):
    # use hhalign to make all-vs-all clades comparison
    comparison_table = pd.DataFrame(index=list(clade_ids), columns=list(clade_ids))
    inter_clade_results = {}
    for c1 in clade_ids:
        inter_clade_results[c1] = {}
        for c2 in clade_ids:
            print(f'hhalign -i {wd}/{c1}_muscle.a3m -t {wd}/{c2}_muscle.a3m')

            runCommand(f'hhalign -i {wd}/{c1}_muscle.a3m -t {wd}/{c2}_muscle.a3m', stdout=f'{wd}/{c1}.{c2}.hhalign')

            inter_clade_results[c1][c2] = extract_hhalign_score(f'{wd}/{c1}.{c2}.hhalign')
            if not keep:
                os.remove(f'{wd}/{c1}.{c2}.hhalign')
    return inter_clade_results

def extract_hhalign_score(results_file):
    results = {}
    hits_flag = False
    with open(results_file) as ifh:
        for line in ifh:
            line = line.rstrip()
            if line.startswith('Query'):
                line_splitted = line.split()
                results['query'] = line_splitted[1]
                if len(line_splitted) > 1:
                    results['query_desc'] = ' '.join(line_splitted[2:])
                else:
                    results['query_desc'] = ''
                print(results['query'])
            if hits_flag and line != '':
                line_splitted = line[0:34].split()
                results['hit'] = line_splitted[1]
                if len(line_splitted) > 1:
                    results['hit_desc'] = ' '.join(line_splitted[2:])
                else:
                    results['hit_desc'] = ''
                results['Prob'] = line[35:40]
                results['E-value'] = line[40:48]
                results['P-value'] = line[48:56]
                results['Score'] = line[56:63]
                results['SS'] = line[63:69]
                results['Cols'] = line[69:74]
            if line.startswith(' No Hit'):
                hits_flag = True
    return results

def write_fasta(msa_object, file):
    with open(file, 'w') as ofh:
        for record in msa_object:
            print('>' + record.description + '\n' + record.seq, file=ofh)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create HHM per nuclease subtype. Detect inter and intra score ranges")
    parser.add_argument("-wd",help="working directory", default='./')
    parser.add_argument("-i",help="input multi fasta file")
    parser.add_argument("-c",help="list of tree files which defines the clades")
    parser.add_argument("-j",help="json file to store alignment results", default='data.json')
    args = parser.parse_args()

    main()

