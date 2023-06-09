import os
import glob
import tempfile
import sys
import numpy as np
import pandas as pd
src_code_path = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1])
sys.path.append(src_code_path)
from utils import HHReader as hh

def get_database_name(db):
    suffix_list = ["_a3m.ffdata", "_cs219.ffdata", "_a3m.ffindex", "_cs219.ffindex", ".ffindex"]
    assert os.path.exists(db), db
    files_list = glob.glob("{}/*{}".format(db, suffix_list[0]))
    ffindex_file = [x for x in files_list]
    assert len(ffindex_file) == 1, ffindex_file
    db_name = ffindex_file[0][:-len(suffix_list[0])]
    for x in ["_a3m.ffdata", "_cs219.ffdata", "_a3m.ffindex", "_cs219.ffindex"]:
        curr_file = "{}{}".format(db_name, x)
        assert os.path.exists(curr_file), curr_file
    return db_name

def run_hhsearch(input_fasta, db_folder, out_file,cpu_num=2):
    assert os.path.exists(input_fasta)
    db = get_database_name(db_folder)
    os.system("hhblits -cpu {} -i {} -d {} > {}".format(cpu_num, input_fasta, db, out_file))
    assert os.path.exists(out_file)
    assert os.path.getsize(out_file) > 0, "{} is empty".format(out_file)


def run_hhsearch_locally(protein_sequence, db, fasta_file, hhr_file,cpu_num=2):
    f = open(fasta_file, "w")
    f.write(">prot\n{}\n".format(protein_sequence))
    f.close()
    run_hhsearch(fasta_file, db, hhr_file,cpu_num)


def run_single_sequence(protein_sequence, db, cpu_num=2):
    with tempfile.TemporaryDirectory() as working_local_dir:
        single_fasta = '{}/single_seq.fasta'.format(working_local_dir)
        hhr_file = '{}/single_seq.hhr'.format(working_local_dir)
        tsv_file = '{}/single_seq.tsv'.format(working_local_dir)
        run_hhsearch_locally(protein_sequence, db, single_fasta, hhr_file, cpu_num)
        hh.convert_file(hhr_file, tsv_file, 100)
        df = pd.read_csv(tsv_file, sep='\t').sort_values('Score', ascending=False)
    return df


def get_hmm_score(hmmer_res, hmm_name):
    min_prob = 0
    I = (hmmer_res["annotation"] == hmm_name) & (hmmer_res["probability"] >= min_prob)
    if I.sum() == 0 :
        return 0
    else:
        return hmmer_res[I]['Score'].max()

def get_all_hmm_scores(folder, hmm_list):
    n = len(hmm_list)
    df = pd.read_csv("{}/hmm_out.tsv".format(folder), sep="\t")
    I = ~df['annotation'].isin(hmm_list)
    assert I.sum() == 0 , list(df[I]['annotation'])
    scores_vec = np.zeros(n)
    for i in range(n):
        scores_vec[i] = get_hmm_score(df, hmm_list[i])
    return scores_vec


def hhalign_command(file1, file2, out):
    hhalign_command = 'hhalign -i {} -t {} -o {} >/dev/null 2>/dev/null'.format(file1, file2, out)
    os.system(hhalign_command)
    assert os.path.exists(out)
    assert os.path.getsize(out) > 0, "{} is empty".format(out)


def run_hhalign_locally(protein_sequence, profile, fasta_file, hhr_file):
    f = open(fasta_file, "w")
    f.write(">prot\n{}\n".format(protein_sequence))
    f.close()
    hhalign_command(fasta_file, profile, hhr_file)


def run_single_sequence_with_profile(protein_sequence, hmm_profile):
    with tempfile.TemporaryDirectory() as working_local_dir:
        single_fasta = '{}/single_seq.fasta'.format(working_local_dir)
        hhr_file = '{}/single_seq.hhr'.format(working_local_dir)
        tsv_file = '{}/single_seq.tsv'.format(working_local_dir)
        run_hhalign_locally(protein_sequence, hmm_profile, single_fasta, hhr_file)
        hh.convert_file(hhr_file, tsv_file, 100)
        df = pd.read_csv(tsv_file, sep='\t').sort_values('Score', ascending=False)
    return df


def run_single_sequence_with_local_profiles_folder(protein_sequence, profiles_folder):
    hhm_list = glob.glob("{}/*.hhm".format(profiles_folder))
    res_list = []
    for h in hhm_list:
        res = run_single_sequence_with_profile(protein_sequence, h)
        res_list.append(res)
    res = pd.concat(res_list).sort_values(by="Score", ascending=False).reset_index(drop=True)
    return res

