import os
import tempfile
import sys
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Align import substitution_matrices
import s3fs

def get_reverse_complement(original_seq):
    seq = Seq(original_seq)
    return str(seq.reverse_complement())


def find_seq_in_fasta(fasta_file, seq):
    fasta_seq = str(get_fasta_seq_from_file(fasta_file).seq).upper()
    return find_seq_in_another_seq(fasta_seq, seq.upper())

def find_seq_in_another_seq(reference_seq, query_seq):
    fasta_ind = reference_seq.find(query_seq)
    if fasta_ind !=- 1:
        return [fasta_ind, fasta_ind+len(query_seq)]
    else:
        rev_seq = get_reverse_complement(query_seq)
        fasta_ind = reference_seq.find(rev_seq)
        if fasta_ind == -1:
            return [-1, -1]
        else:
            return [fasta_ind + len(rev_seq),  fasta_ind]

def get_value_within_bounds(val, min_val, max_val):
    assert min_val < max_val
    if val < min_val:
        return min_val
    elif val > max_val:
        return max_val
    else:
        return val

def get_seq_from_fasta_by_positions(fasta_file, start_pos, end_pos, allow_out_of_bounds = False):
    fasta_seq = str(get_fasta_seq_from_file(fasta_file).seq)
    seq_len = len(fasta_seq)
    if allow_out_of_bounds:
        start_pos = get_value_within_bounds(start_pos, 0, seq_len)
        end_pos = get_value_within_bounds(end_pos, 0, seq_len)
    else:
        assert start_pos == get_value_within_bounds(start_pos, 0, seq_len)
        assert end_pos == get_value_within_bounds(end_pos, 0, seq_len)
    if start_pos == end_pos:
        return ""
    else:
        if start_pos < end_pos:
            return fasta_seq[start_pos:end_pos]
        else:
            rev_seq = fasta_seq[end_pos:start_pos]
            return get_reverse_complement(rev_seq)





def write_to_temp_fasta(query_seq, working_local_dir):
    query_fasta = '{}/temp.fasta'.format(working_local_dir)
    write_fasta_file_from_list(query_fasta, [['temp_fasta', query_seq]])
    return query_fasta

def get_fasta_seq_from_local_file(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    assert len(records) == 1
    return records[0]


def get_fasta_seq_from_file(fasta_file):
    if fasta_file.startswith('s3://'):
        fs = s3fs.S3FileSystem(anon=False)
        with tempfile.TemporaryDirectory() as working_local_dir:
            local_fasta = '{}/sequence.fasta'.format(working_local_dir)
            fs.get(fasta_file, local_fasta)
            fasta_seq = get_fasta_seq_from_local_file(local_fasta)
    else:
        fasta_seq = get_fasta_seq_from_local_file(fasta_file)
    return fasta_seq


def write_fasta_file_from_list(fasta_name, fasta_content_list):
    f = open(fasta_name, "w")
    for fasta_record in fasta_content_list:
        header = fasta_record[0]
        seq = fasta_record[1]
        f.write(">{}\n".format(header))
        f.write("{}\n".format(seq))
    f.close()
    assert os.path.exists(fasta_name)


def create_db(fasta_file, out_dir):
    assert os.path.exists(out_dir)
    assert os.path.exists(fasta_file)
    new_db = '{}/my_db'.format(out_dir)
    out_params = ">{}/out1 2>{}/err1".format(out_dir, out_dir)
    my_cmd = 'makeblastdb  -dbtype nucl -out {} -in {} {}'.format(new_db, fasta_file, out_params)
    os.system(my_cmd)
    assert os.path.exists("{}.nsq".format(new_db))
    assert os.path.exists("{}.nin".format(new_db))
    assert os.path.exists("{}.nhr".format(new_db))
    return new_db


def read_fasta_to_dataframe(fasta_file):
    if fasta_file.startswith('s3://'):
        fs = s3fs.S3FileSystem(anon=False)
        with tempfile.TemporaryDirectory() as working_local_dir:
            local_fasta = '{}/sequence.fasta'.format(working_local_dir)
            fs.get(fasta_file, local_fasta)
            fasta_data = read_local_fasta_to_dataframe(local_fasta)
    else:
        fasta_data = read_local_fasta_to_dataframe(fasta_file)
    return fasta_data


def read_local_fasta_to_dataframe(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    data = {'id':[], 'seq':[]}
    for r in records:
        data['id'].append(r.id)
        data['seq'].append(str(r.seq))
    fasta_data = pd.DataFrame(data)
    return fasta_data


def write_df_to_fasta(df, fasta_file):
    f = open(fasta_file, "w")
    for _, row in df.iterrows():
        f.write(">{}\n{}\n".format(row['id'], row['seq']))
    f.close()


def compare_sequences_with_same_length_by_similarity(seq1, seq2):
    n = len(seq1)
    assert len(seq2) == n
    blosom62 = substitution_matrices.load('BLOSUM62')
    aa_index = {}
    for aa in list(blosom62.alphabet):
        aa_index[aa] = len(aa_index)
    score = 0
    for i in range(n):
        ind1 = aa_index[seq1[i]]
        ind2 = aa_index[seq2[i]]
        score += blosom62[ind1,ind2]
    return score


def compare_sequences_with_same_length_by_identity(seq1, seq2):
    n = len(seq1)
    assert len(seq2) == n
    score = 0
    for i in range(n):
        score += int(seq1[i] == seq2[i])
    return score




