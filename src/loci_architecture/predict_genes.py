"""
import re
"""
import pandas as pd
import subprocess
import os
import glob
import logging
import numpy as np
from Bio import SeqIO
"""
min_score_hmm = 40
"""
max_dist_from_nuclease = 5000
"""
min_gene_length = 200
"""
src_code_path = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1])
logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.INFO)


def parse_hmm_result_line(res_line):
    line_content = [x for x in res_line.split(" ") if x]
    assert len(line_content) == 30
    res = {"ORF": line_content[0],
           "query_accession": line_content[1],
           "tlen": line_content[2],
           #"Hmm": line_content[3],
           "query_accession": line_content[4],
           "qlen": line_content[5],
           "E-value": line_content[6],
           "score": float(line_content[7]),
           "start": int(line_content[23]),
           "end": int(line_content[25]),
           "strand": int(line_content[27])}
    return res


def parse_hmm_result_lines(lines_filtered):
    res_list = [parse_hmm_result_line(res_line) for res_line in lines_filtered]
    res_df = pd.DataFrame(res_list)
    return res_df


def parse_hmm_res_file(hmm_res_file):
    with open(hmm_res_file) as file:
        lines = file.readlines()
    lines_filtered = [x for x in lines if x[0] != "#"]
    if lines_filtered:
        df = parse_hmm_result_lines(lines_filtered)
        df["Hmm"] = hmm_res_file.split("/")[-1].split(".")[0]
    else:
        df = pd.DataFrame(
            columns=["ORF", "query_accession", "tlen", "Hmm", "query_accession", "qlen", "E-value",
                     "score", "start",
                     "end", "strand"])
    return df

def get_cas9_index(res):
    for x in ['Cas9', 'RuvC_III', 'HNH_4','spcas9_818_875','spcas9_975_999']:
        ind = get_index_of_gene(res, x)
        if ind != -1:
            return ind
    return -1


def match_pattern(hmm_res, curr_pattern, cas9_index, cas9_strand):
    if curr_pattern['position'] == 'downstream':
        is_matching =   check_pattern(hmm_res, cas9_index, cas9_strand, curr_pattern['genes'], curr_pattern['direction'], 1)
    else:
        assert curr_pattern['position'] == 'upstream'
        is_matching = check_pattern(hmm_res, cas9_index, cas9_strand, curr_pattern['genes'], curr_pattern['direction'], -1)
    return is_matching


def check_gene(hmm_res, index, expected_gene_name, expected_gene_strand):
    n = hmm_res.shape[0]
    if index < 0 or index >= n:
        return False
    match_gene_name = hmm_res.iloc[index]['Hmm'].startswith(expected_gene_name)
    match_gene_strand = expected_gene_strand == hmm_res.iloc[index]['Strand']
    return match_gene_name and match_gene_strand


def check_pattern(hmm_res, cas9_index, cas9_strand, genes_list_str, genes_direction_str, direction):
    if genes_list_str is None and genes_direction_str is None:
        return True
    direction_int = {'+':1, '-':-1}
    genes_list = genes_list_str.split(',')
    direction_list = [cas9_strand*direction_int[x] for x in genes_direction_str.split(',')]
    n = len(genes_list)
    assert len(direction_list) == n
    for i in range(n):
        curr_index = cas9_index + (i+1)*cas9_strand*direction
        if not check_gene(hmm_res, curr_index, genes_list[i], direction_list[i]):
            return False
    return True


def match_all_patterns(hmm_res, patterns_df):
    if hmm_res.empty:
        return "NA"
    cas9_index = get_cas9_index(hmm_res)
    if cas9_index == -1:
        return "NA"
    cas9_strand = int(hmm_res.iloc[cas9_index]['Strand'])
    for index, curr_pattern in patterns_df.iterrows():
        if match_pattern(hmm_res, curr_pattern, cas9_index, cas9_strand):
            return curr_pattern['pattern_name']
    return "Cas9"


def get_type_II_patterns():
    patterns_file = '{}/loci_architecture/type_II_patterns.tsv'.format(src_code_path)
    patterns_df = pd.read_csv(patterns_file, sep=' ')
    return patterns_df


def get_index_of_gene(res, gene_name):
    cas9_index_list = np.where(res['Hmm'].str.startswith(gene_name))[0]
    if len(cas9_index_list) == 0:
        return -1
    else:
        scores_of_gene = res.loc[cas9_index_list]['score']
        first_orf_name = res.loc[cas9_index_list]['ORF'].iloc[0]
        if len(cas9_index_list) > 1:
            print("{} multiple matches with scores {}".format(first_orf_name, list(scores_of_gene)))
        cas9_index = scores_of_gene.idxmax()
        return cas9_index
"""
def filter_best_match(annotated_genes):
    prev_orf = ''
    best_match = []
    for _,row in annotated_genes.iterrows():
        if row['ORF'] != prev_orf:
            best_match.append(True)
        else:
            best_match.append(False)
        prev_orf = row['ORF']
    best_match_genes = annotated_genes[best_match]
    return best_match_genes.reset_index(drop=True)



def get_genes_hmm_data(curr_folder, start_pos, end_pos):

    genes_file = '{}/genes.tab'.format(curr_folder)
    hmmer_file = '{}/hmmer.tab'.format(curr_folder)
    genes_data = pd.read_csv(genes_file, sep='\t')
    hmmer_data = pd.read_csv(hmmer_file, sep='\t')
    I = (genes_data['Start'] >= start_pos - max_dist_from_nuclease) & (genes_data['End'] <= end_pos + max_dist_from_nuclease) & (genes_data['End']  - genes_data['Start'] >=  min_gene_length)
    genes_data = genes_data[I]
    I = hmmer_data['score'] >= min_score_hmm
    hmmer_data = hmmer_data[I]
    hmmer_data = pd.concat([hmmer_data, pd.DataFrame({'ORF': list(genes_data['ORF']), 'Hmm': 'None', 'score': 0})])
    annotated_genes = pd.merge(genes_data, hmmer_data, on='ORF').sort_values(by=['Start', 'score'],
                                                                             ascending=[True, False])
    return annotated_genes


def read_annotated_genes_from_folder(curr_folder, start_pos=0, end_pos=np.inf):
    annotated_genes = get_genes_hmm_data(curr_folder, start_pos, end_pos)
    if annotated_genes.empty:
        return annotated_genes
    best_match_genes = filter_best_match(annotated_genes)
    assert best_match_genes['score'].min() >=0
    return best_match_genes.reset_index(drop=True)




"""
def read_orf_list_from_prodigal_fasta(prodigal_out_fasta):
    records = list(SeqIO.parse(prodigal_out_fasta, "fasta"))
    if len(records) == 0:
        return pd.DataFrame()
    res_list = []
    for r in records:
        header_info = r.description.split(' # ')
        res = {'ORF': header_info[0], 'Start': int(header_info[1]), 'End': int(header_info[2]),
               'Strand': int(header_info[3]), 'Sequence':str(r.seq)}
        res_list.append(res)
    return pd.DataFrame(res_list)


def write_orf_list(prodigal_out_fasta, orf_list_file):
    required_columns = ['ORF','Start','End','Strand']
    df = read_orf_list_from_prodigal_fasta(prodigal_out_fasta)
    if df.empty:
        return False
    assert all(col in df.columns for col in required_columns)
    df = df[required_columns]
    df.to_csv(orf_list_file, sep='\t', index=False)
    return True


# Load data
def aggregate_hmm_results(hmmer_results_folder, aggregated_hmm_file):
    logging.info('Loading HMMER output')

    # Get files
    hmm_files = glob.glob("{}/*.tab".format(hmmer_results_folder))
    res_list = [parse_hmm_res_file(hmm_res_file) for hmm_res_file in hmm_files]
    non_empty_res_list = [df for df in  res_list if not df.empty]
    if not non_empty_res_list:
        return pd.DataFrame()
    hmm_df = pd.concat(non_empty_res_list)
    hmm_df.reset_index(drop=True, inplace=True)
    return hmm_df.drop_duplicates()
"""
def read_hmm_result(in_file):
    hmm_results = []
    with open(in_file, "r") as f:
        line = f.readline()
        while line:
            if not line.startswith('#'):
                line_parts = re.split('\s+', line.rstrip())
                hmm_results.append({'ORF':line_parts[0],'ORF_accession':line_parts[1], 'HMM':line_parts[3],'HMM_accession':line_parts[4], 'Score':float(line_parts[7])})
            line = f.readline()
    return pd.DataFrame(hmm_results)


"""
def run_prodigal(dna_fasta, proteins_fasta, orf_list_file, prodigal_log_file):
    with open(prodigal_log_file, 'w') as prodigal_log:
        prod_cmd = ['prodigal' ,'-i', dna_fasta, '-a', proteins_fasta, '-p', 'meta']
        subprocess.run(prod_cmd, stdout=subprocess.DEVNULL, stderr=prodigal_log)
        assert os.path.exists(proteins_fasta), "prodigal error: {} wasn't created from {}".format(proteins_fasta, dna_fasta)
        return write_orf_list(proteins_fasta, orf_list_file)


def run_single_hmm_profile(proteins_fasta, hmm_profile, hmm_output_file, log_file):
    with open(log_file, 'a') as hmmer_log:
        hmm_cmd = ['hmmsearch', '--domtblout', hmm_output_file, hmm_profile, proteins_fasta]
        subprocess.run(hmm_cmd, stdout=subprocess.DEVNULL, stderr=hmmer_log)
        assert os.path.exists(hmm_output_file)


def get_profiles_list(hmm_folder):
    profiles_list = glob.glob('{}/*.HMM'.format(hmm_folder)) + glob.glob('{}/*.hmm'.format(hmm_folder))
    return profiles_list



def run_all_profiles_in_directory(hmm_folder, proteins_fasta, work_dir):
    profiles_list = get_profiles_list(hmm_folder)
    hmm_results_folder = '{}/hmmer'.format(work_dir)
    os.mkdir(hmm_results_folder)
    log_file = "{}/hmmer.log".format(work_dir)
    for hmm_profile in profiles_list:
        hmm_out = hmm_profile.split('/')[-1][:-4]
        hmm_output_file = "{}/{}.tab".format(hmm_results_folder, hmm_out)
        run_single_hmm_profile(proteins_fasta, hmm_profile, hmm_output_file, log_file)
    aggregated_hmm_file = "{}/temp_hmmer.tsv".format(work_dir)
    hmm_results = aggregate_hmm_results(hmm_results_folder, aggregated_hmm_file)
    return hmm_results


def run_prodigal_on_contig(fasta_file, working_dir):
    proteins_fasta = '{}/proteins.faa'.format(working_dir)
    prodigal_log = '{}/prodigal.log'.format(working_dir)
    orf_list_file = '{}/genes.tab'.format(working_dir)
    run_prodigal(fasta_file, proteins_fasta, orf_list_file, prodigal_log)
    assert os.path.exists(orf_list_file)
    assert os.path.exists(proteins_fasta)
    return proteins_fasta
"""

def get_hmm_results(proteins_fasta, profiles_folder, nuclease_pos, working_dir):
    hmmer_out_file = '{}/hmmer.tab'.format(working_dir)
    hmm_results = run_all_profiles_in_directory(profiles_folder, proteins_fasta, working_dir)
    hmm_results.to_csv(hmmer_out_file, sep='\t', index=False)
    pos_arr = [int(x) for x in nuclease_pos.split("-")]
    res = read_annotated_genes_from_folder(working_dir, start_pos=pos_arr[0], end_pos=pos_arr[1])
    return res
"""

def get_hmm_results_from_prodigal_fasta(proteins_fasta, position_string, profiles_folder, working_dir, full_genes_list=True,
                                        max_dist=max_dist_from_nuclease):
    """
        Proforms HMM search on a Prodigal fasta file.

        Args:
            fasta_file (str): The path to the fasta file resulting from Prodigal.

            position_string (str): The position of the nuclease in the fasta file.

            profiles_folder (str): The path to the local folder containing the HMM profiles.
            The profiles are located in s3://emendobio-compute/L2/nuclease_classification_db/loci_profiles
            and should be downloaded to a local folder

            full_gened_list (bool): If True, returns all genes found in the region, not only the ones that match the pattern.

            max_dist_from_nuclease (int): The maximum distance from the nuclease to search for genes.


        Returns:
            pandas.DataFrame: A dataframe containing the results of the HMM search.

        """
    pos = [int(x) for x in position_string.split('-')]
    hmm_results = run_all_profiles_in_directory(profiles_folder, proteins_fasta, working_dir)
    I = (hmm_results["start"] >= pos[0] - max_dist) & (
                hmm_results["end"] <= pos[1] + max_dist)
    hmm_results = hmm_results[I].sort_values("score", ascending=False).drop_duplicates("ORF").sort_values("start",
                                                                                                          ascending=True)
    genes_data = pd.read_csv('{}/genes.tab'.format(working_dir), sep='\t')
    if full_genes_list:
        hmm_results = pd.merge(genes_data, hmm_results[["ORF", "Hmm", "score"]], how='left').fillna(
            {'Hmm': "None", 'score': 0})
    else:
        hmm_results = pd.merge(genes_data, hmm_results[["ORF", "Hmm", "score"]])
    return hmm_results
