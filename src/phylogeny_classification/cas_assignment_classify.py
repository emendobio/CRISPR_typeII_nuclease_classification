import json
import subprocess
from Bio import SeqIO
import math
import re
import sys
import os
src_code_path = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1])
sys.path.append(src_code_path)
from loci_architecture import predict_genes as pg
assert sys.version_info >= (3, 5)


def run_command(command, stdout=None, stderr=None):
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

def classify_single_record(clade_vals, record, seq2clade, learned_data_directory, working_directory, single_res, interactive):
    record_results = []
    print(record.name)
    results = search_vs_clade_profiles(learned_data_directory, working_directory, record, clade_vals.keys(), seq2clade, interactive)
    if not results:
        return record_results
    conf = {}
    best_c = None
    best_c_val = 0
    for c in results:
        # print(c, results[c]['d'], clade_vals[c]['internal_max'], clade_vals[c]['external_min'])
        conf[c] = 0
        if results[c]['d'] < clade_vals[c]['internal_max']:
            cas_found = True
            conf[c] = 1
            if best_c_val < 1:
                best_c = c
                best_c_val = 1
            if results[c]['d'] < clade_vals[c]['external_min']:
                conf[c] = 2
                best_c = c
                best_c_val = 2
            if not single_res:
                new_result = {'name': record.name, 'clade': c, 'Koonin_distance': results[c]['d'],
                              'confidence': conf[c], 'P-value': results[c]['P-value'],
                              'Coverage': results[c]['Coverage']}
                record_results.append(new_result)
                # print('\t'.join([]), file=ofh)
    if best_c_val == 0 or single_res:
        # Write the smallest p-value hit if nothing found
        if best_c_val == 0:
            assert len(results.keys()), results
            best_c = sorted(results.keys(), key=lambda item: float(results[item]['P-value']))[0]
        new_result = {'name': record.name, 'clade': best_c, 'Koonin_distance': results[best_c]['d'],
                      'confidence': conf[best_c], 'P-value': results[best_c]['P-value'],
                      'Coverage': results[best_c]['Coverage']}
        record_results.append(new_result)
    return record_results


def search_vs_clade_profiles(learn_dir, wd, record, clades, seq2clade, interactive=True):
    name = re.sub(r'[\\\/\|]', '_', record.name)
    if interactive:
        print(f'classifying {name}')
    with open(f'{wd}/{name}.fa', 'w') as ofh:
        print(f'>{name}\n{record.seq}\n', file=ofh)

    run_command(f'hhalign -i {wd}/{name}.fa -t {wd}/{name}.fa', stdout=f'{wd}/{name}.self.hhalign', stderr=f'{wd}/{name}.self.err')
    self_seq_score = extract_hhalign_score(f'{wd}/{name}.self.hhalign')

    self_clade_scores = {}
    for c in clades:
        self_clade_scores[c] = extract_hhalign_score(f'{learn_dir}/{c}.{c}.hhalign')

    run_command(
        f'hhsearch -i {wd}/{name}.fa -d {learn_dir}/hhdb',
        stdout=f'{wd}/{name}.hhsearch', stderr=f'{wd}/{name}.err')

    hhsearch_results = parse_hhsearch_results(f'{wd}/{name}.hhsearch', seq2clade)

    # Normalize hhsearch results
    for c in hhsearch_results:
        hhsearch_results[c].update({'SelfScore': self_seq_score['Score'], 'SelfCladeScore': self_clade_scores[c]['Score']})
        hhsearch_results[c]['d'] = -math.log(float(hhsearch_results[c]['Score'])/min(float(self_seq_score['Score']), float(self_clade_scores[c]['Score'])))
    return hhsearch_results


def parse_hhsearch_results(results_file, seq2clade):
    hhsearch_results = {}
    table_flag = False
    with open(results_file) as ifh:
        for line in ifh:
            if line == '\n' and table_flag:
                table_flag = False
            # do the actual work here:
            if table_flag:
                line_splitted = line[0:34].split()
                try:
                    hit_name = seq2clade[line_splitted[1]]
                except KeyError:
                    # logging.error("could not find key {} in seq2clade: {}".format(line_splitted[1], seq2clade))
                    raise
                if hit_name not in hhsearch_results:
                    hhsearch_results.update({hit_name: parse_hhresults_line(line)})
            if line.startswith(' No Hit'):
                table_flag = True
    return hhsearch_results


def parse_hhresults_line(line):
    results = {}
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
    results['Coverage'] = line[69:74]
    return results


def get_th_parameters(pars):
    th_params = {}
    seq2clade = {}
    for c in pars:
        th_params[c] = {}
        internal_scores = [s['Score'] for s in pars[c]['internal']]
        internal_dists = [-math.log(float(s['Score'])/min(float(s['SelfCladeScore']), float(s['SelfScore']))) for s in pars[c]['internal']]
        external_scores = [s['Score'] for s in pars[c]['external']]
        external_dists = [-math.log(float(s['Score'])/min(float(s['SelfCladeScore']), float(s['SelfScore']))) for s in pars[c]['external']]
        th_params[c]['internal_min'] = min(internal_dists)
        th_params[c]['internal_max'] = max(internal_dists)
        th_params[c]['external_min'] = min(external_dists)
        th_params[c]['external_max'] = max(external_dists)

        seq2clade.update({s['query'][:30]: c for s in pars[c]['internal']})  # changed for a max of 30 chars - 19.10.2020

    return th_params, seq2clade


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
            if line=='' and hits_flag:
                hits_flag = False
            if hits_flag:
                results.update(parse_hhresults_line(line))
            if line.startswith(' No Hit'):
                hits_flag = True
    return results

def get_record_of_orf_from_fasta(hmm_results, gene, proteins_fasta):
    gene_index = pg.get_index_of_gene(hmm_results, gene)
    orf = hmm_results.iloc[gene_index]["ORF"]
    with open(proteins_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == orf:
                return record
    return None


def classify_sequence_record(record, model_folder, working_dir):
    with open(f'{model_folder}/data.json') as jfh:
        pars = json.load(jfh)
    clade_vals, seq2clade = get_th_parameters(pars)
    classification = classify_single_record(clade_vals, record, seq2clade, model_folder, working_dir, True, True)
    if classification:
        return classification[0]
    else:
        return ""
