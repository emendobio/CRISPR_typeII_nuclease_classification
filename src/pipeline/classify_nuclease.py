import os
import sys
import argparse
import tempfile
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
src_code_path = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1])
sys.path.append(src_code_path)
from loci_architecture import predict_genes as pg
from nuclease_composition import nuclease_composition_info as nci
from phylogeny_classification import cas_assignment_classify as cac
from utils import hmm_utils as hu


def get_classifications_data(fasta_file, nuclease_position, profiles_folder, cas9_model, cas1_model, cas9_seq):
    with tempfile.TemporaryDirectory() as working_dir:
        proteins_fasta = pg.run_prodigal_on_contig(fasta_file, working_dir)
        print(proteins_fasta, profiles_folder, nuclease_position, working_dir)
        hmm_results = pg.get_hmm_results(proteins_fasta, profiles_folder, nuclease_position, working_dir)
        patterns = pg.get_type_II_patterns()
        loci_architecture_classification = pg.match_all_patterns(hmm_results, patterns)
        cas9_record = SeqRecord(Seq(cas9_seq),id="cas9_nuclease", description="cas9_nuclease",name="cas9_nuclease")
        cas1_record = cac.get_record_of_orf_from_fasta(hmm_results, "Cas1", proteins_fasta)
        print(cas1_record.seq)
        cas1_classification = cac.classify_sequence_record(cas1_record, cas1_model, working_dir)
        cas9_classification = cac.classify_sequence_record(cas9_record, cas9_model, working_dir)
    res_map = {'loci_architecture':[loci_architecture_classification]}
    if cas1_classification:
        res_map.update({'cas1_classification': cas1_classification['clade'], 'cas1_dist': cas1_classification['Koonin_distance']})
    if cas9_classification:
        res_map.update({'cas9_classification': cas9_classification['clade'], 'cas9_dist': cas9_classification['Koonin_distance']})
    return res_map

def get_hnh_data(cas9_seq, db):
    min_hnh_score = 55
    res = hu.run_single_sequence_with_profiles_folder(cas9_seq, db)
    if res.empty:
        return {'HNH_Profile':'None'}
    best_match = res.iloc[0]
    if best_match['Score'] < min_hnh_score:
        return {'HNH_Profile': 'None'}
    return {'HNH_Profile': best_match['annotation']}


"""
This is the main function that preforms all the analysis on the nuclease candidate
input parameters: 

protein_sequence - the full protein sequence of the nuclease

fasta_file - contig with the coding DNA sequence of the nuclease

start, end - the coordinates and orientation of the nuclease sequence in the fasta file - 
useful in cases where there is more than one CRISPR per contig, or fixing protein start position

strand - used only for cases when need to fix protein start position

db - path to the downloaded db file 
"""

def run_full_classification_analysis_on_candidate(protein_sequence, fasta_file, start, end, strand, db):
    assert os.path.exists(db)
    nuclease_result = {}
    nuclease_composition_info = nci.get_nuclease_composition_info(protein_sequence,
                                                              "{}/cas9_db".format(db),
                                                              fasta_file, start,
                                                              end, strand)
    nuclease_result.update(nuclease_composition_info)
    if nuclease_result['is_cas9']:
        print("Protein is a nuclease, classifying..")
        classification_data = get_classifications_data(fasta_file, "{}-{}".format(start, end),
                                                       "{}/loci_profiles".format(db),
                                                       "{}/phylogeny_model_cas9".format(db),
                                                       "{}/phylogeny_model_cas1".format(db),
                                                       nuclease_result['protein_seq'])
        nuclease_result.update(classification_data)
        hnh_data = get_hnh_data(nuclease_result['protein_seq'], "{}/hnh_profiles".format(db))
        nuclease_result.update(hnh_data)
    else:
        print("Protein is not a nucelase")
    print(nuclease_result)
    return nuclease_result


def main():
    run_full_classification_analysis_on_candidate(args.protein_sequence, args.fasta_file, args.start, args.end, args.strand, args.db)

if __name__ == "__main__":
    #required parameters
    parser = argparse.ArgumentParser(description="Run CRISPR Type-II Classification")
    parser.add_argument("-protein_sequence", help="The full protein sequence of the nuclease", required=True, type=str)
    parser.add_argument("-fasta_file", help="path to fasta file of the contig containing the nuclease DNA sequence", required=True)
    parser.add_argument("-start", help="start pos of the nuclease in the fasta file", required=True, type=int)
    parser.add_argument("-end", help="end pos of the nuclease in the fasta file", required=True, type=int)
    parser.add_argument("-strand", help="strand of the nuclease in the fasta file, should be 1 or -1", required=True, choices=['1', '-1'])
    parser.add_argument("-db", help="path to location of HMM folder", required=True)
    args = parser.parse_args()
    assert os.path.exists(args.fasta_file)
    assert os.path.exists(args.db)
    main()
