from utils import hmm_utils as hu
from utils import seq_utils as su
D_position_in_hmm = 13
min_len_before_D = 7
query_max_pos = 50
hnh_min_score = 100
ruvc_min_score = 150

def get_score_of_feature(hmm_res, prefix):
    I = hmm_res['annotation'].str.startswith(prefix)
    if I.sum() == 0:
        return 0
    else:
        best_match_score =  hmm_res[I].iloc[0]['Score']
        assert hmm_res[I]['Score'].max() == best_match_score, (hmm_res[I]['Score'])
        return best_match_score


def analyze_scores_of_nuclease(hmm_res):
    nuclease_score = get_score_of_feature(hmm_res, 'CAS-II')
    ruvc_score = get_score_of_feature(hmm_res, 'RuvC')
    hnh_score = get_score_of_feature(hmm_res, 'HNH')
    return [nuclease_score, ruvc_score, hnh_score]


def test_if_cas9_has_correct_start(hmm_res):
    row, _ = get_best_start_match(hmm_res)
    if row is None:
        return False
    else:
        d_curr_pos = D_position_in_hmm - row['TEMPLATEstart'] + row['QUERYstart']
        passed_criteria = row['QUERYstart'] <= query_max_pos and d_curr_pos >= min_len_before_D
        return passed_criteria

def get_best_start_match(hmm_res):
    cas9_start_hmm = '1-24'
    min_hmm_score = 25

    I = (hmm_res['annotation'] == cas9_start_hmm) & (hmm_res['Score'] >= min_hmm_score)
    if I.sum() == 0:
        return [None, None]
    row = hmm_res[I].iloc[0]
    assert row['Score'] == hmm_res[I]['Score'].max()
    conserved_d_pos = row['QUERYstart']  + (D_position_in_hmm - row['TEMPLATEstart'] )
    return [row, conserved_d_pos]


def check_nuclease_essential_parts(protein_seq, db_name):
    hmm_res = hu.run_single_sequence(protein_seq, db_name)
    [nuclease_score, ruvc_score, hnh_score] = analyze_scores_of_nuclease(hmm_res)
    has_correct_start = test_if_cas9_has_correct_start(hmm_res)
    return [has_correct_start, nuclease_score, ruvc_score, hnh_score]

def get_max_len_three_multiple(n):
    return 3*int(n/3)


def get_extended_seq(fasta_name, start_pos, end_pos, strand):
    seq = su.get_fasta_seq_from_file(fasta_name)
    if strand == 1:
        extended_len = get_max_len_three_multiple(start_pos)
        new_pos = start_pos - extended_len
        extended_prot_seq = seq[new_pos:start_pos].translate()
    else:
        assert strand == -1
        extended_len = get_max_len_three_multiple(len(seq)-end_pos)
        new_pos = end_pos + extended_len
        extended_prot_seq = seq[end_pos:new_pos].reverse_complement().translate()
    prot_seq = str(extended_prot_seq.seq)
    stop_codon_ind = prot_seq.rfind("*")
    if stop_codon_ind != -1:
        prot_seq = prot_seq[stop_codon_ind+1:]
    return prot_seq


def get_fixed_version_of_protein(curr_protein_sequence, fasta_file, curr_protein_start, curr_protein_end, strand, db_name):
    prot_seq_ex = get_extended_seq(fasta_file, curr_protein_start, curr_protein_end, strand)
    if not prot_seq_ex or "M" not in prot_seq_ex:
        return "missing_M"
    extended_protein_seq = '{}{}'.format(prot_seq_ex, curr_protein_sequence)
    res1 = hu.run_single_sequence(extended_protein_seq, db_name)
    row, d_pos = get_best_start_match(res1)
    if row is None:
        return "missing_start_match"
    fixed_prot_start = extended_protein_seq[:(d_pos - min_len_before_D)].rfind("M")
    if fixed_prot_start >= len(prot_seq_ex) or fixed_prot_start == -1:
        return "missing_bad_position"
    else:
        fixed_protein = extended_protein_seq[fixed_prot_start:]
        assert len(fixed_protein) > len(curr_protein_sequence), extended_protein_seq[:(d_pos - min_len_before_D)]
        return fixed_protein

"""
This is the main function of nuclease composition
It returns following outputs:

'nuclease_score':  the score of general cas9 profile

'ruvc_score': the score of hhm alignemnt to one of PFAM's RuvC profiles

 'hnh_score': the score of hhm alignemnt to one of PFAM's HNH profiles
 
 'is_cas9':  if the hnh and ruvc scores passed our threshold. We can also add in the future match to Cas12
 
 'has_correct_start': if the protein starts in the RuvC-1 - based on a profile that i have built - in average 84% of the nucleases have correct start
 
  'was_fixed':if it wasn't correct but we were able to restore it - only ~3% of the cases can be restores
  
  'protein_seq': uaually the original seq, unless it was fixed and then its the fixed version
              
"""
def get_nuclease_composition_info(protein_seq, db,  fasta_file, start_pos, end_pos, strand):
    [has_correct_start, nuclease_score, ruvc_score, hnh_score] = check_nuclease_essential_parts(protein_seq, db)
    is_nuclease = ruvc_score > ruvc_min_score and hnh_score > hnh_min_score
    was_fixed = False
    if is_nuclease:
        if not has_correct_start:
            fixed_seq = get_fixed_version_of_protein(protein_seq, fasta_file, start_pos, end_pos, strand, db)
            if not fixed_seq.startswith("missing"):
                was_fixed = True
                protein_seq = fixed_seq
                hmm_res = hu.run_single_sequence(protein_seq, db)
                [nuclease_score, ruvc_score, hnh_score] = analyze_scores_of_nuclease(hmm_res)
    nuclease_composition_info = {'cas9_score': nuclease_score, 'ruvc_score': ruvc_score,
              'hnh_score': hnh_score, 'is_cas9':is_nuclease, 'has_correct_start': has_correct_start,
              'was_fixed':was_fixed, 'protein_seq':protein_seq}
    return nuclease_composition_info
