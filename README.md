# CRISPR_type-II_classification

### Goal
This tool classifies type-II nucleases to subtypes.
Given the a nuclease amino acid sequence, its DNA contig, location in the contig, strand, and DB of CRISPR operon HMM profiles, the tool gives a score for its classification into one of the currently existing type II subtypes.

### Preliminary steps: 
There are three preliminary steps prior to running the code:

1. Download and unzip the db file from here https://zenodo.org/record/7248722
2. Create conda env from the yml file https://github.com/emendobio/CRISPR_typeII_nuclease_classification/blob/main/classification.yml 
3. Activate the conda env 

```
wget https://zenodo.org/record/8013752/files/nuclease_classification_db.tar.gz
conda env create --name classification -f classification.yml
conda activate classification
```

### Input parameters: 

- protein_sequence - the full protein sequence of the nuclease (amino acids)
- fasta_file - contig with the sequence of the nuclease (DNA)
- start, end - the coordinates of the nuclease sequence in the fasta file
- strand - the orientation of the nuclease sequence in the fasta file  1 or -1
- db - path to the downloaded db  

׳
### Running the code: 
Run the code using the following command:
```
python crispr_typeII_classification/src/pipeline/classify_nuclease.py \
-protein_id <prot_id> \
-protein_sequence <protein_amino_acids_sequence> \
-fasta_file <DNA_fasta_file> \
-start <nuclease_start>  \
-end <nuclease_end> \
-strand <strand> \
-db <path to db>
 ```

### Output
The output is a json format containing the following fields:

- cas9_score          : reflects the likelihood of the protein to be a Cas9 nuclease
- ruvc_score          : score of the RuvC domain
- hnh_score           : score of the HNH domain (exists in type II nucleases only)
- is_cas9             : True/False
- has_correct_start   : truncated (True/False)
- was_fixed           : True/False
- protein_seq         : AA sequence
- loci_architecture   : type of architecture
- cas9_classification : The subtype classification
- cas9_dist           : Distance from the closest HMM profile
- HNH_Profile         : The HNH signature


  
