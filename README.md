# CRISPR_type-II_classification


### preliminary steps: 
Classification of type-II nucleases to subtypes

There are three preliminary steps prior to running the code:

1. Download and unzip the db file from here https://zenodo.org/record/7248722
2. Create conda env from the yml file https://github.com/emendobio/CRISPR_typeII_nuclease_classification/blob/main/classification.yml 
3. Activate the conda env 

```
wget https://zenodo.org/record/7248722/files/nuclease_classification_db.tar.gz

command conda create --name classification classification.yml

conda activate classification

```

### input parameters: 

- protein_sequence - the full protein sequence of the nuclease (amino acids)
  
- fasta_file - contig with the sequence of the nuclease (DNA)
  
- start, end - the coordinates of the nuclease sequence in the fasta file

- strand - the orientation of the nuclease sequence in the fasta file  1 or -1
  
- db - path to the downloaded db  


׳
### running the code: 

run the code using the following command:
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
  

  
