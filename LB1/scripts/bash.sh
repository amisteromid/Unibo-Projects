# Keep unique files based on row 4 in pfam file
awk -F "," '!seen[$4]++' pfam00014.csv >pfam00014.csv.2

#change csv of ids in pdb query to newline seperated file and find rows in pfam which has ids in common with pdb query
cat pdb_query_ids.txt |tr ',' '\n' >pdb_query_ids.txt 
python3 scripts/find_in_common.py pfam00014.csv.2 pdb_query_ids.txt 

# Downloading PDB and Fasta files of the 4th coloumn in pfam csv
for i in $(grep -v "^," pfam00014.csv |cut -d, -f 4); do wget -nc https://files.rcsb.org/view/$i.pdb; done
for i in $(grep -v "^," pfam00014.csv |cut -d, -f 4); do wget -nc https://www.rcsb.org/fasta/entry/$i; mv $i $i.fasta; done

# Iterate through all the fasta files in data directory and trim the according to the kunitz domain and put them in all_kunitz.txt
for i in $(grep -v 'PDB' pfam00014.csv.3 |cut -d, -f1); do python3 scripts/trim_filter_seq.py pfam00014.csv.3 data/$i.fasta > $i.fasta.2; done
for i in data/*.fasta.2; do cat $i >>all_kunitz.fasta; done

# Iterate through all the pdb file and trim them according to the pfam and create *.pdb.2
for f in data/*.pdb; do python3 scripts/trim_pdb.py pfam00014.csv.3 $f; done

# Remove redundancy
grep -v None results/all_kunitz.fasta > results/all_kunitz.fasta.2
cd-hit -i results/all_kunitz.fasta.2 -o results/all_kunitz95.fasta -c 0.95 -n 5

# Grep all the PDB ids in the file after removing redundancy
grep -o -P '(?<=>).*(?=_)' all_kunitz95.fasta >final_PDB_ids.txt

# Adding some string to the beginning of the final PDB ids file
awk '{$0=$0".pdb.2"}1' results/final_PDB_ids.txt 
awk '{print "data/"$0}' result/final_PDB_ids.txt   #write them in a file with different name 0.0

# Preforming the structural alignment
../mTM-align/src/mTM-align -i results/final_PDB_ids.txt 

# Building HMM model out of fasta file resulted from structural alignment
hmmbuild kunitz.hmm mTM_result/result.fasta

# Making a list of uniprot ids from our positive and negative dataset
grep ">" data/uniprot-clean-pf00014.fasta |cut -d "|" -f2 >data/list-clean-pfam00014.txt
grep ">" data/uniprot-not-pf00014.fasta |cut -d "|" -f2 >data/list-not-pfam00014.txt

# Shuffle the uniprot ids and split them for cross validation
sort -R data/list-clean-pfam00014.txt >data/shuffle-clean-pfam00014.txt
sort -R data/list-not-pfam00014.txt >data/shuffle-not-pfam00014.txt
head -n 168 data/shuffle-clean-pfam00014.txt >data/positive-2.txt
tail -n 168 data/shuffle-clean-pfam00014.txt >data/positive-1.txt
head -n 283314 data/shuffle-not-pfam00014.txt >data/negative-1.txt
tail -n 283314 data/shuffle-not-pfam00014.txt >data/negative-2.txt
comm -12 <(sort data/negative-1.txt) <(sort data/negative-2.txt) # to check for any item in common


# Make fasta files out of our splited list of ids
python3 scripts/select-seqs.py data/positive-1.txt data/uniprot-clean-pf00014.fasta >data/positive-1.fasta
python3 scripts/select-seqs.py data/positive-2.txt data/uniprot-clean-pf00014.fasta >data/positive-2.fasta
python3 scripts/select-seqs.py data/negative-1.txt data/uniprot-not-pf00014.fasta >data/negative-1.fasta
python3 scripts/select-seqs.py data/negative-2.txt data/uniprot-not-pf00014.fasta >data/negative-2.fasta

# HMMsearch
hmmsearch -Z 1 --domZ 1 --max --noali --tblout negative-1.match kunitz.hmm negative-1.fasta
hmmsearch -Z 1 --domZ 1 --max --noali --tblout negative-2.match kunitz.hmm negative-2.fasta
hmmsearch -Z 1 --domZ 1 --max --noali --tblout positive-1.match kunitz.hmm positive-1.fasta
hmmsearch -Z 1 --domZ 1 --max --noali --tblout positive-2.match kunitz.hmm positive-2.fasta

# Take the uniprot id, e-value from the match, Then add 0 or 1 according to the negative or positive class
grep -v '#' data/negative-1.match |awk '{print $1,$8,0}' >data/negative-1.class
grep -v '#' data/negative-2.match |awk '{print $1,$8,0}' >data/negative-2.class
grep -v '#' data/positive-2.match |awk '{print $1,$8,1}' >data/positive-2.class
grep -v '#' data/positive-1.match |awk '{print $1,$8,1}' >data/positive-1.class

# Add the ones that are not included in the .match file by hmmsearch and give them a random high evalue of 10
comm -23 <(sort negative-1.txt) <(cut -d"" -f 1 negative-1.class |sort) |awk '{print $1,10,0}' >>negative-1.class
comm -23 <(sort negative-1.txt) <(cut -d"" -f 1 negative-2.class |sort) |awk '{print $1,10,0}' >>negative-2.class

# merge sets of negative and positive
cat negative-1.class positive-1.class >set-1.class
cat negative-2.class positive-2.class >set-2.class

# Evaluate the sets using the confusion script and a given threshold
for i in 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11; do python3 scripts/confusion.py overall_set.class $1; done


# TH = 1.8e-09
