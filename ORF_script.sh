#!/bin/bash

# taking in genome ID as input
echo "Enter the genome ID (file) for the genome to be analyzed:"
read GID
# testing GID was 100226.110

# fetching the .fna file from bvbrc
wget ftp://anonymous:xxx@ftp.bvbrc.org/genomes/$GID/$GID.fna

# removing metadata to create the forward string
# replacing new lines with @ | getting the second field (the sequence) | extracting the target string | deleting '@' | keeping everything that doesn't start with $ > saving it into a txt file
tr '\n' '@' < $GID.fna | sed 's/>/\n&/g' | cut -d "@" -f 2- | head -n 2 | tail -n 2 | tr -d '\n' | tr -d '@' | grep -v '^$' > forward.txt
rm $GID.fna

# getting reverse from forward
tr 'a' 'x' < forward.txt | tr 't' 'a' | tr 'c' 'y'| tr 'g' 'c' | tr 'x' 't' | tr 'y' 'g' | rev > rev.txt

# creating reading frames for forward strand
# creating triplets with sed | deleting newlines
sed 's/.\{3\}/& /g' forward.txt | tr -d '\n' > framefor1.txt
cut -c 2- forward.txt | sed 's/.\{3\}/& /g' | tr -d '\n' > framefor2.txt # starting from position 2
cut -c 3- forward.txt | sed 's/.\{3\}/& /g' | tr -d '\n' > framefor3.txt # starting from position 3

# creating reading frames for reverse strand (same process as creating RFs for forward strand)
sed 's/.\{3\}/& /g' rev.txt | tr -d '\n' > framerev1.txt
cut -c 2- rev.txt | sed 's/.\{3\}/& /g' | tr -d '\n' > framerev2.txt
cut -c 3- rev.txt | sed 's/.\{3\}/& /g' | tr -d '\n' > framerev3.txt

# splitting into lines by STOP codons | filtering for lines containing 'atg' | removing sequences before 'atg' and after the stop codons for each line
sed -E "s/tag|taa|tga/&\n/g" framefor1.txt | grep 'atg' | sed -E 's/.*?(atg.*?(taa|tag|tga)).*/\1/' | grep '^atg' | nl > ff1.txt
sed -E "s/tag|taa|tga/&\n/g" framefor2.txt | grep 'atg' | sed -E 's/.*?(atg.*?(taa|tag|tga)).*/\1/' | grep '^atg' | nl > ff2.txt
sed -E "s/tag|taa|tga/&\n/g" framefor3.txt | grep 'atg' | sed -E 's/.*?(atg.*?(taa|tag|tga)).*/\1/' | grep '^atg' | nl > ff3.txt
sed -E "s/tag|taa|tga/&\n/g" framerev1.txt | grep 'atg' | sed -E 's/.*?(atg.*?(taa|tag|tga)).*/\1/' | grep '^atg' | nl > fr1.txt
sed -E "s/tag|taa|tga/&\n/g" framerev2.txt | grep 'atg' | sed -E 's/.*?(atg.*?(taa|tag|tga)).*/\1/' | grep '^atg' | nl > fr2.txt
sed -E "s/tag|taa|tga/&\n/g" framerev3.txt | grep 'atg' | sed -E 's/.*?(atg.*?(taa|tag|tga)).*/\1/' | grep '^atg' | nl > fr3.txt

# labels (one per line) | then replacing tabs with a newline
sed 's/^\ \ */>FR1_ORF_/g' ff1.txt | tr '\t' '\n' > compiled_ORFs0.txt
sed 's/^\ \ */>FR2_ORF_/g' ff2.txt | tr '\t' '\n' >> compiled_ORFs0.txt
sed 's/^\ \ */>FR3_ORF_/g' ff3.txt | tr '\t' '\n' >> compiled_ORFs0.txt
sed 's/^\ \ */>FR4_ORF_/g' fr1.txt | tr '\t' '\n' >> compiled_ORFs0.txt
sed 's/^\ \ */>FR5_ORF_/g' fr2.txt | tr '\t' '\n' >> compiled_ORFs0.txt
sed 's/^\ \ */>FR6_ORF_/g' fr3.txt | tr '\t' '\n' >> compiled_ORFs0.txt

# translation nucleotide -> amino acid (STOP codon is represented by '*')
sed 's/atg/M/g' compiled_ORFs0.txt | sed 's/gca\|gcc\|gcg\|gct/A/g' | sed 's/tgc\|tgt/C/g' | sed 's/gac\|gat/D/g' | sed 's/gaa\|gag/E/g' |  sed 's/ttc\|ttt/F/g' | sed 's/gga\|ggc\|ggg\|ggt/G/g' | sed 's/cac\|cat/H/g' | sed 's/ata\|atc\|att/I/g' | sed 's/aaa\|aag/K/g' | sed 's/cta\|ctc\|ctg\|ctt\|tta\|ttg/L/g' | sed 's/aac\|aat/N/g' | sed 's/cca\|ccc\|ccg\|cct/P/g' | sed 's/caa\|cag/Q/g' | sed 's/aga\|agg\|cga\|cgc\|cgg\|cgt/O/g' | sed 's/agc\|agt\|tca\|tcc\|tcg\|tct/S/g' | sed 's/aca\|acc\|acg\|act/T/g' | sed 's/gta\|gtc\|gtg\|gtt/V/g' | sed 's/tgg/W/g' | sed 's/tac\|tat/Y/g' | sed 's/taa\|tag\|tga/*/g' | tr -d ' ' > translated_ORFs.faa

# formatting compiled_ORFs.txt by removing spaces
cat compiled_ORFs0.txt | tr -d ' ' > compiled_ORFs.fnn

# removing intermediate files
rm forward.txt rev.txt ff1.txt ff2.txt ff3.txt fr1.txt fr2.txt fr3.txt framefor1.txt framefor2.txt framefor3.txt framerev1.txt framerev2.txt framerev3.txt

# counting number of ORFs by line count
let total=$(wc -l compiled_ORFs0.txt)
# dividing by 2 to account for the label line
let n=$(expr $total / 2)
rm compiled_ORFs0.txt

# writing results to the console
echo -e 'fna formatted file outputted to: compiled_ORFs.fna'
echo -e 'faa formaated file outputted to: translated_ORFs.faa'
echo -e "The number of ORFs for the genome ${GID} is ${n}"


