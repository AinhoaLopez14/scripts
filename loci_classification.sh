#Loci counting
##Author: Ainhoa Lopez
##Date: 01 April 2024
#the next command only needs to be executed once
#it's very important to give the output a different name from the genome!!!!
/usr/local/ncbi/blast/bin/makeblastdb -in /Users/ainhoalopez/Documents/pliv/ncbi_dataset/data/GCA_940671915.1/GCA_940671915.1_Pliv_v1_genomic.fna -dbtype nucl -out DB
#the next commands needs to be excuted on the path were the db files are!!!!!
cd ~/Documents/pliv/ncbi_dataset/data/GCA_940671915.1/
/usr/local/ncbi/blast/bin/blastn -query ~/Documents/pliv/ncbi_dataset/data/GCA_940671915.1/sequencies_paracentrotus_clean.fa -db DB -outfmt 6 -evalue 1e-4 -out ~/Documents/pliv/ncbi_dataset/data/GCA_940671915.1/blast_pliv_dhuerta.out

#COUNTING MAPPED AND NOT MAPPED 3730 LOCI
#We have two types of loci: total loci and candidate loci. Candidate loci are the ones under selection, and are included into total
##TOTAL##
#Mapped
cut -f1 blast_pliv_dhuerta.out |sort| uniq|wc -l

#Unmapped
cut -f1 blast_pliv_dhuerta.out |sort| uniq > mapped_dh.list
grep -w -v -f mapped_dh.list sequencies_paracentrotus_clean.fa |grep '>'|wc -l

#Uniq
cut -f1 blast_pliv_dhuerta.out|sort|uniq -c |grep ' 1 ' |wc -l
#Repeated

cut -f1 blast_pliv_dhuerta.out|sort|uniq -c |grep -v ' 1 ' |wc -l 

##CANDIDATE##: 402
#Mapped
cut -f1 blast_pliv_dhuerta.out |sort| uniq| grep -w -f candidates.list |wc -l 

#Unmapped
cut -f1 blast_pliv_dhuerta.out |sort| uniq| grep -w -f candidates.list > mapped_candidatesDH.list
grep -w -v -f mapped_candidatesDH.list candidates.list|wc -l

#Uniq
cut -f1 blast_pliv_dhuerta.out|sort|uniq -c |grep ' 1 ' |grep -w -f candidates.list | wc -l

#Repeated

cut -f1 blast_pliv_dhuerta.out|sort|uniq -c |grep -v ' 1 ' |grep -w -f candidates.list | wc -l

#ClassifyBlastout: this script can be found in: https://github.com/EvolutionaryGenetics-UB-CEAB/restrictionEnzimes.git 
source ~/Documents/gffutils/bin/activate
python3 ~/Documents/pliv/classifyBlastOut.py -b ~/Documents/pliv/ncbi_dataset/data/GCA_940671915.1/blast_pliv_dhuerta.out.txt -g ~/Documents/pliv/ncbi_dataset/data/GCA_940671915.1/plividus_v2.gff3 -o ~/Documents/pliv/ncbi_dataset/data/GCA_940671915.1/PL_huerta.out

#Recounting
#mapped y unmapped
 ##TOTAL
cut -f1 PL_huerta.out |sort| uniq|wc -l


 ##CANDIDATES
cut -f1 PL_huerta.out |sort| uniq| grep -w -f candidates.list |wc -l 


#Uniqs and repeated
 ##TOTAL
cut -f1 PL_huerta.out|sort|uniq -c |grep ' 1 ' |wc -l 
cut -f1 PL_huerta.out|sort|uniq -c |grep -v ' 1 ' |wc -l 
 #extraemos lista de unicos
cut -f1 PL_huerta.out |sort|uniq -c | grep ' 1 ' |sed 's/      1 //' > uniqPL_huerta.list

 ##CANDIDATES
cut -f1 PL_huerta.out|sort|uniq -c |grep ' 1 ' |grep -w -f candidates.list | wc -l 
cut -f1 PL_huerta.out|sort|uniq -c |grep -v ' 1 ' |grep -w -f candidates.list | wc -l 

#Gene-intergenic
 ##TOTAL
#depeding on the reference genome annotation you might not be able to search 'gene', however, you can searc intergenic, and the gene will be the summing loci mapping in exonic and intronic regions
awk 'FNR==NR{a[$1]; next} $1 in a && /intergenic/' uniqPL_huerta.list PL_huerta.out | wc -l 

 ##CANDIDATES
awk 'FNR==NR{a[$1]; next} $1 in a && /intergenic/' uniqPL_huerta.list PL_huerta.out | grep -w -f candidates.list |wc -l 

#Exon-intron
 ##TOTAl
awk 'FNR==NR{a[$1]; next} $1 in a && /exon/' uniqPL_huerta.list PL_huerta.out | wc -l 
awk 'FNR==NR{a[$1]; next} $1 in a && /intron/' uniqPL_huerta.list PL_huerta.out | wc -l
awk 'FNR==NR{a[$1]; next} $1 in a && /exon/ && /intron/' uniqPL_huerta.list PL_huerta.out | wc -l 
#due to the script used for classify them, we could have a loci classified as exon and intron, therefore, you will be counting them twice. So, in this specific case we consider the ones mapping in introns and exons at the same time as exons. 

 ##CANDIDATES
awk 'FNR==NR{a[$1]; next} $1 in a && /exon/' uniqPL_huerta.list PL_huerta.out | grep -w -f candidates.list| wc -l
awk 'FNR==NR{a[$1]; next} $1 in a && /intron/' uniqPL_huerta.list PL_huerta.out | grep -w -f candidates.list| wc -l
awk 'FNR==NR{a[$1]; next} $1 in a && /exon/ && /intron/' uniqPL_huerta.list PL_huerta.out |grep -w -f candidates.list|wc -l



