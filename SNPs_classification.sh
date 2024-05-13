##Classifying SNPs mapping into the reference genome##
##Author: Ainhoa Lopez
##Date: 22 April 2024

#We need to have a gtf/gff file from the reference genome and the bed file with our SNPs
/home/soft/bedtools2-2.30.0/bedtools intersect -a /your/file/path/plividus_annotationDEF_NCBI_v2.gtf -b /your/file/path/SNPs_NCBI.bed |grep -w 'gene' > snps_genes_NCBI_v2.out #to obtain the ones located in genes
/home/soft/bedtools2-2.30.0/bedtools intersect -a /your/file/path/plividus_annotationDEF_NCBI_v2.gtf -b /your/file/path/SNPs_NCBI.bed |grep -w 'exon' > snps_exons_NCBI_v2.out #to obtain the ones located in exons

cut -f1,4,5 snps_genes_NCBI_v2.out > genes_coord_v2.txt
cut -f1,4,5 snps_exons_NCBI_v2.out > exones_coord_v2.txt
cut -f1,4,5 snps_exons_NCBI_v2.out >snps_exons_v2.bed
/home/soft/bedtools2-2.30.0/bedtools intersect -a snps_genes_NCBI_v2.out -b snps_exons_v2.bed -v > snps_introns_NCBI_v2.out 

cut -f1,4,5 snps_introns_NCBI_v2.out  > intrones_coord_v2.txt

#How many SNPs are located in genic regions:
cut -f1,4,5 snps_genes_NCBI_v2.out|sort|uniq|wc -l #totals: 84587

sort genes_coord_v2.txt|uniq -c|grep ' 1 '|wc -l  #60333 uniqs, meaning the ones that map to only one location
sort genes_coord_v2.txt|uniq -c|grep -v ' 1 '|wc -l #24254 repeated, meaning the ones that map to more than one location 

#We generate the uniq list file
sort genes_coord_v2.txt|uniq -c|grep ' 1 '|sed 's/      1 //' > llista_SNPs_UNIQs_in_genes_v2.list

#We classify for the SNPs only mapping into uniq positions
#In this specific case we have two type: StringTie (long-non-coding protein genes) and -v StringTie (protein-coding genes)
grep 'StringTie' snps_genes_NCBI_v2.out |cut -f1,4,5 > genes_coord_StrinTie.txt
grep -v 'StringTie' snps_genes_NCBI_v2.out |cut -f1,4,5 > genes_coord_Augustus.txt


awk 'FNR==NR {array[$1,$2]; next} ($1,$2) in array' llista_SNPs_UNIQs_in_genes_v2.list genes_coord_Augustus.txt |wc -l #to check how many located in protein coding genes
awk 'FNR==NR {array[$1,$2]; next} ($1,$2) in array' llista_SNPs_UNIQs_in_genes_v2.list genes_coord_StrinTie.txt |wc -l #to check how many located in long non coding protein genes
awk 'NR==FNR { snps[$0] = 1; next } $0 in snps' llista_SNPs_UNIQs_in_genes_v2.list exones_coord_v2.txt|sort|uniq|wc -l #to check how many are in exonic regions
awk 'NR==FNR { snps[$0] = 1; next } $0 in snps' llista_SNPs_UNIQs_in_genes.list intrones_coord_v2.txt|sort|uniq|wc -l #to check how many are in intronic regions

#The same process is performed for SNPs under selection detected by RDA
/home/soft/bedtools2-2.30.0/bedtools intersect -a /your/file/path/plividus_annotationDEF_NCBI_v2.gtf -b /your/file/path/SNPs_candidats.bed |grep -w 'gene' > snps_candidates_genes_NCBI.out #to obtain the ones located in genes
/home/soft/bedtools2-2.30.0/bedtools intersect -a /your/file/path/plividus_annotationDEF_NCBI_v2.gtf -b /your/file/path/SNPs_candidats.bed |grep -w 'exon' > snps_candidates_exons_NCBI.out #to obtain the ones located in exons

cut -f1,4,5 snps_candidates_genes_NCBI.out > cand_genes_coord.txt
cut -f1,4,5 snps_candidates_exons_NCBI.out > cand_exons_coord.txt
cut -f1,4,5 snps_candidates_exons_NCBI.out > cand_exons.bed
/home/soft/bedtools2-2.30.0/bedtools intersect -a snps_candidates_genes_NCBI.out -b cand_exons.bed -v > cand_introns.out 

cut -f1,4,5 cand_introns.out > cand_introns_cord.txt

#How many SNPs are located in genic regions:
cut -f1,4,5 snps_candidates_genes_NCBI.out|sort|uniq|wc -l #total

sort cand_genes_coord.txt|uniq -c|grep ' 1 '|wc -l  #uniq
sort cand_genes_coord.txt|uniq -c|grep -v ' 1 '|wc -l #repeated 

#SNPs_uniqs list
sort cand_genes_coord.txt|uniq -c|grep ' 1 '|sed 's/      1 //' > cand_SNPs_UNIQs_in_genes_v2.list

#Classify uniq SNPs
grep 'StringTie' snps_candidates_genes_NCBI.out |cut -f1,4,5 > cand_coord_StrinTie.txt
grep -v 'StringTie' snps_candidates_genes_NCBI.out |cut -f1,4,5 > cand_coord_Augustus.txt


awk 'FNR==NR {array[$1,$2]; next} ($1,$2) in array' cand_SNPs_UNIQs_in_genes_v2.list cand_coord_Augustus.txt |wc -l #SNPs mapping to protein coding genes
awk 'FNR==NR {array[$1,$2]; next} ($1,$2) in array' cand_SNPs_UNIQs_in_genes_v2.list cand_coord_StrinTie.txt |wc -l #SNPs mapping to long non protein coding genes
awk 'NR==FNR { snps[$0] = 1; next } $0 in snps' cand_SNPs_UNIQs_in_genes_v2.list cand_exons_coord.txt|sort|uniq|wc -l #mapping in exons
awk 'NR==FNR { snps[$0] = 1; next } $0 in snps' cand_SNPs_UNIQs_in_genes_v2.list cand_introns_cord.txt|sort|uniq|wc -l #mapping in introns
