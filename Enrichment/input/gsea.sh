#!/bin/bash

# $1 is name of analysis 
# $2 is full path for rnk file

# Home Directory variable
DIR=/data/Jaryd/R/CellSummerProject/Enrichment/input

#This loop will only iterate over data for Coefs, need to add another loop once name for other data types for a given pheno has been decided.

declare -a phenos=("CS" "MAD" "Props" "Fraction") # "p1s1")
declare -a metrics=("Coef" "Psmall" "Pbig" "CoefRank" "PSmallRank" "PBigRank")

for metric in "${metrics[@]}"
do
for pheno in "${phenos[@]}"
do
for i in {1..7}
do
if [ $pheno -eq "Props" ] && [ $i -eq 7 ];
then
break
fi
name=$pheno"_"$i"_"$metric
path=$DIR"/"$phenos"_"$i"_"$metric".rnk"

# command to run gsea Preranked with my specifications
java -cp /home/jaryd/gsea_home/gsea-3.0.jar -Xmx512m xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c2.cp.kegg.v6.1.symbols.gmt -norm meandiv -nperm 1000 -scoring_scheme weighted -rpt_label $name -create_svgs false -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out /data/Jaryd/R/CellSummerProject/Enrichment/results -gui false -rnk $path 
done
done
done
