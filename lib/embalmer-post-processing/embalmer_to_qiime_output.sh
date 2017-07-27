# this script takes an embalmer output file, generates an OTU file, filters then generates all QIIME output files necessary for analyses
# note that embalmer OTU table will NOT contain taxonomy info

# usage:
# bash *.sh otutable.txt tree.txt outputfolder
FP=$1 
TREE=$2
outfolder=$3

# full pipeline
 
# convert to JSON BIOM file
echo "biom convert -i ${FP} --table-type=\"OTU table\" --to-json -o ${outfolder}/otu0.biom"
biom convert -i ${FP} --table-type "OTU table" --to-json -o ${outfolder}/otu0.biom
# for OTU tables with taxonomy info, use:
# echo "biom convert -i ${FP} --process-obs-metadata taxonomy --table-type=\"OTU table\" --to-json -o ${outfolder}/otu0.biom"
# biom convert -i ${FP} --process-obs-metadata taxonomy --table-type "OTU table" --to-json -o ${outfolder}/otu0.biom

# remove OTU present in only one sample
echo "filter_otus_from_otu_table.py -i ${outfolder}/otu0.biom -o ${outfolder}/otu0_s2.biom -s 2"
filter_otus_from_otu_table.py -i ${outfolder}/otu0.biom -o ${outfolder}/otu0_s2.biom -s 2

# get stats
echo "biom summarize-table -i ${outfolder}/otu0_s2.biom -o ${outfolder}/stats.txt"
biom summarize-table -i ${outfolder}/otu0_s2.biom -o ${outfolder}/stats.txt

read -p "Enter filter sample depth (check stats.txt):" depth

# filter low depth samples from OTU table
echo "filter_samples_from_otu_table.py -i ${outfolder}/otu0_s2.biom -o ${outfolder}/otu0_s2_f.biom -n ${depth}"
filter_samples_from_otu_table.py -i ${outfolder}/otu0_s2.biom -o ${outfolder}/otu0_s2_f.biom -n ${depth}

echo "biom convert -i ${outfolder}/otu0_s2_f.biom -o ${outfolder}/otu0_s2_f.txt --to-tsv"
biom convert -i ${outfolder}/otu0_s2_f.biom -o ${outfolder}/otu0_s2_f.txt --to-tsv
# for OTU tables with taxonomy info, use:
# echo "biom convert -i ${outfolder}/otu0_s2_f.biom -o ${outfolder}/otu0_s2_f.txt --to-tsv --header-key=taxonomy"
# biom convert -i ${outfolder}/otu0_s2_f.biom -o ${outfolder}/otu0_s2_f.txt --to-tsv --header-key=taxonomy

#This will convert the OTU table to relative abundance, multiply all relative abundances by X (based on smallest number in OTU table) in order to get all abundances >1 and converted to a whole number (mimicking abs counts).
echo "Rscript /Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/convert.to.relative.abundance.r -i ${outfolder}/otu0_s2_f.txt -o ${outfolder}/final_otu.txt"
Rscript /Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/convert.to.relative.abundance.r -i ${outfolder}/otu0_s2_f.txt -o ${outfolder}/final_otu.txt

# convert back to biom
echo "biom convert -i ${outfolder}/final_otu.txt -o ${outfolder}/final_otu.biom --table-type=\"OTU table\" --to-json"
biom convert -i ${outfolder}/final_otu.txt -o ${outfolder}/final_otu.biom --table-type="OTU table" --to-json
# for OTU tables with taxonomy info, use:
# echo "biom convert -i ${outfolder}/final_otu.txt -o ${outfolder}/final_otu.biom --table-type=\"OTU table\" --process-obs-metadata=\"taxonomy\" --to-json"
# biom convert -i ${outfolder}/final_otu.txt -o ${outfolder}/final_otu.biom --table-type="OTU table" --process-obs-metadata="taxonomy" --to-json

# skip this for embalmer output
#get species-, genus-, and phylum-level taxonomy
#echo "summarize_taxa.py -i ${outfolder}/final_otu.biom -o ./taxa -L 2,6,7"
#summarize_taxa.py -i ${outfolder}/final_otu.biom -o ./taxa -L 2,6,7

#generate beta diversity
echo "beta_diversity.py -i ${outfolder}/final_otu.biom -o ${outfolder}/beta -m bray_curtis,unweighted_unifrac,weighted_unifrac -t ${TREE}"
beta_diversity.py -i ${outfolder}/final_otu.biom -o ${outfolder}/beta -m bray_curtis,unweighted_unifrac,weighted_unifrac -t ${TREE}

#echo "beta_diversity_through_plots.py -i ${outfolder}/final_otu.biom -m ${outfolder}/mapping.txt -o ${outfolder}/beta -t /Users/pvangay/bin/gg_13_8_otus/trees/97_otus.tree"
#beta_diversity_through_plots.py -i ${outfolder}/final_otu.biom -m ${outfolder}/mapping.txt -o ${outfolder}/beta -t /Users/pvangay/bin/gg_13_8_otus/trees/97_otus.tree 

#generate alpha diversity 
echo "alpha_diversity.py -i ${outfolder}/final_otu.biom -o ${outfolder}/alpha.txt -t ${TREE} -m PD_whole_tree,chao1,observed_otus,shannon,simpson"
alpha_diversity.py -i ${outfolder}/final_otu.biom -o ${outfolder}/alpha.txt -t ${TREE} -m PD_whole_tree,chao1,observed_otus,shannon,simpson
# echo "alpha_diversity.py -i ${outfolder}/final_otu.biom -o ${outfolder}/alpha.txt -t /Users/pvangay/bin/gg_13_8_otus/trees/97_otus.tree"
# alpha_diversity.py -i ${outfolder}/final_otu.biom -o ${outfolder}/alpha.txt -t /Users/pvangay/bin/gg_13_8_otus/trees/97_otus.tree


# let's rename a bunch of the files 
mv ${outfolder}/beta/bray_curtis_final_otu.txt ${outfolder}/beta/bray_curtis_dm.txt
mv ${outfolder}/beta/weighted_unifrac_final_otu.txt ${outfolder}/beta/weighted_unifrac_dm.txt
mv ${outfolder}/beta/unweighted_unifrac_final_otu.txt ${outfolder}/beta/unweighted_unifrac_dm.txt
