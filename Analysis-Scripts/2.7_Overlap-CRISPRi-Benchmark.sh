#Jill Moore
#UMass Chan
#October 2026

workingDir=~/Moore-Lab/Papers/3D-Genome

list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
dataDir=/data/projects/encode/data/
ccreDir=$registryDir/Cell-Type-Specific/Individual-Files
#cd $workingDir/1_General-Comparisons

cd ~/Moore-Lab/Papers/3D-Genome/2_Promoter-Centric/CRISPRi/

#https://github.com/EngreitzLab/CRISPR_comparison/blob/v1.0.0/resources/crispr_data/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz 

awk '{if ($11 == "TRUE") print $0}' EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv | \
    awk '{if ($3 < $7) print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $7 "\t" $8 "\t" $4; \
    else print $6 "\t" $7 "\t" $8 "\t" $1 "\t" $2 "\t" $3 "\t" $4}' | 
    awk '{if ($5-$3 >= 10000) print $0}' > True-K562-10k.bedpe
awk '{if ($11 == "FALSE") print $0}' EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv | \
    awk '{if ($3 < $7) print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $7 "\t" $8 "\t" $4; \
    else print $6 "\t" $7 "\t" $8 "\t" $1 "\t" $2 "\t" $3 "\t" $4}' |
    awk '{if ($5-$3 >= 10000) print $0}' > False-K562-10k.bedpe

i=4 #K562
rm -f running.out1 running.out2
cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
echo $cell

rnap2=$(awk '{if (NR == '$i') print $2}' $list)
ctcf=$(awk '{if (NR == '$i') print $4}' $list)
hic=$(awk -F "\t" '{if (NR == '$i') print $19}' $list)


links=$workingDir/2_Promoter-Centric/Promoter-Centric-Annotations/$cell.Promoter-Links.txt

false=$(awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7"__"$8"__"$9}' $links | \
    bedtools pairtopair  -a False-K562-10k.bedpe -b stdin | awk '{print $7}' | sort -u | wc -l | \
    awk '{print $1}')
true=$(awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7"__"$8"__"$9}' $links | \
    bedtools pairtopair  -a True-K562-10k.bedpe -b stdin | awk '{print $7}' | sort -u | wc -l | \
    awk '{print $1}')

echo -e $cell "\t" "all" "\t" $true "\t" $false >> running.out

groups=(P-P P-E P-C)
for group in ${groups[@]}
do
    echo $group
    cat $links | grep $group > tmp.links
    false=$(awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7"__"$8"__"$9}' tmp.links | \
	bedtools pairtopair  -a False-K562-10k.bedpe -b stdin | awk '{print $7}' | sort -u | wc -l | \
	awk '{print $1}')
    true=$(awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7"__"$8"__"$9}' tmp.links | \
        bedtools pairtopair  -a True-K562-10k.bedpe -b stdin | awk '{print $7}' | sort -u | wc -l | \
        awk '{print $1}')
echo -e $cell "\t" $group "\t" $true "\t" $false >> running.out

done

mv running.out Summary-K562-CRISPR-Overlap.txt
rm -r -f tmp*
