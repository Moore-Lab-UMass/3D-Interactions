
workingDir=~/Moore-Lab/Papers/3D-Genome
scriptDir=/home/moorej3/GitHub/3D-Interactions/Analysis-Scripts/Toolkit
list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
ccreDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/Cell-Type-Specific/Individual-Files
chiapetDir=/data/projects/encode/3D-Chromatin/ChIA-PET/Ruan-ENCODE4/
dataDir=/data/projects/encode/data/
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38

source ~/.bashrc

cd $workingDir/3_Dual-Elements
mkdir -p CAPRA-Scores

cell=K562
capra=$workingDir/0_Reference-Set-Up/CAPRA-ENCSR661FOW.Solo-Filtered.txt
dnase=$(grep $cell $list | awk '{print $9}')

grep ELS $ccreDir/$dnase*bed > tmp.els


awk 'FNR==NR {x[$1];next} ($1 in x)' Candidate-Dual-Elements.txt \
    $workingDir/2_Promoter-Centric/Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt | \
    awk 'FNR==NR {x[$1];next} ($5 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.e

awk 'FNR==NR {x[$1];next} ($1 in x)' Candidate-Dual-Elements.txt \
    $workingDir/2_Promoter-Centric/Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | \
    awk 'FNR==NR {x[$1];next} ($5 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.c

awk 'FNR==NR {x[$4];next} ($1 in x)' tmp.e $capra | \
    awk '{print $1 "\t" $3 "\t" "enhancer-state"}' > tmp.output
awk 'FNR==NR {x[$4];next} ($1 in x)' tmp.c $capra | \
    awk '{print $1 "\t" $3 "\t" "ctcf-only-state"}' >> tmp.output

awk 'FNR==NR {x[$4];next} ($5 in x)' tmp.els $registryDir/hg38-cCREs-Unfiltered.bed | \
    awk 'FNR==NR {x[$4];next} ($1 in x)' - $capra | \
    awk '{print $1 "\t" $3 "\t" "K562-enhancer"}' >> tmp.output
sort -R $registryDir/hg38-cCREs-Unfiltered.bed | head -n 100000 | \
    awk 'FNR==NR {x[$4];next} ($1 in x)' - $capra | \
    awk '{print $1 "\t" $3 "\t" "background"}' >> tmp.output

mv tmp.output CAPRA-Scores/K562-CAPRA-Scores.txt

rm tmp.*



