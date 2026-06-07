
workingDir=~/Moore-Lab/Papers/3D-Genome
scriptDir=/home/moorej3/GitHub/3D-Interactions/Analysis-Scripts/Toolkit
list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
ccreDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/Cell-Type-Specific/Individual-Files
chiapetDir=/data/projects/encode/3D-Chromatin/ChIA-PET/Ruan-ENCODE4/
dataDir=/data/projects/encode/data/
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38

source ~/.bashrc

cd $workingDir/3_Dual-Elements
mkdir -p ChromBPNet/TF-Overlap

for i in `seq 3 1 4`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    chrombpnet=$(awk '{if (NR == '$i') print $20}' $list)
    bigWig=$workingDir/3_Dual-Elements/ChromBPNet/Attribution-BigWigs/$chrombpnet.bigWig
    tfSummits=$workingDir/0_Reference-Set-Up/All-$cell-Peak-Summits.bed

    echo $cell

    awk 'FNR==NR {x[$1];next} ($1 in x)' Candidate-Dual-Elements.txt \
	$workingDir/2_Promoter-Centric/Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt | \
	awk 'FNR==NR {x[$1];next} ($5 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.e

    awk 'FNR==NR {x[$1];next} ($1 in x)' Candidate-Dual-Elements.txt \
        $workingDir/2_Promoter-Centric/Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | \
        awk 'FNR==NR {x[$1];next} ($5 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.c

    bedtools intersect -c -a tmp.e -b $tfSummits | \
	awk '{print $4 "\t" $NF "\t" "enhancer-state"}' > tmp.output
    bedtools intersect -c -a tmp.c -b $tfSummits | \
        awk '{print $4 "\t" $NF "\t" "ctcf-only-state"}' >> tmp.output

    mv tmp.output ChromBPNet/TF-Overlap/$cell-TF-Summit-Overlap.txt

done

rm tmp.*



