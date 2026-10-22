
workingDir=~/Moore-Lab/Papers/3D-Genome
scriptDir=/home/moorej3/GitHub/3D-Interactions/Analysis-Scripts/Toolkit
list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
ccreDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/Cell-Type-Specific/Individual-Files
chiapetDir=/data/projects/encode/3D-Chromatin/ChIA-PET/Ruan-ENCODE4/
dataDir=/data/projects/encode/data/
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38

source ~/.bashrc
conda activate deeptools

cd $workingDir/3_Dual-Elements
mkdir -p ChromBPNet/Attribution-Calculations/

for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    chrombpnet=$(awk '{if (NR == '$i') print $20}' $list)
    bigWig=$workingDir/3_Dual-Elements/ChromBPNet/Attribution-BigWigs/$chrombpnet.bigWig

#pull-high-chrombpnet.py

    awk 'FNR==NR {x[$1];next} ($1 in x)' Candidate-Dual-Elements.txt \
	$workingDir/2_Promoter-Centric/Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt | \
	awk 'FNR==NR {x[$1];next} ($5 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.e

    awk 'FNR==NR {x[$1];next} ($1 in x)' Candidate-Dual-Elements.txt \
        $workingDir/2_Promoter-Centric/Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | \
        awk 'FNR==NR {x[$1];next} ($5 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.c
    
    python $scriptDir/pull-high-chrombpnet.py $bigWig tmp.e > out.1
    python $scriptDir/pull-high-chrombpnet.py $bigWig tmp.c > out.2

    ENorm=$(sort -k7,7rg out.1 | \
	awk '{a[i++]=$7/($3-$2);}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')
    CNorm=$(sort -k7,7rg out.2 | \
	awk '{a[i++]=$7/($3-$2);}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')

    ERaw=$(sort -k7,7rg out.1 | \
        awk '{a[i++]=$7;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')
    CRaw=$(sort -k7,7rg out.2 | \
        awk '{a[i++]=$7;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')
    echo -e $cell "\t" $ERaw "\t" $CRaw "\t" $ENorm "\t" $CNorm >> running.output

    awk '{print $0 "\t" "enhancer-state"}' out.1 > tmp.combo
    awk '{print $0 "\t" "ctcf-only-state"}' out.2 >> tmp.combo

    mv tmp.combo ChromBPNet/Attribution-Calculations/$cell-ChromBPNet-EC-Attributions.txt

done

mv running.output ChromBPNet/Summary-Attribution-Calculations.txt



