
workingDir=~/Moore-Lab/Papers/3D-Genome
scriptDir=/home/moorej3/GitHub/3D-Interactions/Analysis-Scripts/Toolkit
list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
ccreDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/Cell-Type-Specific/Individual-Files
chiapetDir=/data/projects/encode/3D-Chromatin/ChIA-PET/Ruan-ENCODE4/
dataDir=/data/projects/encode/data/
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38
tfBed=$workingDir/0_Reference-Set-Up/Representative-TF-rPeak-Motifs.bed
tfList=$workingDir/0_Reference-Set-Up/TF-100-List.txt

regions=$workingDir/3_Dual-Elements/Candidate-Dual-Elements.txt

cd $workingDir/4_TF-Conservation

for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    chrombpnet=$(awk '{if (NR == '$i') print $20}' $list)
    bigWig=$workingDir/3_Dual-Elements/ChromBPNet/Attribution-BigWigs/$chrombpnet.bigWig

    echo $cell

    awk 'FNR==NR {x[$1];next} ($1 in x)' $regions \
	$workingDir/2_Promoter-Centric/Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt | \
	awk 'FNR==NR {x[$1];next} ($5 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.e

    awk 'FNR==NR {x[$1];next} ($1 in x)' $regions \
        $workingDir/2_Promoter-Centric/Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | \
        awk 'FNR==NR {x[$1];next} ($5 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.c


    rm -f running.output
    for j in `seq 1 1 42`
    do
        echo $j
        tf=$(awk -F "\t" '{if (NR == '$j') print $1}' $tfList)

        awk -F "\t" '{if ($4 ~ /'$tf'_/) print $0}' $tfBed | \
            bedtools intersect -u -b tmp.e -a stdin | \
            awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > tmp.tf-e
	awk -F "\t" '{if ($4 ~ /'$tf'_/) print $0}' $tfBed | \
            bedtools intersect -u -b tmp.c -a stdin | \
            awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > tmp.tf-c

        ~/bin/bigWigAverageOverBed $bigWig tmp.tf-e out.tab
        E=$(sort -k5,5g out.tab | \
            awk '{a[i++]=$5;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')
        awk '{print $1 "\t" $5 "\t" "'$tf'" "\t" "Enhancer"}' out.tab >> running.output
	awk '{print $1 "\t" $5 "\t" "'$tf'" "\t" "Enhancer"}' out.tab > tmp.calc

	~/bin/bigWigAverageOverBed $bigWig tmp.tf-c out.tab
	C=$(sort -k5,5g out.tab | \
            awk '{a[i++]=$5;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')
        awk '{print $1 "\t" $5 "\t" "'$tf'" "\t" "CTCF"}' out.tab >> running.output
	awk '{print $1 "\t" $5 "\t" "'$tf'" "\t" "CTCF"}' out.tab >> tmp.calc
        p=$(Rscript $scriptDir/wilcox-test-chrombpnet.R | awk '{print $1}')

        echo -e $cell "\t" $tf "\t" $E "\t" $C "\t" $p >> running.summary
	#echo -e $cell "\t" $tf "\t" $E "\t" $C

    done
    mv running.output Motif-Attribution/Motif-Attribution-$cell.txt

done
mv running.summary Motif-Attribution/Summary-Motif-Attribution.txt
rm tmp.*

