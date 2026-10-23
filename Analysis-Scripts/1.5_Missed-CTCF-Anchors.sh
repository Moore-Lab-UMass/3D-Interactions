#Jill Moore
#UMass Chan
#October 2026


workingDir=~/Moore-Lab/Papers/3D-Genome

list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
ccreDir=$registryDir/Cell-Type-Specific/Individual-Files
dataDir=/data/projects/encode/data/
scriptDir=~/GitHub/3D-Interactions/Analysis-Scripts/Toolkit

cd $workingDir/2_Promoter-Centric

for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    #echo $cell

    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    awk '{if ($2 > 2.32) print $0}' $registryDir/signal-output/$ctcfChIP.txt |
        awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.ctcf

    links=Promoter-Centric-Annotations/$cell.Promoter-Links.txt

    assays=(RNAPII CTCF HiC)
    for assay in ${assays[@]}
    do
	#echo $assay
	grep "__"$assay $links | grep -v "CTCF-Anchor" | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' | bedtools pairtobed -type either -a - -b tmp.ctcf | \
	    awk '{print $1 "\t" $2 "\t" $3 "\n" $4 "\t" $5 "\t" $6}' | bedtools intersect -v -a stdin -b tmp.ctcf | sort -k1,1 -k2,2n | \
	    bedtools closest -d -a stdin -b tmp.ctcf  >> tmp.out
    done
    cat tmp.out | awk '{if ($NF > 2000) sum += 1}END{print "'$cell'" "\t" NR "\t" sum "\t" sum/NR}' >> running.out
    rm tmp.out
done

mv running.out Summary-Missed-CTCF-Anchors.txt
