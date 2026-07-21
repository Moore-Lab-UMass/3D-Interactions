#Jill Moore
#UMass Chan
#July 2026


workingDir=~/Moore-Lab/Papers/3D-Genome

list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
ccreDir=$registryDir/Cell-Type-Specific/Individual-Files
dataDir=/data/projects/encode/data/

cd $workingDir/1_General-Comparisons
source ~/.bashrc

rm -f tmp.out running.output
for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    grep -v Low-DNase $ccreDir/$dnase*bed > tmp.cts

    rnap2=$(awk '{if (NR == '$i') print $2}' $list)
    ctcf=$(awk '{if (NR == '$i') print $4}' $list)
    hic=$(awk -F "\t" '{if (NR == '$i') print $19}' $list)

    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    awk '{if ($2 > 1.64) print $0}' $registryDir/signal-output/$ctcfChIP.txt | awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.ctcf


    assays=(RNAPII CTCF HiC)
    for assay in ${assays[@]}
    do
        links=Assay-Overlap-Details/$cell.$assay.txt
        totalLinks=$(wc -l $links | awk '{print $1}')

        grep PLS tmp.cts > tmp.ccre
	pls=$(bedtools pairtobed -type either -a $links -b tmp.ccre | awk '{print $1$2$3$4$5$6}' | \
	    sort -u | wc -l | awk '{print $1}')
	bedtools pairtobed -type neither -a $links -b tmp.ccre > tmp.next1

	grep ELS tmp.cts > tmp.ccre
	els=$(bedtools pairtobed -type either -a tmp.next1 -b tmp.ccre | awk '{print $1$2$3$4$5$6}' | \
	    sort -u | wc -l | awk '{print $1}')
        bedtools pairtobed -type neither -a tmp.next1 -b tmp.ccre > tmp.next2

	grep CTCF tmp.cts | cat tmp.ctcf - | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > tmp.ccre
        ctcf=$(bedtools pairtobed -type either -a tmp.next2 -b tmp.ccre | awk '{print $1$2$3$4$5$6}' | \
	    sort -u | wc -l | awk '{print $1}')
        bedtools pairtobed -type neither -a tmp.next2 -b tmp.ccre > tmp.next1

	other=$(bedtools pairtobed -type either -a tmp.next1 -b tmp.cts | \
	    awk '{print $1$2$3$4$5$6}' | sort -u | wc -l | awk '{print $1}')
        bedtools pairtobed -type neither -a tmp.next1 -b tmp.cts > tmp.next2

        inactive=$(bedtools pairtobed -type either -a tmp.next2 -b $registryDir/hg38-cCREs-Unfiltered.bed | \
            awk '{print $1$2$3$4$5$6}' | sort -u | wc -l | awk '{print $1}')

	none=$(bedtools pairtobed -type neither -a tmp.next2 -b $registryDir/hg38-cCREs-Unfiltered.bed | \
	    awk '{print $1$2$3$4$5$6}' | sort -u | wc -l | awk '{print $1}')

        echo -e $cell "\t" $assay "\t" $pls "\t" $els "\t" $ctcf "\t" $other "\t" $inactive "\t" $none >> running.output
    done
done
mv running.output Summary-cCRE-Overlap.txt
