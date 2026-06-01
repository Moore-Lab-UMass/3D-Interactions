


workingDir=~/Moore-Lab/Papers/3D-Genome

list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
dataDir=/data/projects/encode/data/

cd $workingDir/1_General-Comparisons


rm -f tmp.out
for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')

    rnap2=$(awk '{if (NR == '$i') print $2}' $list)
    ctcf=$(awk '{if (NR == '$i') print $4}' $list)
    hic=$(awk -F "\t" '{if (NR == '$i') print $19}' $list)

    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    awk '{if ($2 > 2.32) print $0}' $registryDir/signal-output/$ctcfChIP.txt | 
        awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.ctcf

    assays=(RNAPII CTCF HiC)
    for assay in ${assays[@]}
    do
        links=Assay-Overlap-Details/$cell.$assay.txt
	totalLinks=$(wc -l $links | awk '{print $1}')
	A=$(bedtools pairtobed -type either -a $links -b tmp.ctcf | awk '{print $1$2$3$4$5$6}' | \
	    sort -u | wc -l | awk '{print $1}')
	B=$(bedtools pairtobed -type both -a $links -b tmp.ctcf | awk '{print $1$2$3$4$5$6}' | \
	    sort -u | wc -l | awk '{print $1}')
        echo -e $cell "\t" $assay "\t" $A "\t" $B "\t" $totalLinks
    done
done

