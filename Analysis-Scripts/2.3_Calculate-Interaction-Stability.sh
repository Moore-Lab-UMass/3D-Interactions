
workingDir=~/Moore-Lab/Papers/3D-Genome

list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
ccreDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/Cell-Type-Specific/Individual-Files
chiapetDir=/data/projects/encode/3D-Chromatin/ChIA-PET/Ruan-ENCODE4/
dataDir=/data/projects/encode/data/
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38

cd $workingDir/3_Dual-Elements

rm -f running.*
for i in `seq 1 1 6`
do
    cell1=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    echo $cell1

    links1=$workingDir/2_Promoter-Centric/Promoter-Centric-Annotations/$cell1.Promoter-Links.txt
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7"__"$8"__"$9}' $links1 > tmp.links
    grep P-P tmp.links > tmp.links1
    grep P-E tmp.links | awk '{gsub(/P-Ed/,"P-E");print}' | \
	awk '{gsub(/P-Ep/,"P-E");print}' >> tmp.links1
    grep P-C tmp.links >> tmp.links1

    rm -f running.overlap

    for j in `seq 1 1 6`
    do
        cell2=$(awk '{if (NR == '$j') print $1}' $list | awk -F "_" '{print $3}')
        links2=$workingDir/2_Promoter-Centric/Promoter-Centric-Annotations/$cell2.Promoter-Links.txt
        awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7"__"$8"__"$9}' $links2 > tmp.links2

        bedtools pairtopair -is -type both -a tmp.links1 -b tmp.links2 | \
	    awk '{print $7}' | sort -u >> running.overlap

    done
    python ~/bin/count.py running.overlap 0 > tmp.counts

    groups1=(P-P P-E P-C)
    groups2=(CTCF-Anchor No-CTCF)

    for group1 in ${groups1[@]}
    do
        for group2 in ${groups2[@]}
	do
	    grep $group1 tmp.counts | grep $group2 | \
		awk '{if ($NF > 3) high +=1}END{print "'$cell1'" "\t" "'$group1'" "\t" "'$group2'" \
		"\t" NR "\t" high "\t" high/NR}' >> running.summary
	done
    done

    mkdir -p Stability-Cross-Lines
    mv tmp.counts Stability-Cross-Lines/$cell1.Stabilty-Counts.txt
    rm tmp.*
done


mv running.summary Stability-Cross-Lines/Summary-Stability-Counts.txt
rm running.*
