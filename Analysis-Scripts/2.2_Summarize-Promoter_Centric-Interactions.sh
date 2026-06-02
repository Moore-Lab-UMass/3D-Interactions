
workingDir=~/Moore-Lab/Papers/3D-Genome

list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
ccreDir=$registryDir/Cell-Type-Specific/Individual-Files
dataDir=/data/projects/encode/data/
scriptDir=~/GitHub/3D-Interactions/Analysis-Scripts/Toolkit

cd $workingDir/2_Promoter-Centric
mkdir -p Promoter-Centric-Annotations

rm -f Summary-Promoter-Connected-*Breakdown.txt
for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    echo $cell


    links=Promoter-Centric-Annotations/$cell.Promoter-Links.txt
    python $scriptDir/count-cCRE-groups.py  $links 7 | awk '{print $2}' | \
	tr '\n' '\t' | awk '{print "'$cell'" "\t" $0}' >> Summary-Promoter-Connected-cCRE-Breakdown.txt

    assays=(RNAPII CTCF HiC)
    for assay in ${assays[@]}
    do
	echo $assay
	grep "__"$assay $links > tmp.links
	python $scriptDir/count-cCRE-groups.py  tmp.links 7 | \
	    awk '{print $2}' | tr '\n' '\t' | \
	    awk '{print "'$cell'"  "\t" "'$assay'" "\t" $0}' >> Summary-Promoter-Connected-Assay-Breakdown.txt
    done


    modes=(CTCF-Anchors No-CTCF)
    for mode in ${modes[@]}
    do
        echo $mode
        grep $mode $links > tmp.links
        python $scriptDir/count-cCRE-groups.py  tmp.links 7 | \
            awk '{print $2}' | tr '\n' '\t' | \
            awk '{print "'$cell'"  "\t" "'$mode'" "\t" $0}' >> Summary-Promoter-Connected-Architecture-Breakdown.txt
    done

done
