#Jill Moore
#UMass Chan
#October 2026

workingDir=~/Moore-Lab/Papers/3D-Genome/
list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=~/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38
ccreDir=$registryDir/Cell-Type-Specific/Individual-Files
dataDir=/data/projects/encode/data/
tss=~/Lab/Reference/Human/hg38/GENCODE40/TSS.Basic.bed
genes=~/Lab/Reference/Human/hg38/GENCODE40/Genes.Basic.bed
scriptDir=~/GitHub/3D-Interactions/scratch/2

cd $workingDir/2_Promoter-Centric

rm -f running.tmp
for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    groups=(P-P P-E P-C)
    for group in ${groups[@]}
    do
        grep $group Promoter-Centric-Annotations/$cell.Promoter-Links.txt > tmp.links
	a=$(grep CTCF-Anchors tmp.links | awk '{print $5-$3}' | sort -k1,1rg | awk '{a[i++]=$1;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')
	b=$(grep -v CTCF-Anchors tmp.links | awk '{print $5-$3}' | sort -k1,1rg | awk '{a[i++]=$1;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')

	echo -e $cell "\t" $group "\t" $a "\t" $b >> running.tmp
    done
done

mv running.tmp ../3_Candidate-Architectural/Summary-Median-Distances.txt
