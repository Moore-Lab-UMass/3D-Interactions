#Jill Moore
#UMass Chan
#October 2026

workingDir=~/Moore-Lab/Papers/3D-Genome

list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
ccreDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/Cell-Type-Specific/Individual-Files
chiapetDir=/data/projects/encode/3D-Chromatin/ChIA-PET/Ruan-ENCODE4/
dataDir=/data/projects/encode/data/
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38
scriptDir=~/GitHub/3D-Interactions/Analysis-Scripts/Toolkit

cd $workingDir/4_Dual-Elements

source ~/.bashrc

rm -f running.*

for i in `seq 1 1 6`
do
    cell1=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    grep P-E $workingDir/2_Promoter-Centric/Promoter-Centric-Annotations/$cell1.Promoter-Links.txt | \
	grep "CTCF-Anchor" | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7"__"$8"__"$9}' > tmp.1

    for i in `seq 1 1 6`
    do
        cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
        dnase=$(awk '{if (NR == '$i') print $9}' $list)
        ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)

        awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7"__"$8"__"$9}' $workingDir/2_Promoter-Centric/Promoter-Centric-Annotations/$cell.Promoter-Links.txt > tmp.2


        bedtools pairtopair -is -type both -a tmp.2 -b tmp.1 | awk '{print $7}' | \
	    sort -u | awk '{gsub(/__/,"\t");print}'  > tmp.out

        python $scriptDir/count-replicated-arch.py tmp.out | awk '{print "'$cell1'" "\t" "'$cell'" "\t" $0}' >> running.output

    done
    rm tmp.*
done
mv running.output Summary-Arc-P-E-Recurrance.txt
