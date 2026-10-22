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

cd $workingDir/6_Biological-Spectrum

source ~/.bashrc

awk '{print $1}' $workingDir/2_Promoter-Centric/Promoter-Enhancer-Pairs/*.P-E-Pairs.txt > tmp.pe
awk '{print $1 "\n" $2}' $workingDir/2_Promoter-Centric/Promoter-Promoter-Pairs/*.P-P-Pairs.txt >> tmp.pe

for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')

    grep CTCF-Anchor $workingDir/2_Promoter-Centric/Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt > tmp.ctcf
    total=$(wc -l tmp.ctcf | awk '{print $1}')
    awk 'FNR==NR {x[$1];next} !($1 in x)' tmp.pe tmp.ctcf | wc -l |
	awk '{print "'$cell'" "\t" $1 "\t" "'$total'" "\t" $1/'$total'}' >> running.results
done
mv running.results Summary-Remaining-CTCF-Bound.txt

rm tmp.*
