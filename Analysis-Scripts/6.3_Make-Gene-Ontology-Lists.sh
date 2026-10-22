#Jill Moore
#UMass Chan
#October 2026

scriptDir=/data/zusers/moorej3/GitHub/3D-Interactions/Analysis-Scripts/Toolkit
control=/data/projects/encode/Registry/V4/GRCh38/GRCh38-cCREs-All.bed

workingDir=/data/zusers/moorej3/Moore-Lab/Papers/3D-Genome
list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
dataDir=/data/projects/encode/data/

cd $workingDir/6_Biological-Spectrum/Gene-Ontology

for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    echo $cell

    cat $workingDir/2_Promoter-Centric/Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt \
	$workingDir/2_Promoter-Centric/Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | \
	grep CTCF-Anchors | awk '{print $NF}' | sort -u | awk -F "." '{print $1}' > $cell.Arch-Genes.txt

    cat $workingDir/2_Promoter-Centric/Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt \
        $workingDir/2_Promoter-Centric/Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | \
        awk '{print $NF}' | sort -u | awk -F "." '{print $1}' > tmp.x

    cat $workingDir/2_Promoter-Centric/Promoter-Promoter-Pairs/$cell.P-P-Pairs.txt | \
	awk '{print $4 "\n" $5}' | awk -F "." '{print $1}' >> tmp.x

    sort -u tmp.x > $cell.Background-Genes.txt

done
rm tmp.x
