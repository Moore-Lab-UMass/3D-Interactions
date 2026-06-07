
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
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    echo $cell
    grep P-C $workingDir/2_Promoter-Centric/Promoter-Centric-Annotations/$cell.Promoter-Links.txt | grep CTCF-Anchor | \
	awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7"__"$8"__"$9}' >> running.pc
done

for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)

    grep P-E $workingDir/2_Promoter-Centric/Promoter-Centric-Annotations/$cell.Promoter-Links.txt | \
	grep CTCF-Anchor | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7"__"$8"__"$9}' > tmp.pe-arch
    bedtools pairtopair -is -type both -a tmp.pe-arch -b running.pc | awk '{print $7}' > tmp.rep

    awk '{print $1 "\t" $3}' $workingDir/2_Promoter-Centric/Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt | sort -u > tmp.pairs
    python ~/bin/count.py tmp.pairs -1 | awk '{if ($NF == 1) print $0}' | \
	awk 'FNR==NR {x[$1];next} ($3 in x)' - $workingDir/2_Promoter-Centric/Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt | \
	grep CTCF-Anchor | awk 'FNR==NR {x[$1];next} ($3 in x)' tmp.rep - | \
	awk '{print $1}' | sort -u >> running.singleton-e
    rm tmp.*
done


sort -u running.singleton-e > Singleton-Enhancers.txt
cat ../2_Promoter-Centric/Promoter-CTCF-Pairs/*txt | \
    awk 'FNR==NR {x[$1];next} ($1 in x)' - Singleton-Enhancers.txt > Candidate-Dual-Elements.txt

rm running*

