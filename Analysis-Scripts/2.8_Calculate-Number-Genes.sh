
workingDir=~/Moore-Lab/Papers/3D-Genome/
list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=~/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38
ccreDir=$registryDir/Cell-Type-Specific/Individual-Files
dataDir=/data/projects/encode/data/
tss=~/Lab/Reference/Human/hg38/GENCODE40/TSS.Basic.bed
genes=~/Lab/Reference/Human/hg38/GENCODE40/Genes.Basic.bed
scriptDir=~/GitHub/3D-Interactions/scratch/2

cd $workingDir/2_Promoter-Centric

rm -f tmp.out
rm -f running.tmp
for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)

    echo $cell 
    grep PLS $ccreDir/$dnase*bed > tmp.pls
    grep dELS $ccreDir/$dnase*bed > tmp.dels

    grep CTCF-Anchors Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | \
	awk '{print $1}' | sort -u > a
    grep -v CTCF-Anchors Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | \
	awk '{print $1}' | sort -u > b
    awk 'FNR==NR {x[$0];next} ($0 in x)' a b  > tmp.arch
    awk 'FNR==NR {x[$0];next} !($0 in x)' a b  > tmp.noarch
    awk 'FNR==NR {x[$0];next} !($0 in x)' b a  >> tmp.arch
    
    groups=(arch noarch)
    for group in ${groups[@]}
    do
        echo $group

        awk 'FNR==NR {x[$1];next} ($1 in x)' tmp.$group Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | \
            awk '{print $1 "\t" $4}' | sort -u > tmp.pairs
        python ~/bin/count.py tmp.pairs 0 | awk '{print $0 "\t" "'$group'" "\t" "P-C"}' >> running.tmp
    done

    grep CTCF-Anchors Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt | grep P-Ed | \
	awk '{print $1}' | sort -u > a
    grep -v CTCF-Anchors Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt | grep P-Ed | \
	awk '{print $1}' | sort -u > b
    awk 'FNR==NR {x[$0];next} ($0 in x)' a b  > tmp.arch
    awk 'FNR==NR {x[$0];next} !($0 in x)' a b  > tmp.noarch
    awk 'FNR==NR {x[$0];next} !($0 in x)' b a  >> tmp.arch

    groups=(arch noarch)
    for group in ${groups[@]}
    do
        echo $group

        awk 'FNR==NR {x[$1];next} ($1 in x)' tmp.$group Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt | awk '{print $1 "\t" $4}' | sort -u > tmp.pairs
        python ~/bin/count.py tmp.pairs 0 | awk '{print $0 "\t" "'$group'" "\t" "P-E"}' >> running.tmp
    done

    mkdir -p Gene-Counts
    mv running.tmp Gene-Counts/$cell.PE-PC-Gene-Counts.txt
done

