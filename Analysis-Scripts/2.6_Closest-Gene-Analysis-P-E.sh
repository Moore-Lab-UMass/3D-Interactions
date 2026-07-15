
workingDir=~/Moore-Lab/Papers/3D-Genome/
list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=~/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38
ccreDir=$registryDir/Cell-Type-Specific/Individual-Files
dataDir=/data/projects/encode/data/
tss=~/Lab/Reference/Human/hg38/GENCODE40/TSS.Basic.bed
genes=~/Lab/Reference/Human/hg38/GENCODE40/Genes.Basic.bed
scriptDir=~/GitHub/3D-Interactions/scratch/2

cd $workingDir/2_Promoter-Centric

source ~/.bashrc

rm -f tmp.out
for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    grep PLS $ccreDir/$dnase*bed > tmp.pls
    grep dELS $ccreDir/$dnase*bed > tmp.dels

    grep P-Ed Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt | awk '{print $1 "\t" $4}' | sort -u > x
    awk 'FNR==NR {x[$1];next} ($4 in x)' x tmp.dels | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' | \
	sort -k1,1 -k2,2n | bedtools closest -d -a stdin -b $tss | \
	awk '{print $4 "\t" $(NF-1) "\t" $NF}' | sort -u > y
    awk '{if ($NF > 10000) print $0}' y > y.10k

    awk 'FNR==NR {x[$1$2];next} ($1$2 in x)' x y | awk '{print $1}' | sort -u > out.closest-any-linked

    awk 'FNR==NR {x[$1];next} !($1 in x)' out.closest-any-linked y | awk 'FNR==NR {x[$1];next} ($4 in x)' - $ccreDir/$dnase*bed | bedtools closest -d -a stdin -b ~/Lab/Reference/Human/hg38/GENCODE40/TSS.Basic-PC.bed | awk '{print $4 "\t" $(NF-1) "\t" $NF}' | sort -u > yy
    awk 'FNR==NR {x[$1$2];next} ($1$2 in x)' x yy | awk '{print $1}' | sort -u > out.closest-pc-linked  
    #awk 'FNR==NR {x[$1];next} !($1 in x)' out.closest-pc-linked yy | awk '{if ($NF < 10000) print $0}' | awk '{print $1}' | sort -u  > out.10k
    awk '{if ($NF > 10000) print $0}' yy > yy.10k

    cat out.* | awk 'FNR==NR {x[$1];next} !($1 in x)' - x > tmp.round2
    awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.round2 $ccreDir/$dnase*bed | bedtools closest -a stdin -b ~/Lab/Reference/Human/hg38/MANEv1.5/TSS.MANE.bed | awk '{print $4 "\t" $NF}' > tmp.tech
    awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.round2 $ccreDir/$dnase*bed | bedtools closest -a stdin -b ~/Lab/Reference/Human/hg38/MANEv1.5/Genes.MANE.bed | awk '{print $4 "\t" $15}' >> tmp.tech
    awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.round2 $ccreDir/$dnase*bed | bedtools closest -a stdin -b ~/Lab/Reference/Human/hg38/GENCODE40/Genes.Basic.bed | awk '{print $4 "\t" $15}' >> tmp.tech
    awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.round2 $ccreDir/$dnase*bed | bedtools closest -a stdin -b ~/Lab/Reference/Human/hg38/GENCODE40/Genes.Basic-PC.bed | awk '{print $4 "\t" $15}' >> tmp.tech

    awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' tmp.round2 tmp.tech | awk '{print $1}' | sort -u > out.tech

    cat out.* | awk 'FNR==NR {x[$1];next} !($1 in x)' - x > tmp.round3

    awk 'FNR==NR {x[$1];next} ($1 in x)' tmp.round3 yy | awk '{if ($NF < 10000) print $0}' | awk '{print $1}' | sort -u  > out.10k

    cat out.* | awk 'FNR==NR {x[$1];next} !($1 in x)' - x > tmp.round4 

    grep CTCF-Anchors Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt | grep P-Ed | awk '{print $1}' | \
	sort -u | awk 'FNR==NR {x[$1];next} !($1 in x)' out.10k - > a
    grep -v CTCF-Anchors Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt | grep P-Ed | awk '{print $1}' | \
	sort -u | awk 'FNR==NR {x[$1];next} !($1 in x)' out.10k - > b
    awk 'FNR==NR {x[$0];next} ($0 in x)' a b  > tmp.both
    awk 'FNR==NR {x[$0];next} !($0 in x)' a b  > tmp.uniq-b
    awk 'FNR==NR {x[$0];next} !($0 in x)' b a  > tmp.uniq-a

    tboth=$(wc -l tmp.both | awk '{print $1}')
    tb=$(wc -l tmp.uniq-b | awk '{print $1}')
    ta=$(wc -l tmp.uniq-a | awk '{print $1}')


    cat out.closest-any-linked out.closest-pc-linked out.tech | \
	awk 'FNR==NR {x[$1];next} ($1 in x)' - yy | awk '{print $0 "\t" "closest"}' > tmp.histogram
    cat out.closest-any-linked out.closest-pc-linked out.tech | \
        awk 'FNR==NR {x[$1];next} !($1 in x)' - yy | awk '{print $0 "\t" "not-closest"}' >> tmp.histogram

    mv tmp.histogram Closest-Gene-Evaluation/$cell.Distance-Histogram.txt

    a=$(wc -l out.closest-any-linked | awk '{print $1}')
    b=$(wc -l out.closest-pc-linked | awk '{print $1}')
    bb=$(wc -l out.10k | awk '{print $1}')
    c=$(wc -l out.tech | awk '{print $1}')
    d=$(awk '{print $1}' tmp.round4 | sort -u | wc -l | awk '{print $1}')
    both=$(cat tmp.round4 | awk 'FNR==NR {x[$1];next} ($1 in x)' - tmp.both | wc -l | awk '{print $1}')
    ub=$(cat tmp.round4 | awk 'FNR==NR {x[$1];next} ($1 in x)' - tmp.uniq-b | wc -l | awk '{print $1}')
    ua=$(cat tmp.round4 | awk 'FNR==NR {x[$1];next} ($1 in x)' - tmp.uniq-a | wc -l | awk '{print $1}')

    echo -e $cell "\t" $a "\t" $b "\t" $c "\t" $bb "\t" $d "\t" $tboth "\t" $both "\t" $tb "\t" $ub "\t" $ta "\t" $ua

    rnaExp=$(awk '{if (NR == '$i') print $16}' $list)
    tsvs=$(awk '{if (NR == '$i') print $17}' $list | awk '{gsub(/;/,"\t");print}' )

    echo "" > tmp.matrix
    for tsv in ${tsvs[@]}
    do
        awk '{print $1 "\t" $6}' $dataDir/$rnaExp/$tsv.tsv | paste tmp.matrix - > tmp.tmp
        mv tmp.tmp tmp.matrix
    done
    awk '{sum=0; for(i=2;i<=NF;i+=2) sum += $i; print $1 "\t" sum/(NF/2)}' tmp.matrix > tmp.rna-quants

    awk '{if ($NF > 1) print $0}' tmp.rna-quants > tmp.tpm5
    tss2=~/Lab/Reference/Human/hg38/GENCODE40/TSS.Basic.bed
    awk '{print $NF "\t" $0}'  $tss2 | awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' tmp.tpm5 - | awk '{printf "%s", $2; for(i=3;i<=NF;i+=1) printf "\t%s",$i ; print ""}' > tmp.filtered


    awk 'FNR==NR {x[$1];next} ($4 in x)' x tmp.dels | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' | \
	awk 'FNR==NR {x[$1];next} ($4 in x)' y.10k - | sort -k1,1 -k2,2n | \
        bedtools closest -d -a stdin -b tmp.filtered | awk '{print $4 "\t" $(NF-1) "\t" $NF}' | sort -u > y5


    awk 'FNR==NR {x[$1];next} ($1 in x)' y.10k x > x10k

    t1=$(wc -l y5 | awk '{print $1}')
    awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' x10k y5 | wc -l | awk '{print $1 "\t" "'$t1'" "\t" $1/'$t1'}'
    t2=$(wc -l x10k | awk '{print $1}')
    awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' y5 x10k | wc -l | awk '{print $1 "\t" "'$t2'" "\t" $1/'$t2'}'

groups=(tmp.both tmp.uniq-a tmp.uniq-b)
for group in ${groups[@]}
do
echo $group
t1=$(awk 'FNR==NR {x[$1];next} ($1 in x)' $group y5 | wc -l | awk '{print $1}')
awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' x10k y5 | awk 'FNR==NR {x[$1];next} ($1 in x)' $group - | wc -l | awk '{print $1 "\t" "'$t1'" "\t" $1/'$t1'}'
t2=$(awk 'FNR==NR {x[$1];next} ($1 in x)' $group x10k | wc -l | awk '{print $1}')
awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' y5 x10k | awk 'FNR==NR {x[$1];next} ($1 in x)' $group - | wc -l | awk '{print $1 "\t" "'$t2'" "\t" $1/'$t2'}'
done



rm out.* tmp.* x* y*
done

