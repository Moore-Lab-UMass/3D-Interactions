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

cd $workingDir/2_Promoter-Centric

source ~/.bashrc

rm -f tmp.out

echo "Running enhancer..."
for i in `seq 1 1 6`
do
    echo $i
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

    echo -e $cell "\t" $a "\t" $b "\t" $c "\t" $bb "\t" $d "\t" $tboth "\t" $both "\t" $tb "\t" $ub "\t" $ta "\t" $ua >> running.summary

rm out.* tmp.* x* y*
done
mv running.summary ../3_Candidate-Architectural/Summary-Enhancer-Closest-Gene.txt

echo "Running promoter..."
rm -f tmp.out
for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    grep PLS $ccreDir/$dnase*bed > tmp.pls

    cat Promoter-Promoter-Pairs/$cell.P-P-Pairs.txt | awk '{print $1 "\t" $4 "\t" $5 "\n" $2 "\t" $5 "\t" $4}' | sort -u > tmp.assignments
    totalPLS=$(awk '{print $1}' tmp.assignments | sort -u | wc -l | awk '{print $1}')

    #All TSSs
    awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.assignments tmp.pls | \
	awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' | sort -k1,1 -k2,2n | bedtools closest -k 100 -d -a stdin -b $tss > tmp.100
    awk 'FNR==NR {x[$1$2];next} !($4$11 in x)' ~/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/PLS-Gene-List.hg38-V4.txt tmp.100 > tmp.filtered
    python filter.py tmp.filtered > tmp.all-tss-closest
    wc -l tmp.all-tss-closest | awk '{print $1 "\t" "'$totalPLS'"}'
    awk 'FNR==NR {x[$1$3];next} ($1$2 in x)' tmp.assignments tmp.all-tss-closest > out.round-1

    #PC TSSs
    awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.assignments tmp.pls | \
	awk 'FNR==NR {x[$1];next} !($4 in x)' out.round-1 - | \
        awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' | sort -k1,1 -k2,2n | bedtools closest -k 100 -d -a stdin -b ~/Lab/Reference/Human/hg38/GENCODE40/TSS.Basic-PC.bed > tmp.100
    awk 'FNR==NR {x[$1$2];next} !($4$11 in x)' ~/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/PLS-Gene-List.hg38-V4.txt tmp.100 > tmp.filtered
    python filter.py tmp.filtered > tmp.pc-tss-closest
    awk 'FNR==NR {x[$1$3];next} ($1$2 in x)' tmp.assignments tmp.pc-tss-closest > out.round-2

    #Other annotations
    cat out.* | awk 'FNR==NR {x[$1];next} !($1 in x)' - tmp.assignments > tmp.beyond-tss
    awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.beyond-tss tmp.pls | bedtools closest -k 2 -a stdin -b ~/Lab/Reference/Human/hg38/MANEv1.5/Genes.MANE.bed | awk '{print $4 "\t" $15}' > tmp.tech
    awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.beyond-tss tmp.pls | bedtools closest -k 2 -a stdin -b ~/Lab/Reference/Human/hg38/GENCODE40/Genes.Basic.bed | awk '{print $4 "\t" $15}' >> tmp.tech
    awk 'FNR==NR {x[$1];next} ($4 in x)' tmp.beyond-tss tmp.pls | bedtools closest -k 2 -a stdin -b ~/Lab/Reference/Human/hg38/GENCODE40/Genes.Basic-PC.bed | awk '{print $4 "\t" $15}' >> tmp.tech

    awk -F "." 'FNR==NR {x[$1];next} !($1 in x)' tmp.beyond-tss tmp.tech > tmp.filtered
    awk '{print $1 "\t" $3}' tmp.beyond-tss | awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' - tmp.filtered | awk '{print $1}' | sort -u > out.round-3

    cat out.* | awk 'FNR==NR {x[$1];next} !($1 in x)' - tmp.all-tss-closest | awk '{if ($NF < 10000) print $0}' > out.round-4
    cat out.* | awk 'FNR==NR {x[$1];next} !($1 in x)' - tmp.assignments | awk '{print $1}' | sort -u > out.remaining

    grep CTCF-Anchors Promoter-Promoter-Pairs/$cell.P-P-Pairs.txt | awk '{print $1 "\n" $2}' | \
	sort -u | awk 'FNR==NR {x[$1];next} !($1 in x)' out.round-4 - > tmp.can-arch
    grep -v CTCF-Anchors Promoter-Promoter-Pairs/$cell.P-P-Pairs.txt | awk '{print $1 "\n" $2}' | \
	sort -u | awk 'FNR==NR {x[$1];next} !($1 in x)' out.round-4 - > tmp.non-arch
    awk 'FNR==NR {x[$0];next} ($0 in x)' tmp.can-arch tmp.non-arch  > tmp.both
    awk 'FNR==NR {x[$0];next} !($0 in x)' tmp.non-arch tmp.can-arch  > tmp.arch-only
    awk 'FNR==NR {x[$0];next} !($0 in x)' tmp.can-arch tmp.non-arch  > tmp.no-only

    tboth=$(wc -l tmp.both | awk '{print $1}')
    tb=$(wc -l tmp.arch-only | awk '{print $1}')
    ta=$(wc -l tmp.no-only | awk '{print $1}')

    a=$(wc -l out.round-1 | awk '{print $1}')
    b=$(wc -l out.round-2 | awk '{print $1}')
    bb=$(wc -l out.round-4 | awk '{print $1}')
    c=$(wc -l out.round-3 | awk '{print $1}')
    d=$(awk '{print $1}' out.remaining | sort -u | wc -l | awk '{print $1}')
    both=$(cat out.remaining | awk 'FNR==NR {x[$1];next} ($1 in x)' - tmp.both | wc -l | awk '{print $1}')
    ub=$(cat out.remaining | awk 'FNR==NR {x[$1];next} ($1 in x)' - tmp.arch-only | wc -l | awk '{print $1}')
    ua=$(cat out.remaining | awk 'FNR==NR {x[$1];next} ($1 in x)' - tmp.no-only | wc -l | awk '{print $1}')

    echo -e $cell "\t" $a "\t" $b "\t" $c "\t" $bb "\t" $d "\t" $tboth "\t" $both "\t" $tb "\t" $ub "\t" $ta "\t" $ua >> running.summary

    rm out.* tmp.*
done
mv running.summary ../3_Candidate-Architectural/Summary-Promoter-Closest-Gene.txt


echo "Running CTCF ..."
for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    awk '{if ($2 > 2.32) print $0}' $registryDir/signal-output/$ctcfChIP.txt |
        awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.dels

    grep PLS $ccreDir/$dnase*bed > tmp.pls
    #grep dELS $ccreDir/$dnase*bed > tmp.dels

    grep P-C Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | awk '{print $1 "\t" $4}' | sort -u > x
    awk 'FNR==NR {x[$1];next} ($5 in x)' x tmp.dels | awk '{print $1 "\t" $2 "\t" $3 "\t" $5}' | \
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

    grep CTCF-Anchors Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | awk '{print $1}' | \
	sort -u | awk 'FNR==NR {x[$1];next} !($1 in x)' out.10k - > a
    grep -v CTCF-Anchors Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | awk '{print $1}' | \
	sort -u | awk 'FNR==NR {x[$1];next} !($1 in x)' out.10k - > b
    awk 'FNR==NR {x[$0];next} ($0 in x)' a b  > tmp.both
    awk 'FNR==NR {x[$0];next} !($0 in x)' a b  > tmp.uniq-b
    awk 'FNR==NR {x[$0];next} !($0 in x)' b a  > tmp.uniq-a

    tboth=$(wc -l tmp.both | awk '{print $1}')
    tb=$(wc -l tmp.uniq-b | awk '{print $1}')
    ta=$(wc -l tmp.uniq-a | awk '{print $1}')

    a=$(wc -l out.closest-any-linked | awk '{print $1}')
    b=$(wc -l out.closest-pc-linked | awk '{print $1}')
    bb=$(wc -l out.10k | awk '{print $1}')
    c=$(wc -l out.tech | awk '{print $1}')
    d=$(awk '{print $1}' tmp.round4 | sort -u | wc -l | awk '{print $1}')
    both=$(cat tmp.round4 | awk 'FNR==NR {x[$1];next} ($1 in x)' - tmp.both | wc -l | awk '{print $1}')
    ub=$(cat tmp.round4 | awk 'FNR==NR {x[$1];next} ($1 in x)' - tmp.uniq-b | wc -l | awk '{print $1}')
    ua=$(cat tmp.round4 | awk 'FNR==NR {x[$1];next} ($1 in x)' - tmp.uniq-a | wc -l | awk '{print $1}')

    echo -e $cell "\t" $a "\t" $b "\t" $c "\t" $bb "\t" $d "\t" $tboth "\t" $both "\t" $tb "\t" $ub "\t" $ta "\t" $ua >> running.summary
    rm out.* tmp.* x* y*
done
mv running.summary ../3_Candidate-Architectural/Summary-CTCF-Closest-Gene.txt

