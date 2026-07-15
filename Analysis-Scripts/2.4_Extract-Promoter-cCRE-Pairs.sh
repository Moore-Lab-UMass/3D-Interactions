
workingDir=~/Moore-Lab/Papers/3D-Genome/
list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=~/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38
ccreDir=$registryDir/Cell-Type-Specific/Individual-Files
scriptDir=~/GitHub/3D-Interactions/Analysis-Scripts/Toolkit
geneList=$registryDir/PLS-Gene-List.hg38-V4.txt

rm -f tmp.out

cd ~/Moore-Lab/Papers/3D-Genome/2_Promoter-Centric/

mkdir -p Promoter-Promoter-Pairs
mkdir -p Promoter-Enhancer-Pairs
mkdir -p Promoter-CTCF-Pairs

for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    echo $cell

    grep PLS $ccreDir/$dnase*bed > tmp.pls
    grep ELS $ccreDir/$dnase*bed > tmp.els
    awk '{if ($2 > 1.64) print $0}' $registryDir/signal-output/$ctcfChIP.txt | \
	    awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed | \
	    awk '{print $1 "\t" $2 "\t" $3 "\t" $5}' > tmp.ctcf


    grep P-P Promoter-Centric-Annotations/$cell.Promoter-Links.txt > tmp.filtered

    awk '{print $1 "\t" $2 "\t" $3 "\t" $7"__"$8"__"$9 "\n" $4 "\t" $5 "\t" $6 "\t" $7"__"$8"__"$9}' tmp.filtered > tmp.links
    bedtools intersect -wo -a tmp.links -b tmp.pls > tmp.overlap-pls
    python $scriptDir/match-promoter-promoter.py $geneList tmp.overlap-pls > Promoter-Promoter-Pairs/$cell.P-P-Pairs.txt


    grep P-E Promoter-Centric-Annotations/$cell.Promoter-Links.txt > tmp.filtered
    awk '{print $1 "\t" $2 "\t" $3 "\t" $7"__"$8"__"$9 "\n" $4 "\t" $5 "\t" $6 "\t" $7"__"$8"__"$9}' tmp.filtered > tmp.links
    bedtools intersect -wo -a tmp.links -b tmp.pls > tmp.overlap-pls
    bedtools intersect -v -a tmp.links -b tmp.pls | \
        bedtools intersect -wo -a stdin -b tmp.els > tmp.overlap-els

    python $scriptDir/match-promoter-ccre.py $registryDir/PLS-Gene-List.hg38-V4.txt \
        tmp.overlap-pls tmp.overlap-els > Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt


    grep P-C Promoter-Centric-Annotations/$cell.Promoter-Links.txt > tmp.filtered
    awk '{print $1 "\t" $2 "\t" $3 "\t" $7"__"$8"__"$9 "\n" $4 "\t" $5 "\t" $6 "\t" $7"__"$8"__"$9}' tmp.filtered > tmp.links
    bedtools intersect -wo -a tmp.links -b tmp.pls > tmp.overlap-pls
    bedtools intersect -v -a tmp.links -b tmp.pls | \
	bedtools intersect -wo -a stdin -b tmp.ctcf > tmp.overlap-ctcf

    python $scriptDir/match-promoter-ccre.py $registryDir/PLS-Gene-List.hg38-V4.txt \
	tmp.overlap-pls tmp.overlap-ctcf > Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt

    total=$(wc -l tmp.pls | awk '{print $1}')
    pp=$(awk '{print $2}' Promoter-Promoter-Pairs/$cell.P-P-Pairs.txt | sort -u | wc -l | awk '{print $1}')
    pe=$(awk '{print $2}' Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt | sort -u | wc -l | awk '{print $1}')
    pc=$(awk '{print $2}' Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | sort -u | wc -l | awk '{print $1}')

    echo -e $cell "\t" $total "\t" $pp "\t" $pe "\t" $pc >> running.summary

    rm tmp.*
done

mv running.summary Summary-Promoter-Coverage.txt


