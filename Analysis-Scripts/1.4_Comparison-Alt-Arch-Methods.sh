#Jill Moore
#UMass Chan
#October 2026

workingDir=~/Moore-Lab/Papers/3D-Genome

list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
dataDir=/data/projects/encode/data/

cd $workingDir/1_General-Comparisons/


rm -f tmp.out
    
i=4 #K562

cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')

rnap2=$(awk '{if (NR == '$i') print $2}' $list)
ctcf=$(awk '{if (NR == '$i') print $4}' $list)
hic=$(awk -F "\t" '{if (NR == '$i') print $19}' $list)

ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
awk '{if ($2 > 2.32) print $0}' $registryDir/signal-output/$ctcfChIP.txt | 
    awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.ctcf

assays=(RNAPII CTCF HiC)
for assay in ${assays[@]}
do
    links=Assay-Overlap-Details/$cell.$assay.txt
    totalLinks=$(wc -l $links | awk '{print $1}')
    bedtools pairtobed -type both -a $links -b tmp.ctcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | \
	sort -u > tmp.$assay.CTCF
    bedtools pairtobed -type both -a $links -b Cohesin-Comparison/ENCFF930WPG.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | \
        sort -u > tmp.$assay.RAD21

    a=$(wc -l tmp.$assay.CTCF | awk '{print $1}')
    b=$(wc -l tmp.$assay.RAD21 | awk '{print $1}')
    awk 'FNR==NR {x[$0];next} ($0 in x)' tmp.$assay.CTCF tmp.$assay.RAD21 | wc -l | awk '{print "'$cell'" "\t" "'$assay'" "\t" "'$a'" "\t" "'$b'"  "\t" $1}' > running.out

done


i=3 #HepG2

cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')

rnap2=$(awk '{if (NR == '$i') print $2}' $list)
ctcf=$(awk '{if (NR == '$i') print $4}' $list)
hic=$(awk -F "\t" '{if (NR == '$i') print $19}' $list)

ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
awk '{if ($2 > 2.32) print $0}' $registryDir/signal-output/$ctcfChIP.txt |
    awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.ctcf

assays=(RNAPII CTCF HiC)
for assay in ${assays[@]}
do
    links=Assay-Overlap-Details/$cell.$assay.txt
    totalLinks=$(wc -l $links | awk '{print $1}')
    bedtools pairtobed -type both -a $links -b tmp.ctcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | \
        sort -u > tmp.$assay.CTCF
    bedtools pairtobed -type both -a $links -b Cohesin-Comparison/ENCFF304NRB.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | \
        sort -u > tmp.$assay.RAD21

    a=$(wc -l tmp.$assay.CTCF | awk '{print $1}')
    b=$(wc -l tmp.$assay.RAD21 | awk '{print $1}')
    awk 'FNR==NR {x[$0];next} ($0 in x)' tmp.$assay.CTCF tmp.$assay.RAD21 | wc -l | awk '{print "'$cell'" "\t" "'$assay'" "\t" "'$a'" "\t" "'$b'"  "\t" $1}' >> running.out
done

mv running.out Summary-Comparison-RAD21-Overlap.txt

rm -f running.out
for i in `seq 1 1 6`
do
    rm -f running.out1 running.out2
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    echo $cell

    rnap2=$(awk '{if (NR == '$i') print $2}' $list)
    ctcf=$(awk '{if (NR == '$i') print $4}' $list)
    hic=$(awk -F "\t" '{if (NR == '$i') print $19}' $list)

    dnase=$(awk '{if (NR == '$i') print $9}' $list)

    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    awk '{if ($2 > 2.32) print $0}' $registryDir/signal-output/$ctcfChIP.txt | 
        awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.ctcf
    ctcfTotal=$(wc -l tmp.ctcf | awk '{print $1}')

    cts=$ccreDir/$dnase*bed


    assays=(RNAPII CTCF HiC) #Update to stratify by assay
    for assay in ${assays[@]}
    do
        echo -e "\t" $assay
        links=$workingDir/1_General-Comparisons/Assay-Overlap-Details/$cell.$assay.txt

        bedtools pairtobed -type both -a $links -b tmp.ctcf >> running.out1

    done

    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' running.out1 | \
        bedtools pairtobed -type both -a stdin -b tmp.ctcf > tmp.links

    totalLinks=$(wc -l tmp.links | awk '{print $1}')

    awk '{print $1 "\t" $2 "\t" $3 "\n" $4 "\t" $5 "\t" $6}' tmp.links | \
        bedtools intersect -u -a tmp.ctcf -b stdin | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > tmp.bed

    ~/bin/twoBitToFa ~/Lab/Reference/Human/hg38/hg38.2bit -bed=tmp.bed tmp.fa
    /data/common/tools/meme-5.1.1/bin/fimo --oc tmp-fimo Cohesin-Comparison/CTCF.H14CORE.0.P.B_meme_format.meme tmp.fa

    awk -F "\t" '{if ($6 == "+") print $3}' tmp-fimo/fimo.tsv | \
        awk 'FNR==NR {x[$1];next} ($4 in x)' - tmp.ctcf > tmp.pos 
    awk -F "\t" '{if ($6 == "-") print $3}' tmp-fimo/fimo.tsv | \
        awk 'FNR==NR {x[$1];next} ($4 in x)' - tmp.ctcf > tmp.neg

    bothAnchor=$(bedtools intersect -u -a tmp.links -b tmp.pos tmp.neg | awk '{print $4 "\t" $5 "\t" $6}' | \
	bedtools intersect -u -a - -b tmp.pos tmp.neg | wc -l | awk '{print $1}')

    bedtools intersect -u -a tmp.links -b tmp.pos | awk '{print $4 "\t" $5 "\t" $6}' | \
        bedtools intersect -u -a - -b tmp.neg | wc -l | \
	awk '{print "'$cell'" "\t" "'$totalLinks'" "\t" "'$bothAnchor'" "\t" $1}' >> running.out

done

mv running.out Summary-Convergent-Motifs.txt

rm -f -r tmp*

