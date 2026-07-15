#Jill Moore
#UMass Chan
#July 2026

workingDir=~/Moore-Lab/Papers/3D-Genome

list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
dataDir=/data/projects/encode/data/

cd $workingDir/1_General-Comparisons/Cohesin-Comparison


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
    links=../Assay-Overlap-Details/$cell.$assay.txt
    totalLinks=$(wc -l $links | awk '{print $1}')
    bedtools pairtobed -type both -a $links -b tmp.ctcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | \
	sort -u > tmp.$assay.CTCF
    bedtools pairtobed -type both -a $links -b ENCFF930WPG.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | \
        sort -u > tmp.$assay.RAD21

    a=$(wc -l tmp.$assay.CTCF | awk '{print $1}')
    b=$(wc -l tmp.$assay.RAD21 | awk '{print $1}')
    c=$(awk 'FNR==NR {x[$0];next} ($0 in x)' tmp.$assay.CTCF tmp.$assay.RAD21 | wc -l | awk '{print $1}') 

    bedtools intersect -u -a tmp.$assay.CTCF -b ctcf-motifs.plus.bed | awk '{print $4 "\t" $5 "\t" $6 "\t" $7}' | \
	bedtools intersect -u -a stdin -b ctcf-motifs.minus.bed | wc -l | \
	awk '{print "'$cell'" "\t" "'$assay'" "\t" "'$a'" "\t" "'$b'"  "\t" "'$c'" "\t" $1}' >> running.out

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
    links=../Assay-Overlap-Details/$cell.$assay.txt
    totalLinks=$(wc -l $links | awk '{print $1}')
    bedtools pairtobed -type both -a $links -b tmp.ctcf | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | \
        sort -u > tmp.$assay.CTCF
    bedtools pairtobed -type both -a $links -b ENCFF304NRB.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | \
        sort -u > tmp.$assay.RAD21

    a=$(wc -l tmp.$assay.CTCF | awk '{print $1}')
    b=$(wc -l tmp.$assay.RAD21 | awk '{print $1}')
    c=$(awk 'FNR==NR {x[$0];next} ($0 in x)' tmp.$assay.CTCF tmp.$assay.RAD21 | wc -l | awk '{print $1}')

    bedtools intersect -u -a tmp.$assay.CTCF -b ctcf-motifs.plus.bed | awk '{print $4 "\t" $5 "\t" $6 "\t" $7}' | \
        bedtools intersect -u -a stdin -b ctcf-motifs.minus.bed | wc -l | \
        awk '{print "'$cell'" "\t" "'$assay'" "\t" "'$a'" "\t" "'$b'"  "\t" "'$c'" "\t" $1}' >> running.out
done

mv running.out Summary-Comparison-Arch-Methods.txt
rm -f tmp*

