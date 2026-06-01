
workingDir=~/Moore-Lab/Papers/3D-Genome

list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
ccreDir=$registryDir/Cell-Type-Specific/Individual-Files
dataDir=/data/projects/encode/data/

cd $workingDir/2_Promoter-Centric
mkdir -p Promoter-Centric-Annotations

for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    echo $cell
    grep PLS $ccreDir/$dnase*bed > tmp.pls

    awk '{if ($2 > 1.64) print $0}' $registryDir/signal-output/$ctcfChIP.txt | \
	awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.ctcf

    rm -f tmp.combo tmp.all tmp.summary
    assays=(RNAPII CTCF HiC)
    for assay in ${assays[@]}
    do
        cat ../1_General-Comparisons/Assay-Overlap-Details/$cell.$assay.txt | \
	    awk '{print $0 "\t" "'$assay'"}' >> tmp.all
    done

    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" "id-"NR"_"$7"__"$8}' tmp.all > tmp.combo

    bedtools pairtobed -type both -a tmp.combo -b tmp.pls | \
	awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7  "\t" "P-P"}' | \
	sort -u > tmp.summary

    bedtools pairtobed -type xor -a tmp.combo -b tmp.pls | \
	awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' | \
	sort -u > tmp.one

    awk '{print $1 "\t" $2 "\t" $3 "\t" $7 "\n" $4 "\t" $5 "\t" $6 "\t" $7}' tmp.one > tmp.links

     bedtools intersect -u -a tmp.links -b tmp.pls > tmp.overlap
     awk 'FNR==NR {x[$1$2$3$4];next} !($1$2$3$4 in x)' tmp.overlap tmp.links > tmp.candidates
     
     grep pELS $ccreDir/$dnase*bed | bedtools intersect -u -a tmp.candidates -b stdin | \
	 awk 'FNR==NR {x[$4];next} ($7 in x)' - tmp.combo | \
	 awk '{print $0 "\t" "P-Ep"}' >> tmp.summary
     grep pELS $ccreDir/$dnase*bed | bedtools intersect -v -a tmp.candidates -b stdin > tmp.no1
     
     grep dELS $ccreDir/$dnase*bed | bedtools intersect -u -a tmp.no1 -b stdin | \
         awk 'FNR==NR {x[$4];next} ($7 in x)' - tmp.combo | \
         awk '{print $0 "\t" "P-Ed"}' >> tmp.summary
     grep dELS $ccreDir/$dnase*bed | bedtools intersect -v -a tmp.no1 -b stdin > tmp.no2

     cat tmp.ctcf | bedtools intersect -u -a tmp.no2 -b stdin | \
	 awk 'FNR==NR {x[$4];next} ($7 in x)' - tmp.combo | \
         awk '{print $0 "\t" "P-C"}' >> tmp.summary
     cat tmp.ctcf | bedtools intersect -v -a tmp.no2 -b stdin > tmp.no1

     grep -v Low-DNase $ccreDir/$dnase*bed | bedtools intersect -u -a tmp.no1 -b stdin | \
	 awk 'FNR==NR {x[$4];next} ($7 in x)' - tmp.combo | \
         awk '{print $0 "\t" "P-O"}' >> tmp.summary
     grep -v Low-DNase $ccreDir/$dnase*bed | bedtools intersect -v -a tmp.no1 -b stdin > tmp.no2

     grep Low-DNase $ccreDir/$dnase*bed | bedtools intersect -u -a tmp.no2 -b stdin | \
         awk 'FNR==NR {x[$4];next} ($7 in x)' - tmp.combo | \
         awk '{print $0 "\t" "P-I"}' >> tmp.summary
     awk 'FNR==NR {x[$7];next} !($7 in x)' tmp.summary tmp.one | awk '{print $0 "\t" "P-N"}' > tmp
     cat tmp >> tmp.summary


     awk '{if ($2 > 2.32) print $0}' $registryDir/signal-output/$ctcfChIP.txt | \
        awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.highctcf
     bedtools pairtobed -type both -a tmp.summary -b tmp.highctcf | \
        awk 'FNR==NR {x[$7];next} ($7 in x)' - tmp.summary > tmp.bothanchor
     bedtools pairtopair -is -type both -a tmp.summary -b tmp.bothanchor | \
	awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "CTCF-Anchors"}' | \
	sort -u > tmp.out
     bedtools pairtopair -is -type notboth -a tmp.summary -b tmp.bothanchor | \
	awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" "No-CTCF"}' | \
	sort -u >> tmp.out

     mv tmp.out Promoter-Centric-Annotations/$cell.Promoter-Links.txt
     rm tmp.*
done
