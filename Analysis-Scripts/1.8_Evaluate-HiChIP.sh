#Jill Moore
#UMass Chan
#July 2026


workingDir=~/Moore-Lab/Papers/3D-Genome

list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
ccreDir=$registryDir/Cell-Type-Specific/Individual-Files
dataDir=/data/projects/encode/data/

cd $workingDir/REVISION/HiChIP
source ~/.bashrc

i=4
cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
dnase=$(awk '{if (NR == '$i') print $9}' $list)
grep -v Low-DNase $ccreDir/$dnase*bed > tmp.cts

ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
links=$1
key=$2

totalLinks=$(wc -l $links | awk '{print $1}')


    awk '{if ($2 > 2.32) print $0}' $registryDir/signal-output/$ctcfChIP.txt |
        awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.ctcf

        totalLinks=$(wc -l $links | awk '{print $1}')
        A=$(bedtools pairtobed -type either -a $links -b tmp.ctcf | awk '{print $1$2$3$4$5$6}' | \
            sort -u | wc -l | awk '{print $1}')
        B=$(bedtools pairtobed -type both -a $links -b tmp.ctcf | awk '{print $1$2$3$4$5$6}' | \
            sort -u | wc -l | awk '{print $1}')
echo -e "CTCF-Occupancy" "\t" $links "\t" $A "\t" $B "\t" $totalLinks 


awk '{if ($2 > 1.64) print $0}' $registryDir/signal-output/$ctcfChIP.txt | awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.ctcf

links=$1

totalLinks=$(wc -l $links | awk '{print $1}')

grep PLS tmp.cts > tmp.ccre
pls=$(bedtools pairtobed -type either -a $links -b tmp.ccre | awk '{print $1$2$3$4$5$6}' | \
	    sort -u | wc -l | awk '{print $1}')
	bedtools pairtobed -type neither -a $links -b tmp.ccre > tmp.next1

	grep ELS tmp.cts > tmp.ccre
	els=$(bedtools pairtobed -type either -a tmp.next1 -b tmp.ccre | awk '{print $1$2$3$4$5$6}' | \
	    sort -u | wc -l | awk '{print $1}')
        bedtools pairtobed -type neither -a tmp.next1 -b tmp.ccre > tmp.next2

	grep CTCF tmp.cts | cat tmp.ctcf - | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > tmp.ccre
        ctcf=$(bedtools pairtobed -type either -a tmp.next2 -b tmp.ccre | awk '{print $1$2$3$4$5$6}' | \
	    sort -u | wc -l | awk '{print $1}')
        bedtools pairtobed -type neither -a tmp.next2 -b tmp.ccre > tmp.next1

	other=$(bedtools pairtobed -type either -a tmp.next1 -b tmp.cts | \
	    awk '{print $1$2$3$4$5$6}' | sort -u | wc -l | awk '{print $1}')
        bedtools pairtobed -type neither -a tmp.next1 -b tmp.cts > tmp.next2

        inactive=$(bedtools pairtobed -type either -a tmp.next2 -b $registryDir/hg38-cCREs-Unfiltered.bed | \
            awk '{print $1$2$3$4$5$6}' | sort -u | wc -l | awk '{print $1}')

	none=$(bedtools pairtobed -type neither -a tmp.next2 -b $registryDir/hg38-cCREs-Unfiltered.bed | \
	    awk '{print $1$2$3$4$5$6}' | sort -u | wc -l | awk '{print $1}')

        echo -e "cCRE-Breakdown" "\t" $links "\t" $pls "\t" $els "\t" $ctcf "\t" $other "\t" $inactive "\t" $none 


    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    echo $cell
    grep PLS $ccreDir/$dnase*bed > tmp.pls

    awk '{if ($2 > 1.64) print $0}' $registryDir/signal-output/$ctcfChIP.txt | \
	awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.ctcf

    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $1"_"$2"_"$3"_"$5"_"$6}' $links > tmp.combo

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

mv tmp.out $key.Promoter-Centric-Interactions.txt

rm tmp.*
