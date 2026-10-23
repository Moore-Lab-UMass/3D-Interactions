#Jill Moore
#UMass Chan
#October 2026

workingDir=~/Moore-Lab/Papers/3D-Genome

list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
dataDir=/data/projects/encode/data/

source ~/.bashrc

cd $workingDir/1_General-Comparisons
totalLinks=0

rm -f tmp.out running.*
for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    echo $cell
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    grep -v Low-DNase $registryDir/Cell-Type-Specific/Individual-Files/$dnase*bed | awk '{print $1 "\t" $2 "\t" $3}' > tmp.ccre

    awk '{if ($2 > 1.64) print $0}' $registryDir/signal-output/$ctcfChIP.txt | \
        awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed | awk '{print $1 "\t" $2 "\t" $3}'  >> tmp.ccre

    assays=(RNAPII CTCF HiC)
    for assay in ${assays[@]}
    do
	echo -e "\t" $assay
        links=Assay-Overlap-Details/$cell.$assay.txt
	totalLinks=$(wc -l $links | awk '{print $1 + "'$totalLinks'"}')
	awk '{print $1 "\t" $2 "\t" $3 "\t" $7 "\n" $4 "\t" $5 "\t" $6 "\t" $7}' $links > tmp.links
	bedtools intersect -u -a tmp.ccre -b tmp.links >> running.ccre

	awk '{print $3-$2 "\n" $6-$5}' $links | sort -k1,1rg | \
	    awk 'BEGIN{min=1000000; max=0}{a[i++]=$1;if ($1 < min) min = $1; if ($1 > max) max = $1}\
	    END{x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2 "\t" min "\t" max; else print a[x-1] "\t" min "\t" max;}' | \
	    awk '{print "'$cell'" "\t" "'$assay'" "\t" $0}' >> running.anchor

	awk '{print ($5+$6)/2-($2+$3)/2}' $links | sort -k1,1rg | \
            awk 'BEGIN{min=1000000; max=0}{a[i++]=$1;if ($1 < min) min = $1; if ($1 > max) max = $1}\
            END{x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2 "\t" min "\t" max; else print a[x-1] "\t" min "\t" max;}' | \
            awk '{print "'$cell'" "\t" "'$assay'" "\t" $0}' >> running.distance

    done
done

sort -u running.ccre | wc -l | awk '{print "Total interactions: '$totalLinks'" "\n" "Total cCREs: " $1}'
mv running.anchor Summary-Anchor-Stats.txt
mv running.distance Summary-Distance-Stats.txt
rm -f running.* tmp*
