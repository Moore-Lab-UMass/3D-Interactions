#Jill Moore
#UMass Chan
#October 2026

workingDir=~/Moore-Lab/Papers/3D-Genome

list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
dataDir=/data/projects/encode/data/

#cd $workingDir/1_General-Comparisons

rm -f running.out
for i in `seq 1 1 6`
do
    rm -f running.out1 running.out2
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    echo $cell

    rnap2=$(awk '{if (NR == '$i') print $2}' $list)
    ctcf=$(awk '{if (NR == '$i') print $4}' $list)
    hic=$(awk -F "\t" '{if (NR == '$i') print $19}' $list)
    dnase=$(awk '{if (NR == '$i') print $8"-"$9}' $list)

    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    awk '{if ($2 > 2.32) print $0}' $registryDir/signal-output/$ctcfChIP.txt | 
        awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.ctcf
    ctcfTotal=$(wc -l tmp.ctcf | awk '{print $1}')

    assays=(RNAPII CTCF HiC)
    for assay in ${assays[@]}
    do
	echo -e "\t" $assay
        links=$workingDir/1_General-Comparisons/Assay-Overlap-Details/$cell.$assay.txt

        arch=$(bedtools pairtobed -type both -a $links -b tmp.ctcf | \
            awk '{print $1 "\t" $2 "\t" $3 "\n" $4 "\t" $5 "\t" $6}' | \
	    bedtools intersect -u -a tmp.ctcf -b stdin | wc -l | awk '{print $1}')

	bedtools pairtobed -type both -a $links -b tmp.ctcf | \
	    awk '{print $1 "\t" $2 "\t" $3 "\n" $4 "\t" $5 "\t" $6}' | \
	     bedtools intersect -u -a tmp.ctcf -b stdin | awk '{print $4 "\t" "'$assay'"}' >> running.out1

        any=$(awk '{print $1 "\t" $2 "\t" $3 "\n" $4 "\t" $5 "\t" $6}' $links | \
            bedtools intersect -u -a tmp.ctcf -b stdin | wc -l | awk '{print $1}')

	awk '{print $1 "\t" $2 "\t" $3 "\n" $4 "\t" $5 "\t" $6}' $links | \
	    bedtools intersect -u -a tmp.ctcf -b stdin | awk '{print $4 "\t" "'$assay'"}' >> running.out2 

	echo -e $cell "\t" $assay "\t" $ctcfTotal "\t" $arch "\t" $any >> running.out
    done
    
    arch=$(awk '{print $1}' running.out1 | sort -u | wc -l | awk '{print $1}')
    any=$(awk '{print $1}' running.out2 | sort -u | wc -l | awk '{print $1}')

    echo -e $cell "\t" "all" "\t" $ctcfTotal "\t" $arch "\t" $any >> running.out

    #DNase signal
    awk 'FNR==NR {x[$1];next} ($1 in x)' running.out1 $registryDir/signal-output/$dnase.txt > tmp.yes
    awk 'FNR==NR {x[$1];next} !($4 in x)' running.out1 tmp.ctcf | awk 'FNR==NR {x[$4];next} ($1 in x)' - $registryDir/signal-output/$dnase.txt > tmp.no

    yZ=$(sort -k2,2g tmp.yes | awk '{a[i++]=$2;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')
    ySig=$(sort -k3,3g tmp.yes | awk '{a[i++]=$3;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')

    nZ=$(sort -k2,2g tmp.no | awk '{a[i++]=$2;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')
    nSig=$(sort -k3,3g tmp.no | awk '{a[i++]=$3;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')

    p=$(Rscript $scriptDir/signal-p.R | awk '{print $1}')

    echo -e $cell "\t" $yZ "\t" $nZ "\t" $ySig "\t" $nSig "\t" $p >> running.dnase

    #CTCF signal
    awk 'FNR==NR {x[$1];next} ($1 in x)' running.out1 $registryDir/signal-output/$ctcfChIP.txt > tmp.yes
    awk 'FNR==NR {x[$1];next} !($4 in x)' running.out1 tmp.ctcf | awk 'FNR==NR {x[$4];next} ($1 in x)' - $registryDir/signal-output/$ctcfChIP.txt > tmp.no

    yZ=$(sort -k2,2g tmp.yes | awk '{a[i++]=$2;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')
    ySig=$(sort -k3,3g tmp.yes | awk '{a[i++]=$3;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')


    nZ=$(sort -k2,2g tmp.no | awk '{a[i++]=$2;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')
    nSig=$(sort -k3,3g tmp.no | awk '{a[i++]=$3;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')

    p=$(Rscript $signalDir/signal-p.R | awk '{print $1}')

    echo -e $cell "\t" $yZ "\t" $nZ "\t" $ySig "\t" $nSig "\t" $p >> running.ctcf

done

mv running.out Summary-CTCF-Bound-Interactions.txt
mv running.ctcf Summary-CTCF-Bound-Interactions-CTCF-Signal.txt
mv running.dnase Summary-CTCF-Bound-Interactions-DNase-Signal.txt


i=4 #K562

rad21=$workingDir/1_General-Comparisons/Cohesin-Comparison/ENCFF930WPG.bed

rm -f running.out1
cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
echo $cell

rnap2=$(awk '{if (NR == '$i') print $2}' $list)
ctcf=$(awk '{if (NR == '$i') print $4}' $list)
hic=$(awk -F "\t" '{if (NR == '$i') print $19}' $list)

ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
awk '{if ($2 > 2.32) print $0}' $registryDir/signal-output/$ctcfChIP.txt | 
    awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.ctcf
ctcfTotal=$(wc -l tmp.ctcf | awk '{print $1}')

assays=(RNAPII CTCF HiC)
for assay in ${assays[@]}
do
    echo -e "\t" $assay
    links=$workingDir/1_General-Comparisons/Assay-Overlap-Details/$cell.$assay.txt

    bedtools pairtobed -type both -a $links -b tmp.ctcf | \
        awk '{print $1 "\t" $2 "\t" $3 "\n" $4 "\t" $5 "\t" $6}' | \
	bedtools intersect -u -a tmp.ctcf -b stdin | awk '{print $4 "\t" "'$assay'"}' >> running.out1

done
    
arch=$(awk '{print $1}' running.out1 | sort -u | wc -l | awk '{print $1}')
yO=$(awk 'FNR==NR {x[$1];next} ($4 in x)' running.out1 tmp.ctcf | bedtools intersect -u -a stdin -b $rad21 | wc -l | awk '{print $1}')

noarch=$(awk 'FNR==NR {x[$1];next} !($4 in x)' running.out1 tmp.ctcf | wc -l | awk '{print $1}')
nO=$(awk 'FNR==NR {x[$1];next} !($4 in x)' running.out1 tmp.ctcf | bedtools intersect -u -a stdin -b $rad21 | wc -l | awk '{print $1}')

echo -e $cell "\t" $arch "\t" $yO "\t" $noarch "\t" $nO > running.rad21



i=3 #HepG2

rad21=$workingDir/1_General-Comparisons/Cohesin-Comparison/ENCFF304NRB.bed

rm -f running.out1
cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
echo $cell

rnap2=$(awk '{if (NR == '$i') print $2}' $list)
ctcf=$(awk '{if (NR == '$i') print $4}' $list)
hic=$(awk -F "\t" '{if (NR == '$i') print $19}' $list)

ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
awk '{if ($2 > 2.32) print $0}' $registryDir/signal-output/$ctcfChIP.txt |
    awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.ctcf
ctcfTotal=$(wc -l tmp.ctcf | awk '{print $1}')

assays=(RNAPII CTCF HiC)
for assay in ${assays[@]}
do
    echo -e "\t" $assay
    links=$workingDir/1_General-Comparisons/Assay-Overlap-Details/$cell.$assay.txt

    bedtools pairtobed -type both -a $links -b tmp.ctcf | \
        awk '{print $1 "\t" $2 "\t" $3 "\n" $4 "\t" $5 "\t" $6}' | \
        bedtools intersect -u -a tmp.ctcf -b stdin | awk '{print $4 "\t" "'$assay'"}' >> running.out1

done

arch=$(awk '{print $1}' running.out1 | sort -u | wc -l | awk '{print $1}')
yO=$(awk 'FNR==NR {x[$1];next} ($4 in x)' running.out1 tmp.ctcf | bedtools intersect -u -a stdin -b $rad21 | wc -l | awk '{print $1}')

noarch=$(awk 'FNR==NR {x[$1];next} !($4 in x)' running.out1 tmp.ctcf | wc -l | awk '{print $1}')
nO=$(awk 'FNR==NR {x[$1];next} !($4 in x)' running.out1 tmp.ctcf | bedtools intersect -u -a stdin -b $rad21 | wc -l | awk '{print $1}')

echo -e $cell "\t" $arch "\t" $yO "\t" $noarch "\t" $nO >> running.rad21

mv running.output Summary-CTCF-Bound-Interactions-RAD21-Overlap.txt


rm -f tmp* running.*
