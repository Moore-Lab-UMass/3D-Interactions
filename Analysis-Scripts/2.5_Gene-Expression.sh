
workingDir=~/Moore-Lab/Papers/3D-Genome

list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
ccreDir=$registryDir/Cell-Type-Specific/Individual-Files
dataDir=/data/projects/encode/data/
scriptDir=~/GitHub/3D-Interactions/Analysis-Scripts/Toolkit

cd $workingDir/2_Promoter-Centric
mkdir -p Promoter-Centric-Annotations

for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    grep PLS $ccreDir/$dnase*bed > tmp.pls

    rnaExp=$(awk '{if (NR == '$i') print $16}' $list)
    tsvs=$(awk '{if (NR == '$i') print $17}' $list | awk '{gsub(/;/,"\t");print}' )

    echo "" > tmp.matrix
    for tsv in ${tsvs[@]}
    do
        awk '{print $1 "\t" $6}' $dataDir/$rnaExp/$tsv.tsv | paste tmp.matrix - > tmp.tmp
        mv tmp.tmp tmp.matrix
    done
    awk '{sum=0; for(i=2;i<=NF;i+=2) sum += $i; print $1 "\t" sum/(NF/2)}' tmp.matrix > tmp.rna-quants


    awk '{print $4 "\n" $5}' Promoter-Promoter-Pairs/$cell.P-P-Pairs.txt | sort -u > tmp.pp
    awk '{print $4}' Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | \
	awk 'FNR==NR {x[$1];next} !($1 in x)' tmp.pp - | sort -u > tmp.pec
    awk 'FNR==NR {x[$4];next} ($1 in x)' tmp.pls $registryDir/PLS-Gene-List.hg38-V4.txt | awk '{print $2}' | \
	awk 'FNR==NR {x[$1];next} !($1 in x)' tmp.pp - | awk 'FNR==NR {x[$1];next} !($1 in x)' tmp.pec - | sort -u > tmp.none

    awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' tmp.pp tmp.rna-quants | awk '{print $0 "\t" "P-P"}' > running.quants
    P=$(awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' tmp.pp tmp.rna-quants | sort -k2,2rg | \
	awk '{a[i++]=$2;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')
    awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' tmp.pec tmp.rna-quants | awk '{print $0 "\t" "P-E-C"}' >> running.quants
    EC=$(awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' tmp.pec tmp.rna-quants | sort -k2,2rg | \
	awk '{a[i++]=$2;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')
    awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' tmp.none tmp.rna-quants | awk '{print $0 "\t" "None"}' >> running.quants
    N=$(awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' tmp.none tmp.rna-quants | sort -k2,2rg | \
        awk '{a[i++]=$2;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')

    echo -e $cell "\t" $P "\t" $EC "\t" $N

    mv running.quants Expression/$cell.Expression.txt

    rm tmp*

done

