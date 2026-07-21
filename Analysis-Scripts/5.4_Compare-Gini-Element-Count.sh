gini=/zata/public_html/users/elhajjajys/cCRE_paper/tissue_specificity/BRUseq_matrix.gene_sense.avg.quantile.gini.tsv
scriptDir=/data/zusers/moorej3/GitHub/3D-Interactions/Analysis-Scripts/Toolkit

workingDir=/data/zusers/moorej3/Moore-Lab/Papers/3D-Genome
list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
ccreDir=$registryDir/Cell-Type-Specific/Individual-Files
dataDir=/data/projects/encode/data/

cd $workingDir/5_Biological-Spectrum/

source ~/.bashrc
conda activate clinfun
mkdir -p Gini

for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    dnase=$(awk '{if (NR == '$i') print $9}' $list)
    h3k27ac=$(awk '{if (NR == '$i') print $12"-"$13}' $list)
    ctcfChIP=$(awk '{if (NR == '$i') print $14"-"$15}' $list)
    rm -f tmp.gini*

    awk '{if ($2 > 1.64) print $0}' $registryDir/signal-output/$ctcfChIP.txt | \
	awk 'FNR==NR {x[$1];next} ($4 in x)' - $registryDir/hg38-cCREs-Unfiltered.bed > tmp.ctcf

    cat $workingDir/2_Promoter-Centric/Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt \
	$workingDir/2_Promoter-Centric/Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | 
	grep CTCF-Anchors | awk 'FNR==NR {x[$5];next} ($1 in x)' tmp.ctcf - | \
	awk '{print $1 "\t" $NF}' | sort -u > tmp.x

    python ~/bin/count.py tmp.x 1 > tmp.counts

    for j in `seq 1 1 4`
    do
        awk '{if ($NF == "'$j'") print $1}' tmp.counts | \
            awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' - $gini | \
	    awk '{print $0 "\t" "'$j'"}' >> tmp.gini
	awk '{if ($NF == "'$j'") print $1}' tmp.counts | \
            awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' - $gini | sort -k2,2rg |
            awk '{a[i++]=$2;}END {x=int((i+1)/2); \
            if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}' | \
	    awk '{print "'$cell'" "\t" "'$j'" "\t" $1}' >> running.summary
    done

    awk '{if ($NF >= 5) print $1}' tmp.counts | \
        awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' - $gini | awk '{print $0 "\t" "5"}' >> tmp.gini

    awk '{if ($NF >= 5) print $1}' tmp.counts | \
        awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' - $gini | sort -k2,2rg |
        awk '{a[i++]=$2;}END {x=int((i+1)/2); \
        if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}' | \
        awk '{print "'$cell'" "\t" "5+" "\t" $1}' >> running.summary


    grep CTCF-Anchor $workingDir/2_Promoter-Centric/Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt | \
        awk '{print $NF}' > tmp.anchor-genes
    grep CTCF-Anchor $workingDir/2_Promoter-Centric/Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | \
        awk '{print $NF}' >> tmp.anchor-genes
    sort -u tmp.anchor-genes > tmp
    mv tmp tmp.anchor-genes

    cat $workingDir/2_Promoter-Centric/Promoter-Enhancer-Pairs/$cell.P-E-Pairs.txt \
        $workingDir/2_Promoter-Centric/Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt | \
        awk '{print $NF}' | awk 'FNR==NR {x[$1];next} !($1 in x)' tmp.anchor-genes - > tmp.other-genes


    awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' tmp.other-genes $gini | awk '{print $0 "\t" "0"}' >> tmp.gini



    awk -F "." 'FNR==NR {x[$1];next} ($1 in x)' tmp.other-genes $gini | sort -k2,2rg |
        awk '{a[i++]=$2;}END {x=int((i+1)/2); \
        if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}' | \
        awk '{print "'$cell'" "\t" "0" "\t" $1}' >> running.summary


    Rscript $scriptDir/jonckheere.R

    mv tmp.gini Gini/$cell-By-Count-Gini.txt
    rm tmp.*
done
mv running.summary Gini/Summary-Median-Gini-By-Counts.txt
