
workingDir=~/Moore-Lab/Papers/3D-Genome

tfBed=$workingDir/0_Reference-Set-Up/Representative-TF-rPeak-Motifs.bed
tfList=$workingDir/0_Reference-Set-Up/TF-100-List.txt
ccres=$workingDir/3_Dual-Elements/Candidate-Dual-Elements.bed
bigWig=~/Lab/Reference/Human/hg38/Conservation/hg38_240mammal_phyloP.bw
scriptDir=~/GitHub/3D-Interactions/Analysis-Scripts/Toolkit
alignment=/zata/zippy/andrewsg/projects/zoonomia/hdf5/zoonomia_241_v2.h5

source ~/.bashrc
conda activate notebook

cd $workingDir/4_TF-Conservation

rm -f running.output

bedtools intersect -wo -a Beds/test.CTCF -b $workingDir/3_Dual-Elements/Candidate-Dual-Elements.bed | awk '!seen[$10]++'> tmp.CTCF

for j in `seq 1 1 42`
do
    tf=$(awk -F "\t" '{if (NR == '$j') print $1}' $tfList)

    bedtools intersect -wo -a Beds/test.$tf -b $workingDir/3_Dual-Elements/Candidate-Dual-Elements.bed | awk '!seen[$10]++' > tmp.$tf

    awk 'FNR==NR {x[$10];next} ($10 in x)' tmp.CTCF tmp.$tf | sort -k10,10 > match.tf
    awk 'FNR==NR {x[$10];next} ($10 in x)' tmp.$tf tmp.CTCF | sort -k10,10 > match.ctcf

    python $scriptDir/count_species_alignment.py $alignment match.tf tf.tsv --mode both --min-fraction=0.75 \
	--include-species-indices --metadata-json species_metadata.json --reference-species Homo_sapiens \
	--matrix-output tf_matrix.tsv --matrix-metric sequence
    python $scriptDir/count_species_alignment.py $alignment match.ctcf ctcf.tsv --mode both --min-fraction=0.75 \
	--include-species-indices --metadata-json species_metadata.json --reference-species Homo_sapiens \
	--matrix-output ctcf_matrix.tsv --matrix-metric sequence


    awk '{print $12}' tf.tsv > tmp.b
    awk '{print $12}' ctcf.tsv > tmp.a

    a=$(sort -k1,1rg tmp.a | awk '{a[i++]=$1;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')
    b=$(sort -k1,1rg tmp.b | awk '{a[i++]=$1;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}')

    p=$(Rscript $scriptDir/paired_wilcox_test.R | awk '{print $1}')

    c=$(paste tmp.a tmp.b | awk 'BEGIN{sum=0}{if ($1 < $2) sum+=1}END{print sum/NR}')

    paste tmp.a tmp.b | awk '{print $1-$2}' | sort -k1,1rg | \
	awk '{a[i++]=$1;}END {x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}' | \
	awk '{print "CTCF" "\t" "'$tf'" "\t" $1 "\t", "'$a'" "\t" "'$b'" "\t" "'$c'" "\t" "'$p'"}' >> running.output

    mv tf.tsv Similarity-Scores/$tf-CTCF.75.tsv
    mv ctcf.tsv Similarity-Scores/CTCF-$tf.75.tsv

    mv tf_matrix.tsv Similarity-Scores/Matrix.$tf-CTCF.75.tsv
    mv ctcf_matrix.tsv Similarity-Scores/Matrix.CTCF-$tf.75.tsv

done


mv running.output Summary-Paired-Similarity-Scores.txt
