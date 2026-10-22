#Jill Moore
#UMass Chan
#October 2026
#
workingDir=~/Moore-Lab/Papers/3D-Genome

tfBed=$workingDir/0_Reference-Set-Up/Representative-TF-rPeak-Motifs.bed
tfList=$workingDir/0_Reference-Set-Up/TF-100-List.txt
ccres=$workingDir/4_Dual-Elements/Candidate-Dual-Elements.bed
bigWig=~/Lab/Reference/Human/hg38/Conservation/hg38_240mammal_phyloP.bw
scriptDir=~/GitHub/3D-Interactions/Analysis-Scripts/Toolkit
alignment=/zata/zippy/andrewsg/projects/zoonomia/hdf5/zoonomia_241_v2.h5

source ~/.bashrc
conda activate notebook

cd $workingDir/4_TF-Conservation

### Example 1

ccre=EH38E2226676
tf=JUN
name1=$(grep $ccre $workingDir/4_Dual-Elements/Candidate-Dual-Elements.bed | bedtools intersect -u -a Similarity-Scores/CTCF-$tf.75.tsv -b stdin | awk '{print $4}')
name2=$(grep $ccre $workingDir/4_Dual-Elements/Candidate-Dual-Elements.bed | bedtools intersect -u -a Similarity-Scores/$tf-CTCF.75.tsv -b stdin | awk '{print $4}')


python $scriptDir/plot_comparison_tree.py 241-mammalian-2020v2.phast-242.nh Similarity-Scores/CTCF-$tf.75.tsv Similarity-Scores/$tf-CTCF.75.tsv \
    $ccre-CTCF-$tf-Comparison.pdf --metadata-json species_metadata.json \
    --tf1-matrix Similarity-Scores/Matrix.CTCF-$tf.75.tsv --tf2-matrix Similarity-Scores/Matrix.$tf-CTCF.75.tsv \
    --report-output $ccre-CTCF-$tf-Comparison.tsv --tf1-region-name $name1 --tf2-region-name $name2 \
    --tf1-label CTCF --tf2-label $tf --metric sequence --color-group-ring

ccre=EH38E1475832
tf=CEBPB
name1=$(grep $ccre $workingDir/4_Dual-Elements/Candidate-Dual-Elements.bed | bedtools intersect -u -a Similarity-Scores/CTCF-$tf.75.tsv -b stdin | awk '{print $4}')
name2=$(grep $ccre $workingDir/4_Dual-Elements/Candidate-Dual-Elements.bed | bedtools intersect -u -a Similarity-Scores/$tf-CTCF.75.tsv -b stdin | awk '{print $4}')


python $scriptDir/plot_comparison_tree.py 241-mammalian-2020v2.phast-242.nh Similarity-Scores/CTCF-$tf.75.tsv Similarity-Scores/$tf-CTCF.75.tsv \
    $ccre-CTCF-$tf-Comparison.pdf --metadata-json species_metadata.json \
    --tf1-matrix Similarity-Scores/Matrix.CTCF-$tf.75.tsv --tf2-matrix Similarity-Scores/Matrix.$tf-CTCF.75.tsv \
    --report-output $ccre-CTCF-$tf-Comparison.tsv --tf1-region-name $name1 --tf2-region-name $name2 \
    --tf1-label CTCF --tf2-label $tf --metric sequence --color-group-ring


