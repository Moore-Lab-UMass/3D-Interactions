
workingDir=~/Moore-Lab/Papers/3D-Genome

tfBed=$workingDir/0_Reference-Set-Up/Representative-TF-rPeak-Motifs.bed
tfList=$workingDir/0_Reference-Set-Up/TF-100-List.txt
ccres=$workingDir/3_Dual-Elements/Candidate-Dual-Elements.bed
bigWig=~/Lab/Reference/Human/hg38/Conservation/hg38_240mammal_phyloP.bw

cd $workingDir/4_Conservation

for j in `seq 1 1 42`
do
    echo $j
    tf=$(awk -F "\t" '{if (NR == '$j') print $1}' $tfList)

    awk -F "\t" '{if ($4 ~ /'$tf'_/) print $0}' $tfBed | \
        bedtools intersect -u -b $ccres -a stdin | \
        awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > test.$tf

    computeMatrix reference-point --binSize 1 --referencePoint center -b 25 -a 25 \
	-R test.$tf -S $bigWig -o $tf.gz --outFileSortedRegions $tf.bed
    plotProfile -m $tf.gz -out $tf.png
    ~/bin/bigWigAverageOverBed $bigWig test.$tf out.tab
    awk '{print $1 "\t" $5 "\t" "'$tf'"}' out.tab >> running.output
done

