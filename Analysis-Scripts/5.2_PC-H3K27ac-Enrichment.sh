



scriptDir=/data/zusers/moorej3/GitHub/3D-Interactions/Analysis-Scripts/Toolkit
control=/data/projects/encode/Registry/V4/GRCh38/GRCh38-cCREs-All.bed

workingDir=/data/zusers/moorej3/Moore-Lab/Papers/3D-Genome
list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38/
dataDir=/data/projects/encode/data/

cd $workingDir/5_Biological-Spectrum/H3K27ac-Enrichment


cat Raw-Runs/K562/*txt | awk 'BEGIN{print "Biosample"}{print $1}' > tmp.matrix

for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')
    echo $cell
#    group=$workingDir/4_Dual-Elements/H3K27ac-Enrichment/cCRE-Sets/$cell-Candidates.bed
#    outputDir=$workingDir/4_Dual-Elements/H3K27ac-Enrichment/Raw-Runs/$cell
#    pc=$workingDir/2_Promoter-Centric/Promoter-CTCF-Pairs/$cell.P-C-Pairs.txt

#    grep CTCF-Anchor $pc | awk 'FNR==NR {x[$1];next} ($5 in x)' - $control > $group

#    num=562
#    jobid=$(sbatch --nodes 1 --array=1-$num%25 --mem=10G --time=00:30:00 \
#        --output=/zata/zippy/moorej3/Job-Logs/jobid_%A_%a.output \
#        --error=/zata/zippy/moorej3/Job-Logs/jobid_%A_%a.error \
#        $scriptDir/Batch-Run-H3K27ac-Enrichment.sh $group $control $outputDir | awk '{print $4}')
#    echo $jobid
#    sleep 20
#    jobList=100
#    while [ $jobList -gt 1 ]
#    do
#        jobList=$(squeue -j $jobid | wc -l | awk '{print $1}')
#        sleep 10
#    done

    cat Raw-Runs/$cell/*txt | awk 'BEGIN{print "'$cell'"}{print ($2/$3)/($4/$5)}' > tmp.col
    paste tmp.matrix tmp.col > tmp.tmp
    mv tmp.tmp tmp.matrix

done

mv tmp.matrix Matrix-H3K27ac-Enrichment.txt

awk '{if (NR > 1) print $1 "\t" ($2 + $4 + $5 + $6 + $7)/5}' Matrix-H3K27ac-Enrichment.txt | grep -v HepG2 | grep -v K562 | grep -v Caco-2 | grep -v MCF > tmp.ave

