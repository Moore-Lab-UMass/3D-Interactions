#Jill Moore
#UMass Chan
#July 2026

workingDir=~/Moore-Lab/Papers/3D-Genome/
list=$workingDir/0_Reference-Set-Up/Master-Biosample-Experiment-List.txt
registryDir=/home/moorej3/Lab/ENCODE/Encyclopedia/V7/Registry/V7-hg38
ccreDir=$registryDir/Cell-Type-Specific/Individual-Files
chiaDir=/data/projects/encode/3D-Chromatin/ChIA-PET/Ruan-ENCODE4
hicDir=$registryDir/3D-Links/HiC
dataDir=/data/projects/encode/data/
folder=Assay-Overlap-Details

cd $workingDir/1_General-Comparisons
mkdir -p $folder

source ~/.bashrc

for i in `seq 1 1 6`
do
    cell=$(awk '{if (NR == '$i') print $1}' $list | awk -F "_" '{print $3}')

    rnap2=$(awk '{if (NR == '$i') print $2}' $list)
    ctcf=$(awk '{if (NR == '$i') print $4}' $list)
    hic=$(awk -F "\t" '{if (NR == '$i') print $19}' $list)

    RNAPII=$chiaDir/$rnap2.DAC-Links.V2.bedpe
    CTCF=$chiaDir/$ctcf.DAC-Links.V2.bedpe
    HiC=$hicDir/$hic.bedpe

    allLinks=(RNAPII CTCF HiC)
    for links1 in ${allLinks[@]}
    do
	awk '{if ($1 !~ /#/) print $0}' "${!links1}" > tmp.1
	total=$(wc -l tmp.1 | awk '{print $1}')
        #echo $links1
        for links2 in ${allLinks[@]}
    	do
	    #echo -e "\t" $links2
	    awk '{if ($1 !~ /#/) print $0}' "${!links2}" > tmp.2
	    bedtools pairtopair -is -type both -a tmp.1 -b tmp.2 > tmp.$links2
	done
        a=$(awk 'FNR==NR {x[$1$2$3$4$5$6];next} ($1$2$3$4$5$6 in x)' tmp.RNAPII tmp.CTCF |
		awk 'FNR==NR {x[$1$2$3$4$5$6];next} ($1$2$3$4$5$6 in x)' tmp.HiC - | 
		awk '{print $1$2$3$4$5$6}' | sort -u | wc -l | awk '{print $1}')
	b=$(awk 'FNR==NR {x[$1$2$3$4$5$6];next} ($1$2$3$4$5$6 in x)' tmp.RNAPII tmp.CTCF |
                awk '{print $1$2$3$4$5$6}' | sort -u | wc -l | awk '{print $1 - '$a'}')
	c=$(awk 'FNR==NR {x[$1$2$3$4$5$6];next} ($1$2$3$4$5$6 in x)' tmp.RNAPII tmp.HiC |
                awk '{print $1$2$3$4$5$6}' | sort -u | wc -l | awk '{print $1 - '$a'}')
	d=$(awk 'FNR==NR {x[$1$2$3$4$5$6];next} ($1$2$3$4$5$6 in x)' tmp.HiC tmp.CTCF |
                awk '{print $1$2$3$4$5$6}' | sort -u | wc -l | awk '{print $1 - '$a'}')
	e=$(wc -l tmp.1 | awk '{print $1 - ('$a' + '$b' + '$c' + '$d')}')

	echo -e $cell "\t" $links1 "\t" $total "\t" $a "\t" $b "\t" $c "\t" $d "\t" $e >> running.output

	awk 'FNR==NR {x[$1$2$3$4$5$6];next} ($1$2$3$4$5$6 in x)' tmp.RNAPII tmp.CTCF |
                awk 'FNR==NR {x[$1$2$3$4$5$6];next} ($1$2$3$4$5$6 in x)' tmp.HiC - |
                awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' | sort -u | 
		awk '{print $0 "\t" "ALL"}' > out.$links1
	awk 'FNR==NR {x[$1$2$3$4$5$6];next} ($1$2$3$4$5$6 in x)' tmp.RNAPII tmp.CTCF |
                awk 'FNR==NR {x[$1$2$3$4$5$6];next} !($1$2$3$4$5$6 in x)' tmp.HiC - |
                awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' | sort -u |
                awk '{print $0 "\t" "RNAPII-CTCF"}' >> out.$links1
	awk 'FNR==NR {x[$1$2$3$4$5$6];next} ($1$2$3$4$5$6 in x)' tmp.RNAPII tmp.HiC |
                awk 'FNR==NR {x[$1$2$3$4$5$6];next} !($1$2$3$4$5$6 in x)' tmp.CTCF - |
                awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' | sort -u |
                awk '{print $0 "\t" "RNAPII-HiC"}' >> out.$links1
	awk 'FNR==NR {x[$1$2$3$4$5$6];next} ($1$2$3$4$5$6 in x)' tmp.CTCF tmp.HiC |
                awk 'FNR==NR {x[$1$2$3$4$5$6];next} !($1$2$3$4$5$6 in x)' tmp.RNAPII - |
                awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' | sort -u |
                awk '{print $0 "\t" "CTCF-HiC"}' >> out.$links1
	awk 'FNR==NR {x[$1$2$3$4$5$6];next} !($1$2$3$4$5$6 in x)' out.$links1 tmp.1 |
		awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' | sort -u |
                awk '{print $0 "\t" "Single"}' > tmp.tmp
	cat tmp.tmp >> out.$links1
	sort -k1,1 -k2,2n -k3,3n -k5,5n -k6,6n out.$links1 > $folder/$cell.$links1.txt
	rm out.$links1
    done
done

rm tmp.*

mv running.output Summary-Assay-Overlap.txt
