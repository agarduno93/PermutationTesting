#! /bin/bash
#SBATCH -J 3_Multivariate_Analysis-Restrict4-Primary-RandomShuffleNDcases_RUNNER
#SBATCH -o 3_Multivariate_Analysis-Restrict4-Primary-RandomShuffleNDcases_RUNNER.out
#SBATCH --mem=10G


rand=$RANDOM
temp_dir="temp_alexis_runs_NDcases_"$rand""
echo "Temporary Directory: "$temp_dir

# Check if temp directory is creatable (See if there's a directory of that name)
if [ -d $temp_dir ]
then
    echo "Temp directory already exists, rerun to generate a new random number"
    exit
else
    echo "Temp directory valid"
    mkdir $temp_dir
fi


echo "#! /bin/bash
#SBATCH --mem=4G
#SBATCH --array=1-1000%1
#SBATCH -o "$temp_dir"/RandomShuffleNDcases_run_%A.%a.out
#SBATCH -J RandomShuffleNDcases_"$rand"

date

Rscript 3_Multivariate_Analysis-Restrict4-Primary-RandomShuffleNDcases.R \$SLURM_ARRAY_TASK_ID $temp_dir

date" > $temp_dir/random_shuffle_RUNNER.sh

sbatch $temp_dir/random_shuffle_RUNNER.sh



# Wait for random shuffling to finish running
echo "Waiting for analysis to finish running..."
c=1
while [ $c -ne 0 ]
    do
    c=$(squeue -o "%.100i %.100j" | awk '{print $2}' | grep "RandomShuffleNDcases_"$rand"" | wc -l)
    sleep 60
done



# Merge all outputs together
head -1 $temp_dir/All_10May21_NDShuffle_it1.txt > $temp_dir/All_10May21_NDShuffle.txt
ls $temp_dir/All_10May21_NDShuffle_it* | sed 's/^/sed 1d /g' | sed "s]$] >> $temp_dir/All_10May21_NDShuffle.txt]g" > $temp_dir/All_10May21_NDShuffle_combine.sh
bash $temp_dir/All_10May21_NDShuffle_combine.sh 

mv $temp_dir/All_10May21_NDShuffle.txt All_10May21_NDShuffle.txt

# Merge all outputs together
head -1 $temp_dir/ND_10May21_NDShuffle_it1.txt > $temp_dir/ND_10May21_NDShuffle.txt
ls $temp_dir/ND_10May21_NDShuffle_it* | sed 's/^/sed 1d /g' | sed "s]$] >> $temp_dir/ND_10May21_NDShuffle.txt]g" > $temp_dir/ND_10May21_NDShuffle_combine.sh
bash $temp_dir/ND_10May21_NDShuffle_combine.sh 

mv $temp_dir/ND_10May21_NDShuffle.txt ND_10May21_NDShuffle.txt

# Merge all outputs together
head -1 $temp_dir/T2D_10May21_NDShuffle_it1.txt > $temp_dir/T2D_10May21_NDShuffle.txt
ls $temp_dir/T2D_10May21_NDShuffle_it* | sed 's/^/sed 1d /g' | sed "s]$] >> $temp_dir/T2D_10May21_NDShuffle.txt]g" > $temp_dir/T2D_10May21_NDShuffle_combine.sh
bash $temp_dir/T2D_10May21_NDShuffle_combine.sh 

mv $temp_dir/T2D_10May21_NDShuffle.txt T2D_10May21_NDShuffle.txt


# Remove temporary directory
#rm -r $temp_dir



