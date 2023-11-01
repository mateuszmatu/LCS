#$ -S /bin/bash
#$ -l h_rt=10:00:00
#$ -q research-r8.q
#$ -l h_rss=8G
#$ -l mem_free=8G 
#$ -l h_data=8G
#$ -t 1-32
#$ -wd /lustre/storeB/users/mateuszm/scripts/logs


bash -l
conda activate opendrift
python /home/mateuszm/LCS/LCS/main.py $SGE_TASK_ID


