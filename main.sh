#$ -S /bin/bash
#$ -l h_rt=10:00:00
#$ -q research-r8.q
#$ -l h_rss=8G
#$ -l mem_free=8G
#$ -t 1-31
#$ -wd /lustre/storeB/users/mateuszm/scripts/logs

source /modules/bionic/conda/Feb2021/etc/profile.d/conda.sh
bash -l
conda activate opendrift
python /lustre/storeB/users/mateuszm/scripts/main.py $SGE_TASK_ID


