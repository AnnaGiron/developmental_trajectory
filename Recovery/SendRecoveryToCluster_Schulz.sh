#here you need to specify how many instances you want to run simulations on. 
#currently it si set up as such that for any "real" collected data from our pilot there will be Simulations.
#this is not to fit a model yet. this makes extensive Simulations. I fit the model back at another place.
#HowManyTimes=3
HowManySubjects=100
HowManyAqus=3
HowManyLearners=2

#j=2;
#i=1;
#k=5;
#NJobs=100
#j=2;
#i=1;
#k=5;
PATH_LOG_E="./RecoverylogsE/$(date '+%Y%m%d_%H%M%S')/slurm-%j.err"
PATH_LOG_O="./RecoverylogsY/$(date '+%Y%m%d_%H%M%S')/slurm-%j.out"
# path to the text file with all subject ids:
#PATH_SUB_LIST="${PATH_SCRIPT}/subs.txt"
# CREATE RELEVANT DIRECTORIES:
# ==============================================================================
# create output directory:
if [ ! -d ${PATH_LOG_E} ]; then
	mkdir -p ${PATH_LOG_E}
fi
# create directory for log files:
if [ ! -d ${PATH_LOG_O} ]; then
	mkdir -p ${PATH_LOG_O}
fi

for z in `seq 1 $HowManyTimes`;
do
  for i in `seq 100 $HowManySubjects`;
  do
  	for j in `seq 1 $HowManyAqus`;
  	do
  		for k in `seq 1 $HowManyLearners`;
  	   do
  		echo "$i  $Subjects";
  		#here i specify the job for the cluster.
  		#for input_file in INPUT/* ; do
  			#echo "#PBS -m n"                         > job.pbs
  		  echo '#!/bin/bash'                                > job.slurm
       		  echo "#SBATCH --job-name S_Sims_GPUCB_CV_Recov${i}_${j}_${k}_${idx}"     >> job.slurm
        	  echo "#SBATCH --partition long"                   >> job.slurm
        	  echo "#SBATCH --mem 10GB"                         >> job.slurm
        	  echo "#SBATCH --cpus-per-task 8"                  >> job.slurm
        	  echo "#SBATCH --time 100:0:0"                     >> job.slurm
        	  echo "#SBATCH --workdir ."                        >> job.slurm
        	  echo "#SBATCH --error ${PATH_LOG_E}"              >> job.slurm
        	  echo "#SBATCH --output ${PATH_LOG_O}"             >> job.slurm
  		  echo "module load R/4; Rscript Model_Recovery_Cluster_Schulz.R $i $j $k $z"	>>  job.slurm
      		  sbatch job.slurm
      		  rm -f job.slurm
  		  done
  		done
  	done
done
