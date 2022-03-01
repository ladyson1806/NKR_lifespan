#!/bin/bash
#SBATCH --account=rrg-hussinju
#SBATCH --time=4:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-500%50
#SBATCH --output=./tmp/%x_%A_%a_out.txt
#SBATCH --error=./tmp/%x_%A_%a_err.txt

DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" HG_list)
cd $DIR

SEQFILE=$(find . -type f -name *.seq)
/home/sbesse/projects/def-hussinju/shared/bin/Tango/tango_x86_64_release -inputfile=$SEQFILE > log.out
cp *_aggregation.txt /home/sbesse/scratch/computational_mutagenesis/RESULTS/"${SEQFILE//.\/}"_aggregation.txt
rm *.txt


#!/bin/bash
#SBATCH --account=rrg-hussinju
#SBATCH --time=4:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-500%50
#SBATCH --output=./tmp/%x_%A_%a_out.txt
#SBATCH --error=./tmp/%x_%A_%a_err.txt

DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" MM_list)
cd $DIR

SEQFILE=$(find . -type f -name *.seq)
/home/sbesse/projects/def-hussinju/shared/bin/Tango/tango_x86_64_release -inputfile=$SEQFILE > log.out
cp *_aggregation.txt /home/sbesse/scratch/computational_mutagenesis/RESULTS/"${SEQFILE//.\/}"_aggregation.txt
rm *.txt


#!/bin/bash
#SBATCH --account=rrg-hussinju
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=1
#SBATCH --array=701-773%6
#SBATCH --output=./tmp/%a_%A_out.txt
#SBATCH --error=./tmp/%a_%A_err.txt


DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" rerun_list)
cd $DIR

SEQFILE=$(find . -type f -name *.seq)
echo "${SEQFILE}"
echo "${SEQFILE//.\/}"
/home/sbesse/projects/def-hussinju/shared/bin/Tango/tango_x86_64_release -inputfile=$SEQFILE > log.out
cp *_aggregation.txt* /home/sbesse/scratch/computational_mutagenesis/RESULTS_RERUN/"${SEQFILE//.\/}"_aggregation.txt
rm *.txt
