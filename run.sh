#!/usr/bin/bash
#SBATCH --job-name=run
#SBATCH --output=sl_run_%j
#SBATCH --partition=test
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --time=1:59:59
#SBATCH --mail-type=ALL

module load comp/intel/2020.0

export LC_NUMERIC="en_US.UTF-8"

fun ()
{
for k in data$1/[0-9]*; do a=$(basename $k); echo "ITERATION $a"; ./cg.exe comp 0 test_d$1\f$a\_ $k/particles_binary.bin timesteps.txt ; sleep 2; echo "END ITERATION $a"; done; 
}

for m in $(seq 1 6)
do
echo "MMMAIN Cyle" $m
fun $m &
done

wait

for m in $(seq 7 12)
do
echo "MMMAIN Cyle" $m
fun $m &
done

wait

for m in $(seq 13 18)
do
echo "MMMAIN Cyle" $m
fun $m &
done

wait

sleep 10
