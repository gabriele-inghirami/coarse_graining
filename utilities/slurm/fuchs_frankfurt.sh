#!/usr/bin/bash
#SBATCH --job-name=cgE8A
#SBATCH --output=sl_cgE8A_%j
#SBATCH --partition=fuchs
#SBATCH --account=xxxxx
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --time=19:30:00
#SBATCH --mail-type=ALL

module load comp/intel/2020.0

export LC_NUMERIC="en_US.UTF-8"

suffix_plus=E8A

topdir=$PWD
urqmd_dir=$topdir/RUN_$suffix_plus
scratch_dir=$topdir/TFO_scratch$suffix_plus
cg_dir=$topdir/cg_1.5.3b
results=$cg_dir/results_$suffix_plus
executable=cg.exe

if [ ! -d $cg_dir ]; then
   echo "Unrecoverable error, $cg_dir does not exists!!!"
   exit 1
fi

if [ ! -d $scratch_dir ]; then
   mkdir -p $scratch_dir
fi

if [ ! -d $results ]; then
   mkdir -p $results
else
  cp -R $results $results\_backup_job_$SLURM_JOB_ID
fi

#cd $urqmd_dir
#make clean
#make
#cd $cg_dir
#make clean
#make 
#cd $topdir

for kk in $(seq 1 19)
do
  echo Iteration $kk
  cd $urqmd_dir

  rnds_base=$((108159999+$kk*50)) #the factor multiplying $kk must be greater than the number of processors

  suffix=$suffix_plus$kk

  nev=200

#  for i in $(seq 0 39)
  for i in $(seq 0 19) #usually it is the number of processors available in the system
  do

    sed -e "s/nev 0/nev $nev/" inputfile > inputfile_tmp

    rnds=$(($rnds_base+$i))

    sed -e "s/rsd 0/rsd $rnds/" inputfile_tmp > inputfile$suffix\_$i

    rm inputfile_tmp

    sleep 1

    export ftn09=inputfile$suffix\_$i
    export ftn14=$scratch_dir/out$suffix\_$i.f14

    ./urqmd.$(uname -m) &

  done
  wait
  sleep 1
  echo "DISK_usage"
  du -hs $scratch_dir
  cd $cg_dir

  odir=$scratch_dir/"out"$suffix
  current=partial_new$suffix_plus
  old=partial_old$suffix_plus

  if [ ! -d $current ]; then
     mkdir -p $current
  fi

  if [ ! -d $old ]; then
     mkdir -p $old
  fi
  mv $results/* $old/

  if [ ! -d $odir ]; then
     mkdir -p $odir
  fi


  for i in $(seq 0 19)
  do
  ./$executable comp 0 $odir/out$i $scratch_dir/out$suffix\_$i.f14 timesteps.txt&
  done
  wait

#  for i in $(seq 10 19)
#  do
#  ./$executable comp 0 $odir/out$i $scratch_dir/out$suffix\_$i.f14 timesteps.txt&
#  done
#  wait

#  for i in $(seq 20 29)
#  do
#  ./$executable comp 0 $odir/out$i $scratch_dir/out$suffix\_$i.f14 timesteps.txt&
#  done
#  wait

#  for i in $(seq 30 39)
#  do
#  ./$executable comp 0 $odir/out$i $scratch_dir/out$suffix\_$i.f14 timesteps.txt&
#  done
#  wait
  sleep 2

  mytask () {
    num_steps=$2
    step_size=$3
    start_interval=$(($1+$step_size))
    end_interval=$(($1+$num_steps*$step_size))
    for xqx in $(seq $start_interval $step_size $end_interval)
    do
      b=$(echo "$xqx/100" | bc -l)
      c=$(printf "%07.3f" $b)
      d=${c/./_}
      ./$executable avg 0 $current/tensor $odir/out*Tmunu_$d
      if [ -f "$old/tensor_Tmunu_$d" ]; then
         ./$executable avg 0 $results/tensor $current/tensor_Tmunu_$d $old/tensor_Tmunu_$d
         echo $d done 
      else
         mv $current/tensor_Tmunu_$d $results/
         echo "Moved $current/tensor_Tmunu_$d into $results/"
      fi
    done
  }

  for i in $(seq 0 200 3800)
  do
     mytask $i 8 25&
  done
  wait
  sleep 1

  rm -rf $scratch_dir/out*
  rm -rf $old/* $current/*
  wait
  sleep 1

  cd $topdir
done
wait
sleep 2

