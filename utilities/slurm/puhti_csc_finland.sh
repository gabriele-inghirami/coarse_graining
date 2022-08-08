#!/usr/bin/bash
#SBATCH --job-name=cgAu14_A2
#SBATCH --output=sl_cgAu14_A2_%j
#SBATCH --partition=small
#SBATCH --account=xxxxxxxx
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=1-03:59:59
#SBATCH --mem=100G
#SBATCH --gres=nvme:1600
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@institute.edu

module load intel/19.0.4

export LC_NUMERIC="en_US.UTF-8"

suffix_plus=14A

topdir=$PWD
urqmd_local=$topdir/RUN_$suffix_plus
urqmd_dir=$LOCAL_SCRATCH/RUN_$suffix_plus
cp -R $urqmd_local $urqmd_dir
sleep 2
wait
scratch_dir=$LOCAL_SCRATCH/TFO_scratch$suffix_plus
cg_dir=$topdir/cg_1.5.3
results_final=$cg_dir/results_Au_$suffix_plus
results=$LOCAL_SCRATCH/results_Au_$suffix_plus

if [ ! -d $cg_dir ]; then
   echo "Unrecoverable error, $cg_dir does not exists!!!"
   exit 1
fi

if [ ! -d $urqmd_dir ]; then
   echo "Unrecoverable error, $urqmd_dir does not exists!!!"
   exit 1
fi

if [ ! -d $scratch_dir ]; then
   mkdir -p $scratch_dir
   echo created $scratch_dir
fi

if [ ! -d $results ]; then
   mkdir -p $results
   echo created $results
fi

if [ -d $results_final ]; then
  mv $results_final $results_final\_backup_job_$SLURM_JOB_ID
  echo moved $results_final to $results_final\_backup_job_$SLURM_JOB_ID
fi

#cd $urqmd_dir
#make clean
#make
#cd $cg_dir
#make clean
#make 
#cd $topdir

#for kk in $(seq 100 599)
for kk in $(seq 10 29)
do
  echo Iteration $kk
  cd $urqmd_dir

  rnds_base=$((4010000+$kk*50))

  suffix=$suffix_plus$rnds_base$suffix_plus

  nev=75

  inc=$(($nev+1))
  for i in $(seq 0 39)
  do

    sed -e "s/nev 0/nev $nev/" inputfile > inputfile_tmp

    rnds=$(($rnds_base+$inc*$i))

    sed -e "s/rsd 0/rsd $rnds/" inputfile_tmp > inputfile$suffix\_$i

    rm inputfile_tmp

    sleep 1

    export ftn09=inputfile$suffix\_$i
    export ftn14=$scratch_dir/out$suffix\_$i.f14

    ./urqmd.$(uname -m) &

  done
  wait
  sleep 1
  echo DISK_SPACE_CHECK
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
  
  if [ "$(ls -A $results)" ]; then
     mv $results/* $old/
     echo moved the content of $results to $old
  fi

  if [ ! -d $odir ]; then
     mkdir -p $odir
  fi

  echo "entering 0-9 at" $(date)
  df
  for i in $(seq 0 9)
  do
  if [ -f $scratch_dir/out$suffix\_$i.f14 ]; then
     echo "found what needed " $i
     ./cg_new.exe comp 0 $odir/out$i $scratch_dir/out$suffix\_$i.f14 timesteps.txt&
  else
     echo $scratch_dir/out$suffix\_$i.f14 waiting
     wait 15
     if [ -f $scratch_dir/out$suffix\_$i.f14 ]; then
        ./cg_new.exe comp 0 $odir/out$i $scratch_dir/out$suffix\_$i.f14 timesteps.txt&
     else
        echo $scratch_dir/out$suffix\_$i.f14 definitely missing
     fi
  fi
  done
  wait
  df
  echo "exiting 0-9 at" $(date)

  echo "entering 10-19 at" $(date)
  df
  for i in $(seq 10 19)
  do
  if [ -f $scratch_dir/out$suffix\_$i.f14 ]; then
     echo "found what needed " $i
     ./cg_new.exe comp 0 $odir/out$i $scratch_dir/out$suffix\_$i.f14 timesteps.txt&
  else
     echo $scratch_dir/out$suffix\_$i.f14 waiting
     wait 15
     if [ -f $scratch_dir/out$suffix\_$i.f14 ]; then
        ./cg_new.exe comp 0 $odir/out$i $scratch_dir/out$suffix\_$i.f14 timesteps.txt&
     else
        echo $scratch_dir/out$suffix\_$i.f14 definitely missing
     fi
  fi
  done
  wait
  df
  echo "exiting 10-19 at" $(date)

  echo "entering 20-29 at" $(date)
  df
  for i in $(seq 20 29)
  do
  if [ -f $scratch_dir/out$suffix\_$i.f14 ]; then
     echo "found what needed " $i
     ./cg_new.exe comp 0 $odir/out$i $scratch_dir/out$suffix\_$i.f14 timesteps.txt&
  else
     echo $scratch_dir/out$suffix\_$i.f14 waiting
     wait 15
     if [ -f $scratch_dir/out$suffix\_$i.f14 ]; then
        ./cg_new.exe comp 0 $odir/out$i $scratch_dir/out$suffix\_$i.f14 timesteps.txt&
     else
        echo $scratch_dir/out$suffix\_$i.f14 definitely missing
     fi
  fi
  done
  wait
  df
  echo "exiting 20-29 at" $(date)

  echo "entering 30-39 at" $(date)
  df
  for i in $(seq 30 39)
  do
  if [ -f $scratch_dir/out$suffix\_$i.f14 ]; then
     echo "found what needed " $i
     ./cg_new.exe comp 0 $odir/out$i $scratch_dir/out$suffix\_$i.f14 timesteps.txt&
  else
     echo $scratch_dir/out$suffix\_$i.f14 waiting
     wait 15
     if [ -f $scratch_dir/out$suffix\_$i.f14 ]; then
        ./cg_new.exe comp 0 $odir/out$i $scratch_dir/out$suffix\_$i.f14 timesteps.txt&
     else
        echo $scratch_dir/out$suffix\_$i.f14 definitely missing
     fi
  fi
  done
  wait
  df
  echo "exiting 30-39 at" $(date)

  sleep 1

  mytask () {
    num_steps=$2
    step_size=$3
    start_interval=$(($1+$step_size))
    end_interval=$(($1+$num_steps*$step_size))
    echo "This is mytask from $start_interval to $end_interval with step size"
    for xqx in $(seq $start_interval $step_size $end_interval)
    do
      b=$(echo "$xqx/100" | bc -l)
      c=$(printf "%07.3f" $b)
      d=${c/./_} 
      
      ./cg_new.exe avg 0 $current/tensor $odir/out*Tmunu_$d
      sleep 1
      if [ -f $current/tensor_Tmunu_$d ]; then
         echo "found in mytask $current/tensor_Tmunu_$d"
         ls -l $current/tensor_Tmunu_$d
      else
         echo "waiting in mytask $current/tensor_Tmunu_$d"
         sleep 15
         if [ -f $current/tensor_Tmunu_$d ]; then
            echo "appeared in mytask $current/tensor_Tmunu_$d"
            ls -l $current/tensor_Tmunu_$d
         else
            echo "missing in mytask $current/tensor_Tmunu_$d. FATAL ERROR."
            exit 3
         fi
      fi
      if [ -f "$old/tensor_Tmunu_$d" ]; then
          echo "averaging in mytask with older results"
         ./cg_new.exe avg 0 $results/tensor $current/tensor_Tmunu_$d $old/tensor_Tmunu_$d
         echo $d done 
      else
         #we already checked that this file exists 
         echo "Moving $current/tensor_Tmunu_$d into $results/"
         mv $current/tensor_Tmunu_$d $results/
         echo "Moved $current/tensor_Tmunu_$d into $results/"
      fi

    done
  }

  for i in $(seq 0 400 3600)
  do
     mytask $i 16 25&
  done
  wait
  sleep 1

  rm -rf $scratch_dir/out*
  rm -rf $old/* $current/*
  wait

  echo "cycle ended"
  date
  df
  sleep 2
  cd $topdir
done
wait

#we check that we obtained all the results
echo "*** FINAL CHECK ***"
for i in $(seq 25 25 4000)
do
    b=$(echo "$i/100" | bc -l)
    c=$(printf "%07.3f" $b)
    d=${c/./_} 
    if [ -f $results/tensor_Tmunu_$d ]; then
       ls -l $results/tensor_Tmunu_$d
    else 
       echo "waiting $results/tensor_Tmunu_$d"
       sleep 15
       if [ -f $results/tensor_Tmunu_$d ]; then
          ls -l $results/tensor_Tmunu_$d
       else
          echo "$results/tensor_Tmunu_$d probably definitely missing"
       fi
    fi
done

cp -R $results $results_final && rm -rf $results
wait   
ls -l $results_final
wait
sleep 2

