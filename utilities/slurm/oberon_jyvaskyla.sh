#!/usr/bin/bash
#SBATCH --job-name=lrd40
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=7G
#SBATCH --time=1-18:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@institute.edu

module load OpenMPI/3.1.3-GCC-8.2.0-2.31.1

export LC_NUMERIC="en_US.UTF-8"

suffix_arr=(A B C D E F G H I J K L M N O P Q R S T U V X Y Z a b c d e f g h i j k l m n o p q rr s t u v x y z)
suffix_plus=RR

topdir=$PWD
urqmd_dir=$topdir/RUN_Au_7.7_export
cg_dir=$topdir/cg_new_v.1.2.0
results=$cg_dir/results_Au_7_7
scratch_dir=$topdir/scratch

if [ ! -d $results ]; then
   mkdir -p $results
else
  cp -R $results $results\_backup_job_$SLURM_JOB_ID
fi

cd $urqmd_dir
make clean
make
cd $cg_dir
make clean
make 
cd $topdir

for kk in $(seq 60 99)
do
  cd $urqmd_dir

  rnds_base=$((70000000+$kk*100000))

  suffix=$suffix_plus${suffix_arr[$kk]}

  nev=500

  inc=$(($nev+2))
  for i in $(seq 1 40)
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
  cd $cg_dir

  odir=$scratch_dir/"out"$suffix
  current=partial_new
  old=partial_old

  if [ ! -d $current ]; then
     mkdir -p $current
  fi

  if [ ! -d $old ]; then
     mkdir -p $old
  else
    mv $results/* $old/
  fi

  if [ ! -d $odir ]; then
     mkdir -p $odir
  fi


  for i in $(seq 1 40)
  do
  ./cg_new.exe comp 0 $odir/out$i $scratch_dir/out$suffix\_$i.f14 timesteps.txt&
  done
  wait

  sleep 2

  mytask () {
    b=$(echo "$1/10" | bc -l)
    c=$(printf "%07.3f" $b)
    d=${c/./_}
    ./cg_new.exe avg 0 $current/tensor $odir/out*Tmunu_$d
    if [ -f "$old/tensor_Tmunu_$d" ]; then
       ./cg_new.exe avg 0 $results/tensor $current/tensor_Tmunu_$d $old/tensor_Tmunu_$d
       echo $d done 
    else
       mv $current/tensor_Tmunu_$d $results/
       echo "Moved $current/tensor_Tmunu_$d into $results/"
    fi
    #b=$(echo "$2/10" | bc -l)
    #c=$(printf "%07.3f" $b)
    #d=${c/./_}
    #./cg_new.exe avg 0 $current/tensor $odir/out*Tmunu_$d
    #if [ -f "$old/tensor_Tmunu_$d" ]; then
    #   ./cg_new.exe avg 0 $results/tensor $current/tensor_Tmunu_$d $old/tensor_Tmunu_$d
    #   echo $d done 
    #else
    #   mv $current/tensor_Tmunu_$d $results/
    #   echo "Moved $current/tensor_Tmunu_$d into $results/"
    #fi
  }

  #for i in $(seq 0 10 390)
  for i in $(seq 10 10 400)
  do
     #mytask $(($i+5))  $(($i+10))&
     mytask $i&
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

