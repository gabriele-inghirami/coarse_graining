outdir=res_77
mkdir -p $outdir
export LC_NUMERIC="en_US.UTF-8"
for k in $(seq 25 25 3000)
do
    b=$(echo "$k/100" | bc -l)
    c=$(printf "%07.3f" $b)
    d=${c/./_}
    ./cg.exe avg 1 $outdir/tensor results_Au_*/tensor_Tmunu_$d
done

