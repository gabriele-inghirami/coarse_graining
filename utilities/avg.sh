export LC_NUMERIC="en_US.UTF-8"
for k in $(seq 500 500 2500)
do
    b=$(echo "$k/100" | bc -l)
    c=$(printf "%07.3f" $b)
    d=${c/./_}
    for k in $(seq 1 18)
    do
      ./cg.exe avg 1 res_smash/piece_$k\_tensor rest/test_d$k\f*_Tmunu_$d
    done
done

