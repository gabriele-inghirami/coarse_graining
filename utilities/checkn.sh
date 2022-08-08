export LC_NUMERIC="en_US.UTF-8"
#this program checks that all the processed files in <directory> are based on the same number of events
#this program exploits the auxiliary program check_nev.py

#usage: bash checkn.sh directory

#here we assume that the output files are called tensor_Tmunu_001_000... and so on
#the first number in the for k loop is the starting time mulitplied by 100
#the second number is the width of the timestep*100
#the third number is the final timestep*100

for k in $(seq 25 25 4000)
do
    b=$(echo "$k/100" | bc -l)
    c=$(printf "%07.3f" $b)
    d=${c/./_}
    if [ ! -e $1/tensor_Tmunu_$d ]
    then
    echo $1/tensor_Tmunu_$d missing
    exit
    fi 
done
echo "All files are there. Now checking the number of events in each of them"
python3 check_nev.py $1/tensor_Tmunu_0*
