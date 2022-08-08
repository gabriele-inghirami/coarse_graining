for i in $(seq 1 3); do for k in data$i/*; do a=$(basename $k); ./cg.exe comp 0 test_d$i\f$a\_ $k/particles_binary.bin timesteps.txt ; done ; done
