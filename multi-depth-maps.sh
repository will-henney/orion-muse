source activate py27
for bin in 001 002 004 008 016 032 064 128 256; do
    for emline in S_II-6731 N_II-6583 S_III-9069 Cl_III-5538 O_II-7330 H_I-6563 O_III-5007; do
        time python correct-for-extinction.py $emline -bin$bin
        time python muse-depth.py $emline -bin$bin
    done
done
