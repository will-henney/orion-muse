# for bin in 001 ; do
#     for fuzz in 009; do
for bin in 001 002 004 008 016 032 064 128 256; do
    for fuzz in 000 001 002 003 004 005 006 007 008 009; do
        python muse-Te-Ne-fuzz-histograms.py $fuzz $bin full
        python muse-Te-Ne-fuzz-histograms.py $fuzz $bin sweet
    done
done
