for bin in 001 002 004 008 016 032 064 128 256; do
#for bin in 008 016 032 064 128 256; do
    time python muse-Te-Ne-histograms.py $1 $bin
done
