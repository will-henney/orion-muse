shopt -s extglob
files=LineMaps/linesum-+([^-])-[0-9][0-9][0-9][0-9].fits
for file in $files; do
    prefix=$(basename $file .fits)
    echo "Calculating S/N for $prefix"
    python multibin-signal-to-noise.py $prefix
done
