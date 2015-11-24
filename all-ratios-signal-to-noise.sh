files=LineMaps/ratio-[0-9][0-9][0-9][0-9]-[0-9][0-9][0-9][0-9]-bin001.fits
for file in $files; do
    prefix=$(basename $file -bin001.fits)
    echo "Calculating S/N for $prefix"
    python multibin-signal-to-noise.py $prefix
done
