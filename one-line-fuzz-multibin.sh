linelist=LineMaps/linesum-$1-fuzz???.fits
for line in $linelist; do
    echo "Processing $line"
    time python multibin-map.py $line
done
