linelist=LineMaps/continuum-*[0-9][0-9][0-9][0-9].fits
for line in $linelist; do
    echo "Processing $line"
    time python multibin-map.py $line > ${line}-multibin.log
done
