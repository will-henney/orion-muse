linelist=LineMaps/continuum-*-[0-9][0-9][0-9][0-9].fits
D=../multibin-maps
for line in $linelist; do
    echo "Processing $line"
    time python $D/multibin-map.py $line > ${line}-multibin.log 2>&1
done
