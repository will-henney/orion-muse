linelist=LineMaps/linesum-*[0-9][0-9][0-9][0-9]-fuzz$1.fits
for line in $linelist; do
    echo "Processing $line"
    time python multibin-map.py $line > ${line}-multibin.log 2>&1
done
