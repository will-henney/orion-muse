SN=$1
files=LineMaps/ratio-[0-9][0-9][0-9][0-9]-[0-9][0-9][0-9][0-9]-bin001.fits
for file in $files; do
    prefix=$(basename $file -bin001.fits)
    if [ -f LineMaps/${prefix}-SN-bin256.fits ]; then
        echo "+++ Calculating S/N masks for $prefix at S/N = $SN"
        python multibin-mask-s-n.py $prefix $SN
        echo "+++ Calculating combined image for $prefix at S/N = $SN"
        python multibin-combine-s-n.py $prefix $SN
        echo "+++"
    fi
done
