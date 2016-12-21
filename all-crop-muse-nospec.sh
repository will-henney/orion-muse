FILTERS="f487n f547m fq575n f656n f658n fq672n f673n fq674n"
for f in $FILTERS; do
    time python crop_muse.py $f nospec
done
