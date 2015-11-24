for file in LineMaps/linesum*[0-9][0-9][0-9][0-9]$1-bin*.fits; do
    python muse-ew.py $file
done
