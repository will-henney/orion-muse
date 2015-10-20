for i in $(seq 7); do
    python rebin_datacube.py muse-hr-data-wavsec$i.fits 5 5
done
