F=$1
MDIR=~/Source/Montage/bin
TDIR=~/Work/RubinWFC3/Tsquared
$MDIR/mProjectPP -h 0 -X $TDIR/full_${F}-s070.fits wfc3-resample-muse-$F.fits muse-full-frame.hdr
