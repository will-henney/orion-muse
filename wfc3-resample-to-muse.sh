F=$1
MDIR=~/Source/Montage/bin
TDIR=~/Work/RubinWFC3/Tsquared
$MDIR/mProjectPP -h 0 -X $TDIR/newsmooth-${F}-s080.fits wfc3-new-resample-muse-$F.fits muse-full-frame.hdr
