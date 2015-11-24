from astropy.table import Table
from misc_utils import sanitize_string
from muse_ew_utils import save_ew
linetab = Table.read('basic-line-list.tab', format='ascii.tab')
for row in linetab:
    wav = row['wav0']
    wavid = str(int(wav+0.5))
    species = sanitize_string(row['Ion'])
    fname = 'LineMaps/linesum-{}-{}.fits'.format(species, wavid)
    save_ew(fname)
