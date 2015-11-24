from __future__ import print_function
from astropy.table import Table
from misc_utils import sanitize_string
linetab = Table.read('basic-line-list.tab', format='ascii.tab')
for row in linetab:
    if row['strength'] < 0.0:
        wav = row['wav0']
        wavid = str(int(wav+0.5))
        species = sanitize_string(row['Ion'])
        print('ew-{}-{}'.format(species, wavid))
