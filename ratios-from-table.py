import sys
from astropy.table import Table
from muse_line_ratio import save_line_ratio_map

ratiotab = Table.read('line-ratio-list.tab', format='ascii.tab')

# Only do the lines, not continuum
prefix = 'linesum'

try:
    suffix = '-' + sys.argv[1]
except IndexError:
    suffix = ''

try:
    group = sys.argv[2]
except IndexError:
    group = None

for row in ratiotab:
    if group is None or group in row['Group']:
        save_line_ratio_map(row['Numerator'], row['Denominator'],
                            prefix, suffix)
