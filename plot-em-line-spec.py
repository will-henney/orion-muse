from astropy.table import Table
from misc_utils import sanitize_string
from matplotlib import pyplot as plt
import seaborn as sns

linetab = Table.read('basic-line-list.tab', format='ascii.tab')

for row in linetab:
    wav = row['wav0']
    wavid = str(int(wav+0.5))
    species = sanitize_string(row['Ion'])
    sname = 'Linemaps/spec1d-{}-{}.tab'.format(species, wavid)
    spec = Table.read(sname, format='ascii.tab')
    for xkey, xlabel in [['vhel', 'Heliocentric velocity, km/s'],
                         ['wav', 'Observed air wavelength, Angstrom']]:
        fig, ax = plt.subplots(1, 1)
        ax.plot(spec[xkey], spec['flux'])
        ax.plot(spec[xkey], spec['cont'])
        if xkey == 'wav':
            ax.set_xlim(row['wav0']-8.0, row['wav0']+8.0)
        else:
            ax.set_xlim(-300.0, 300.0)
        ax.set_ylim(0.0, None)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Mean flux per pixel')
        ax.set_title('{} {:.2f}'.format(row['Ion'], row['wav0']))
        fig.set_size_inches(5, 5)
        fig.savefig(sname.replace('.tab', '-{}.pdf'.format(xkey)))
        del(fig)
        del(ax)
