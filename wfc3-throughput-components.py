import matplotlib.pyplot as plt
import pysynphot
bp = pysynphot.ObsBandpass('wfc3,uvis1,f547m')
comps = {}
for p in bp.obsmode._throughput_filenames: 
    if p != 'clear':
        if 'uvis_ccd1' in p:
            filterfile_ccd2 = p.replace('ccd1', 'ccd2')
        k = p.split('/')[-1].split('.')[0]
        comps[k] = pysynphot.FileBandpass(p)
for k, v in comps.items():
    plt.plot(v.wave, v.throughput, label=k)

plt.plot(bp.wave, bp.throughput, label='F469N')
bp2 = pysynphot.FileBandpass(filterfile_ccd2)
plt.plot(bp2.wave, bp2.throughput, label='CCD 2')
plt.xlim(4500, 8000)
plt.legend(fontsize='small')
plt.savefig('wfc3-throughput-components.pdf')
