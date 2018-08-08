

import matplotlib.pyplot
import sys
import os



plot_name = sys.argv[1]
sim_names = sys.argv[2:]
#sim_colors_labels = dict({
#    's3-c1': ['#90EE90','S3 + 30S-S3-2i2p'],
#    's3-c2': ['#87CEFA','S3 + 30S-S3S5-2i2p'],
#    's5-c2': ['#FFDEAD','S5 + 30S-S3S5-2i2p'],
#    's5-c3': ['#FF6347','S5 + 30S-S5-2i2p'],
#    's3-c1-4': ['#3CB371','S3 + 30S-S3-4adv'],
#    's3-c2-4': ['#6495ED','S3 + 30S-S3S5-4adv'],
#    's5-c2-4': ['#D2691E','S5 + 30S-S3S5-4adv'],
#    's5-c3-4': ['#DC143C','S5 + 30S-S5-4adv']
#})

#python plot_rates.py s3-same s3s5-4adv.s3-4adv s3-4adv.s3-4adv s3-2i2p.s3-2i2p s3s5-2i2p.s3-2i2p 

#python plot_rates.py s3-md s3-2i2p.s3-md s3s5-2i2p.s3-md s3s5-4adv.s3-md s3-4adv.s3-md

#python plot_rates.py s3-md-trm s3-2i2p.s3-md-trm s3s5-2i2p.s3-md-trm s3s5-4adv.s3-md-trm s3-4adv.s3-md-trm


#python plot_rates.py s5-same s5-2i2p.s5-2i2p s3s5-2i2p.s5-2i2p s3s5-4adv.s5-4adv s5-4adv.s5-4adv

#python plot_rates.py s5-md s5-2i2p.s5-md s3s5-2i2p.s5-md s3s5-4adv.s5-md s5-4adv.s5-md

#python plot_rates.py s5-md-trm s5-2i2p.s5-md-trm s3s5-2i2p.s5-md-trm s3s5-4adv.s5-md-trm s5-4adv.s5-md-trm



sim_colors_labels = dict({
    's3-2i2p.s3-2i2p': ['#90EE90','S3-2i2p + 30S-S3-2i2p'],
    's3s5-2i2p.s3-2i2p': ['#87CEFA','S3-2i2p + 30S-S3S5-2i2p'],
    's3s5-4adv.s3-4adv': ['#3CB371','S3-4adv + 30S-S3S5-4adv'],
    's3-4adv.s3-4adv': ['#6495ED','S3-4adv + 30S-S3-4adv'],

    's3-2i2p.s3-md': ['#90EE90','S3-md + 30S-S3-2i2p'],
    's3s5-2i2p.s3-md': ['#87CEFA','S3-md + 30S-S3S5-2i2p'],
    's3s5-4adv.s3-md': ['#3CB371','S3-md + 30S-S3S5-4adv'],
    's3-4adv.s3-md': ['#6495ED','S3-md + 30S-S3-4adv'],

    's3-2i2p.s3-md-trm': ['#90EE90','S3-md-trm + 30S-S3-2i2p'],
    's3s5-2i2p.s3-md-trm': ['#87CEFA','S3-md-trm + 30S-S3S5-2i2p'],
    's3s5-4adv.s3-md-trm': ['#3CB371','S3-md-trm + 30S-S3S5-4adv'],
    's3-4adv.s3-md-trm': ['#6495ED','S3-md-trm + 30S-S3-4adv'],

    's5-2i2p.s5-2i2p': ['#90EE90','S5-2i2p + 30S-S3-2i2p'],
    's3s5-2i2p.s5-2i2p': ['#87CEFA','S5-2i2p + 30S-S3S5-2i2p'],
    's3s5-4adv.s5-4adv': ['#3CB371','S5-4adv + 30S-S3S5-4adv'],
    's5-4adv.s5-4adv': ['#6495ED','S5-4adv + 30S-S3-4adv'],

    's5-2i2p.s5-md': ['#90EE90','S5-md + 30S-S3-2i2p'],
    's3s5-2i2p.s5-md': ['#87CEFA','S5-md + 30S-S3S5-2i2p'],
    's3s5-4adv.s5-md': ['#3CB371','S5-md + 30S-S3S5-4adv'],
    's5-4adv.s5-md': ['#6495ED','S5-md + 30S-S3-4adv'],

    's5-2i2p.s5-md-trm': ['#90EE90','S5-md-trm + 30S-S5-2i2p'],
    's3s5-2i2p.s5-md-trm': ['#87CEFA','S5-md-trm + 30S-S3S5-2i2p'],
    's3s5-4adv.s5-md-trm': ['#3CB371','S5-md-trm + 30S-S3S5-4adv'],
    's5-4adv.s5-md-trm': ['#6495ED','S5-md-trm + 30S-S5-4adv'],
})


keys = ['s3-2i2p.s3-2i2p','s3-2i2p.s3-md','s3-2i2p.s3-md-trm','s3s5-2i2p.s3-2i2p','s3s5-2i2p.s3-md','s3s5-2i2p.s3-md-trm','s3s5-2i2p.s5-2i2p','s3s5-2i2p.s5-md','s3s5-2i2p.s5-md-trm','s5-2i2p.s5-2i2p','s5-2i2p.s5-md','s5-2i2p.s5-md-trm','s5-4adv.s5-2i2p','s5-4adv.s5-md','s5-4adv.s5-md-trm','s3s5-4adv.s5-2i2p','s3s5-4adv.s5-md','s3s5-4adv.s5-md-trm','s3s5-4adv.s3-2i2p','s3s5-4adv.s3-md','s3s5-4adv.s3-md-trm','s3-4adv.s3-2i2p','s3-4adv.s3-md','s3-4adv.s3-md-trm']

rate_data = dict()
fig = matplotlib.pyplot.figure()
ax = fig.add_subplot(1,1,1)
ax.set_xlim([0.,24.])
ax.set_ylim([100000,1000000000])
for sim_name in sim_names:
    with open('../runs/{0}.rates_bootstrap'.format(sim_name),'r') as data_in:
        next(data_in)
        rate_data.update({sim_name: list()})
        for line in data_in:
            entries = line.split()
            if float(entries[0]) < 24:
                rate_data[sim_name].append([float(entries[0]), float(entries[1]), float(entries[2])])


for key in keys:
    if key in rate_data.keys():
        x_angstroms,y_rate,s_dev = zip(*rate_data[key])
        ax.errorbar(x_angstroms, y_rate, yerr=s_dev, c=sim_colors_labels[key][0],label=sim_colors_labels[key][1])


# get handles
handles, labels = ax.get_legend_handles_labels()
# remove the errorbars
handles = [h[0] for h in handles]
ax.legend(handles,labels, loc=2,prop={'size':10})


ax.set_yscale('log')
ax.set_xlabel('Native-Site Separation ($\AA$)')
ax.set_ylabel('Rate Constant ($M^{-1}s^{-1}$)')
ax.set_title('Rates of r-Protein S3 Assocation from BD Simulations')
fig.savefig('../plots/{0}'.format(plot_name),format='pdf')

