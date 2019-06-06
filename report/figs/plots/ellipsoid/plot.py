import numpy as np
import matplotlib.pyplot as plt


n = np.array([      4,      16,      64,     256,    1024,    4096,   16384,
         65536,  262144, 1048576])


pref = np.array([[ 1.59817,  2.79499,  2.51165,  5.12104],
       [ 1.65723,  3.0477 ,  2.70167,  7.67245],
       [ 1.67332,  3.15797,  2.78886,  9.45952],
       [ 1.66822,  3.18045,  2.81027, 10.0711 ],
       [ 1.67489,  3.19695,  2.8088 , 10.1898 ],
       [ 1.61189,  3.08604,  2.79606,  9.0032 ],
       [ 1.55117,  2.99381,  2.78885,  8.51551],
       [ 1.53784,  2.98256,  2.7844 ,  8.44633],
       [ 1.50715,  2.93291,  2.7506 ,  8.38612],
       [ 1.26723,  2.45531,  2.42534,  6.13949]])
pref = np.transpose(pref)

funcs = ['fit_ellipsoid', 'fit_ellipsoid_opt', 'fit_ellipsoid_simd_points_shuffles', 'fit_ellipsoid_simd_points']

# plot performance
plotname = "ellipsoid_performance.png"
fig = plt.figure()
ax = fig.add_subplot(111)

# cache bound
ax.axvline(1365, color="tab:gray" )
ax.axvline(10922, color="tab:gray" )
ax.axvline(349525, color="tab:gray" )
#
ax.plot(n, pref[0][:], "-k*", label= "baseline")
ax.plot(n, pref[1][:], "-g*", label= "optimized")
ax.plot(n, pref[2][:], "-b*", label= "simd points shuffles")
ax.plot(n, pref[3][:], "-r*", label= "simd points")
#
ax.set_xscale('log', basex=2)
ax.set_ylim([0,12])
ax.yaxis.grid()
#

#
plt.savefig(plotname)
plt.clf()
