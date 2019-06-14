import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as ticker

# performance
perfs = np.array([
[16, 0.892558, 3.05728],
[48, 0.292956, 0.910627],
[80, 0.265725, 0.689911],
[112, 0.19738, 0.595828],
[144, 0.150848, 0.583671],
])

perfs_t = perfs.T
n = perfs_t[0]
perf_noopt = perfs_t[1]
perf_opt1 = perfs_t[2]
speedup = [ b/a for a, b in zip(perf_noopt, perf_opt1)]
for i, j in enumerate(n):
   print(i, speedup[i],j )

# plot performance
plotname = "overall_performance.png"
fig = plt.figure()
ax = fig.add_subplot(111)

# cache bound
ax.axvline(31, color="tab:gray" )
ax.axvline(115, color="tab:gray" )
#
ax.plot(n, perf_noopt, "-ko", label= "baseline")
ax.plot(n, perf_opt1, "-ro", label= "optimized")
#
ax.set_ylim([0, 3.5])
#
tick_spacing = 0.5
ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
#
ax.yaxis.grid()
#
plt.savefig(plotname)
plt.clf()
