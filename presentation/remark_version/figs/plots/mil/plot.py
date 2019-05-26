import numpy as np
import matplotlib.pyplot as plt

with open("mil_benchmark.csv","r") as f:
   data = f.readlines()

data.pop(0)

n = []
simd_nonblocking = []
simd_blocking = []
blocking = []
nonblocking = []
baseline = []

for line in data:
    values = line.strip().split()
    n.append(int(values[0]) )
    simd_nonblocking.append(float(values[1]))
    simd_blocking.append(float(values[2]))
    blocking.append(float(values[3]))
    nonblocking.append(float(values[4]))
    baseline.append( float(values[5]))


p = [0.0, 2.0]
n_l2 = [31, 31]
n_l3 = [115, 115]

# plot performance
plotname = "mil_performance.png"
fig = plt.figure()
ax = fig.add_subplot(111)
#
ax.axvline(31, color="tab:gray" )
ax.axvline(115, color="tab:gray" )
#
ax.plot(n, baseline, "-ko", label= "baseline")
# ax.plot(n, simd_nonblocking, "-ro", label= "simd nonblocking")
ax.plot(n, simd_blocking, "-ro", label= "simd blocking")
ax.plot(n, blocking, "-go", label= "blocking")
ax.plot(n, nonblocking, "-bo", label= "nonblocking")
#
ax.yaxis.grid()
ax.set_ylim([0, 2.0])
#
plt.savefig(plotname)
plt.clf()

# Roofline
beta_l1 = 12.0*8.0
beta_l2 = 8.0*8.0
beta_l3 = 4.0*8.0
beta_mem = 2.0

xi = np.arange(0.0,20.0,0.03125)
yi_l1 = beta_l1*xi 
yi_l2 = beta_l2*xi 
yi_l3 = beta_l3*xi 
yi_mem = beta_mem*xi 

alen = len(xi)
p_wo_simd = np.full(alen, 2.0)

# bounds
plotname = "mil_roofline.png"
fig = plt.figure()
ax = fig.add_subplot(111)

#ax.plot(xi, p_wo_simd, "b")
ax.plot(xi, yi_mem, "k")

ax.axhline(8, color="b" )
ax.axhline(2, color="b" )

ax.scatter(0.5, 0.15, marker='o') # baseline/simd nonblocking
ax.scatter(6.5, 0.25, marker='o') # baseline
ax.scatter(6.5, 0.65, marker='o') # blocking with accumulators
ax.scatter(6.5, 1.0, marker='o') # simd blocking large n
ax.scatter(6.5, 1.7, marker='o') # simd blocking max performance

# I_n = 1.0/24.0
# for i in perf_noopt:
#     ax.scatter(I_n, i, marker='o')
# for i in perf_opt1:
#     ax.scatter(I_n, i, marker='*')

ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)

ax.set_xlim([2**(-2),2**5])
ax.set_ylim([2**(-4),2**4])

# ax.set_xlabel("Operational Intensity [flops/bytes]")
# ax.set_ylabel("Performance [flops/cycle]")

plt.savefig(plotname)
plt.clf()
