import numpy as np
import matplotlib.pyplot as plt

with open("mil_benchmark.csv","r") as f:
   data = f.readlines()

data.pop(0)

# n = []
# simd_nonblocking = []
# simd_blocking = []
# blocking = []
# nonblocking = []
# baseline = []

# for line in data:
#     values = line.strip().split()
#     n.append(int(values[0]) )
#     simd_nonblocking.append(float(values[1]))
#     simd_blocking.append(float(values[2]))
#     blocking.append(float(values[3]))
#     nonblocking.append(float(values[4]))
#     baseline.append( float(values[5]))

n = [16, 32, 48, 64, 80, 96, 112, 128, 144, 160, 176, 192, 208, 224, 240, 256, 272, 288, 304, 320, 336, 352]
baseline = [0.253126, 0.273263,0.285723,0.268263,0.292533,0.22245,0.190012,0.141454,0.154034,0.14583,0.130745,0.131468,0.132991,0.13068,0.129215,0.11497,0.125344,0.12235,0.117174,0.115357,0.122026,0.120587]
nonblocking = [ 0.507966,0.518339,0.442633,0.357604,0.272747,0.267925,0.204402,0.166334,0.178962,0.173255,0.178665,0.152704,0.168882,0.149292,0.140991,0.115668,0.127536,0.120879,0.123546,0.0945074,0.118324,0.0977236]
blocking = [ 0.484023,0.471778,0.481992,0.322162,0.445875,0.42407,0.417652,0.28337,0.426249,0.404632,0.404399,0.293365,0.397039,0.386905,0.387979,0.277467,0.397114,0.37891,0.379555,0.298291,0.378078,0.379609]
simd = [ 1.30327,1.13983,1.24452,0.616068,1.07943,0.95392,0.581806,0.474591,0.945199,0.874368,0.857527,0.612672,0.836425,0.834621,0.806127,0.511187,0.809358,0.780495,0.791103,0.598147,0.778034,0.748431]



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
ax.plot(n, simd, "-ro", label= "simd blocking")
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
# ax.scatter(6.5, 0.25, marker='o') # baseline
# ax.scatter(6.5, 0.65, marker='o') # blocking with accumulators
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
