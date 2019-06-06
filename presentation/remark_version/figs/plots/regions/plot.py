import numpy as np
import matplotlib.pyplot as plt

# performance
# n = np.arange(16, 130, 16)
#n = [ 4, 8, 12, 16, 20, 32,  48,  64,  80,  96, 112, 128]
n = [8, 16, 32,  48,  64,  80,  96, 112, 128]
size = [24*(i**3) for i in n]
l1 = 32000
l2 = 256000
l3 = 8000000
for i, j in zip(n, size):
   print("n = {:d}, size = {:d}, ".format(i, j))
   if j < l1:
      print("l1 \n")
   elif (j < l2):
      print("l2 \n")
   elif (j < l3):
      print("l3 \n")
   elif (j > l3):
      print("mem \n")
     

# performance for [ 8 16 32 48 64 80 96 112 128]
# with high res image in RAM
perf_noopt=[0.775293, 0.581542, 0.400057, 0.344852, 0.17986, 0.17542, 0.185002, 0.186485, 0.192796]
perf_opt1=[0.765889, 0.603839, 0.488816, 0.356639, 0.19169, 0.181904, 0.180183, 0.190401, 0.189839]

# plot performance
plotname = "regions_performance.png"
fig = plt.figure()
ax = fig.add_subplot(111)

# cache bound
ax.axvline(11, color="tab:gray" )
ax.axvline(22, color="tab:gray" )
ax.axvline(69, color="tab:gray" )
# # performance bound
# ax.axhline(1.0, 0, 69, color="m" )
# ax.axhline(0.5, 69,150, color="m" )
# ax.axvline(69, 0.5,1.0, color="m" )
#
ax.plot(n, perf_noopt, "-ko", label= "baseline")
ax.plot(n, perf_opt1, "-ro", label= "loop unrolling")
#
ax.set_ylim([0, 1.2])
ax.yaxis.grid()
plt.savefig(plotname)
plt.clf()

# Roofline
beta_l1 = 12.0*8.0
beta_l2 = 8.0*8.0
beta_l3 = 4.0*8.0
beta_mem = 2.0*8.0

xi = np.arange(0.0,20.0,0.03125)
yi_l1 = beta_l1*xi 
yi_l2 = beta_l2*xi 
yi_l3 = beta_l3*xi 
yi_mem = beta_mem*xi 

alen = len(xi)
p_wo_simd = np.full(alen, 1.0)

# bounds
plotname = "regions_roofline_bounds.png"
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(xi, p_wo_simd, "b")
# ax.plot(xi, yi_l1, "r--")
# ax.plot(xi, yi_l2, "b--")
ax.plot(xi, yi_l3, "g--")
ax.plot(xi, yi_mem, "k--")

I_n = 1.0/24.0
for i in perf_noopt:
    ax.scatter(I_n, i, marker='o')
for i in perf_opt1:
    ax.scatter(I_n, i, marker='*')

ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)

ax.set_xlim([2**(-5),2**3])
ax.set_ylim([2**(-4),2**5])

# ax.set_xlabel("Operational Intensity [flops/bytes]")
# ax.set_ylabel("Performance [flops/cycle]")

plt.savefig(plotname)
plt.clf()
