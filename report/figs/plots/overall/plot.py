import numpy as np
import matplotlib.pyplot as plt

# performance
n = [16, 32, 80]
# performance for 
# with high res image in RAM
perf_noopt=[1.71839, 0.862286,  0.660607 ]
perf_opt1=[7.35537, 2.6, 1.83053]

# plot performance
plotname = "ovreall_performance.png"
fig = plt.figure()
ax = fig.add_subplot(111)

#
ax.plot(n, perf_noopt, "-ko", label= "baseline")
ax.plot(n, perf_opt1, "-ro", label= "optimized")
#
ax.set_ylim([0, 10])
ax.yaxis.grid()
#
plt.savefig(plotname)
plt.clf()

# # Roofline
# beta_l1 = 12.0*8.0
# beta_l2 = 8.0*8.0
# beta_l3 = 4.0*8.0
# beta_mem = 2.0*8.0

# xi = np.arange(0.0,20.0,0.03125)
# yi_l1 = beta_l1*xi 
# yi_l2 = beta_l2*xi 
# yi_l3 = beta_l3*xi 
# yi_mem = beta_mem*xi 

# alen = len(xi)
# p_wo_simd = np.full(alen, 1.0)

# # bounds
# plotname = "regions_roofline_bounds.png"
# fig = plt.figure()
# ax = fig.add_subplot(111)

# ax.plot(xi, p_wo_simd, "b")
# # ax.plot(xi, yi_l1, "r--")
# # ax.plot(xi, yi_l2, "b--")
# ax.plot(xi, yi_l3, "g--")
# ax.plot(xi, yi_mem, "k--")

# I_n = 1.0/24.0
# for i in perf_noopt:
#     ax.scatter(I_n, i, marker='o')
# for i in perf_opt1:
#     ax.scatter(I_n, i, marker='*')

# ax.set_xscale('log', basex=2)
# ax.set_yscale('log', basey=2)

# ax.set_xlim([2**(-5),2**3])
# ax.set_ylim([2**(-4),2**5])

# # ax.set_xlabel("Operational Intensity [flops/bytes]")
# # ax.set_ylabel("Performance [flops/cycle]")

# plt.savefig(plotname)
# plt.clf()
