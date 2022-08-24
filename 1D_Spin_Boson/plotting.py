import matplotlib.pyplot as plt
import numpy as np
fig, ax = plt.subplots() # 创建图实例
from matplotlib.pyplot import MultipleLocator, tick_params

lw = 2.0
legendsize = 32         # size for legend
font_legend = {'family':'Times New Roman', 'weight': 'roman', 'size':18}
color1 = 'black'            # 
color2 = 'r'

# ==============================================================================================
#                                      Fig 1a    
# ==============================================================================================

Unitlen = 5
fig = plt.figure(figsize=(4.5 * Unitlen, Unitlen), dpi = 512)

plt.subplot(1, 3, 1)

data = np.loadtxt("Pt_Model1_3.txt", dtype = float)                                  
plt.plot(data[:,0], data[:,1], "--", linewidth = lw,  color = 'green', label = "tiers = 2")

data = np.loadtxt("Pt_Model1_5.txt", dtype = float)                                  
plt.plot(data[:,0], data[:,1], "--", linewidth = lw,  color = 'blue', label = "tiers = 4")

data = np.loadtxt("Pt_Model1_10.txt", dtype = float)                                  
plt.plot(data[:,0], data[:,1], "-", linewidth = lw,  color = color1, label = "tiers = 10")

data0 = np.loadtxt("DVR_1.txt", dtype = float)
plt.plot(data0[:,0], data0[:,1], ".", color = color2, label = "DVR")

# ==============================================================================================

time = 30.0             # x-axis range: (0, time)
y1, y2 = -1.5, 2.0     # y-axis range: (y1, y2)

# scale for major and minor locator
x_major_locator = MultipleLocator(10)
x_minor_locator = MultipleLocator(2)
y_major_locator = MultipleLocator(0.5)
y_minor_locator = MultipleLocator(0.1)

# x-axis and LHS y-axis
ax = plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.yaxis.set_minor_locator(y_minor_locator)
ax.tick_params(which = 'major', length = 8, labelsize = 10)
ax.tick_params(which = 'minor', length = 4)

x1_label = ax.get_xticklabels()
[x1_label_temp.set_fontname('Times New Roman') for x1_label_temp in x1_label]
y1_label = ax.get_yticklabels()
[y1_label_temp.set_fontname('Times New Roman') for y1_label_temp in y1_label]

plt.tick_params(labelsize = 20, which = 'both', direction = 'in')
plt.xlim(0.0, time)
plt.ylim(y1, y2)

# RHS y-axis
ax2 = ax.twinx()
ax2.yaxis.set_major_locator(y_major_locator)
ax2.yaxis.set_minor_locator(y_minor_locator)
ax2.tick_params(which = 'major', length = 8)
ax2.tick_params(which = 'minor', length = 4)
ax2.axes.yaxis.set_ticklabels([])

plt.tick_params(which = 'both', direction = 'in')
plt.ylim(y1, y2)

# name of x, y axis and the panel
ax.set_xlabel('tΔ', font = 'Times New Roman', size = 18)
ax.set_ylabel('P(t)', font = 'Times New Roman', size = 18)

# legend location, font & markersize
plt.legend(title = '(a)', frameon = False, title_fontsize = legendsize)
ax.legend(loc = 'upper center', frameon = False, prop = font_legend) # lower right

# ==============================================================================================
#                                      Fig 1b     
# ==============================================================================================

plt.subplot(1, 3, 2)

# data = np.loadtxt("Pt_Model2_3.txt", dtype = float)                                  
# plt.plot(data[:,0], data[:,1], "--", linewidth = lw,  color = 'green', label = "tiers = 2")

data = np.loadtxt("Pt_Model2_5.txt", dtype = float)                                  
plt.plot(data[:,0], data[:,1], "--", linewidth = lw,  color = 'green', label = "tiers = 4")

data = np.loadtxt("Pt_Model2_11.txt", dtype = float)                                  
plt.plot(data[:,0], data[:,1], "--", linewidth = lw,  color = 'blue', label = "tiers = 10")

data = np.loadtxt("Pt_Model2_20.txt", dtype = float)                                  
plt.plot(data[:,0], data[:,1], "-", linewidth = lw,  color = color1, label = "tiers = 20")

data0 = np.loadtxt("DVR_2.txt", dtype = float)
plt.plot(data0[:,0], data0[:,1], ".", color = color2, label = "DVR")

# ==============================================================================================

time = 20.0             # x-axis range: (0, time)
y1, y2 = -1.5, 2.5     # y-axis range: (y1, y2)

# scale for major and minor locator
x_major_locator = MultipleLocator(5)
x_minor_locator = MultipleLocator(1)
y_major_locator = MultipleLocator(0.5)
y_minor_locator = MultipleLocator(0.1)

# x-axis and LHS y-axis
ax = plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.yaxis.set_minor_locator(y_minor_locator)
ax.tick_params(which = 'major', length = 8, labelsize = 10)
ax.tick_params(which = 'minor', length = 4)

x1_label = ax.get_xticklabels()
[x1_label_temp.set_fontname('Times New Roman') for x1_label_temp in x1_label]
y1_label = ax.get_yticklabels()
[y1_label_temp.set_fontname('Times New Roman') for y1_label_temp in y1_label]

plt.tick_params(labelsize = 20, which = 'both', direction = 'in')
plt.xlim(0.0, time)
plt.ylim(y1, y2)

# RHS y-axis
ax2 = ax.twinx()
ax2.yaxis.set_major_locator(y_major_locator)
ax2.yaxis.set_minor_locator(y_minor_locator)
ax2.tick_params(which = 'major', length = 8)
ax2.tick_params(which = 'minor', length = 4)
ax2.axes.yaxis.set_ticklabels([])

plt.tick_params(which = 'both', direction = 'in')
plt.ylim(y1, y2)

# name of x, y axis and the panel
ax.set_xlabel('tΔ', font = 'Times New Roman', size = 18)

plt.legend(title = '(b)', frameon = False, title_fontsize = legendsize)
ax.legend(loc = 'upper center', frameon = False, prop = font_legend) # lower right

# ==============================================================================================
#                                      Fig 1c 
# ==============================================================================================

plt.subplot(1, 3, 3)

# data = np.loadtxt("Pt_Model3_3.txt", dtype = float)                                  
# plt.plot(data[:,0], data[:,1], "--", linewidth = lw,  color = 'green', label = "tiers = 2")
# 
data = np.loadtxt("Pt_Model3_51.txt", dtype = float)                                  
plt.plot(data[:,0], data[:,1], "--", linewidth = lw,  color = 'blue', label = "tiers = 50")

data = np.loadtxt("Pt_Model3_75.txt", dtype = float)                                  
plt.plot(data[:,0], data[:,1], "-", linewidth = lw,  color = color1, label = "tiers = 75")

data0 = np.loadtxt("DVR_3.txt", dtype = float)
plt.plot(data0[:,0], data0[:,1], ".", color = color2, label = "DVR")

# ==============================================================================================

time = 20.0             # x-axis range: (0, time)
y1, y2 = -0.6, 1.0     # y-axis range: (y1, y2)

# scale for major and minor locator
x_major_locator = MultipleLocator(5)
x_minor_locator = MultipleLocator(1)
y_major_locator = MultipleLocator(0.5)
y_minor_locator = MultipleLocator(0.1)

# x-axis and LHS y-axis
ax = plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.yaxis.set_minor_locator(y_minor_locator)
ax.tick_params(which = 'major', length = 8, labelsize = 10)
ax.tick_params(which = 'minor', length = 4)

x1_label = ax.get_xticklabels()
[x1_label_temp.set_fontname('Times New Roman') for x1_label_temp in x1_label]
y1_label = ax.get_yticklabels()
[y1_label_temp.set_fontname('Times New Roman') for y1_label_temp in y1_label]

plt.tick_params(labelsize = 20, which = 'both', direction = 'in')
plt.xlim(0.0, time)
plt.ylim(y1, y2)

# RHS y-axis
ax2 = ax.twinx()
ax2.yaxis.set_major_locator(y_major_locator)
ax2.yaxis.set_minor_locator(y_minor_locator)
ax2.tick_params(which = 'major', length = 8)
ax2.tick_params(which = 'minor', length = 4)
ax2.axes.yaxis.set_ticklabels([])

plt.tick_params(which = 'both', direction = 'in')
plt.ylim(y1, y2)

# name of x, y axis and the panel
ax.set_xlabel('tΔ', font = 'Times New Roman', size = 18)

plt.legend(title = '(c)', frameon = False, title_fontsize = legendsize)
ax.legend(loc = 'upper center', frameon = False, prop = font_legend) # lower right

plt.savefig("figure_1.png", bbox_inches='tight')
# plt.show()