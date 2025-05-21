import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colours
import matplotlib as mpl

plt.rcParams['text.usetex'] = True

data = np.arange(-10.0, 10.01, 0.01)

def gen_colour(cmap, n):
    c_map = mpl.colormaps[str(cmap)]
    arr = np.linspace(0, 1, n)
    colourlist = list()
    for c in arr:
        rgba=c_map(c)
        clr = colours.rgb2hex(rgba)
        colourlist.append(str(clr))

    return colourlist

def mag_onevertex(x):
    return (2*np.exp(x) - np.exp(-x)) / (2*np.exp(x) + np.exp(-x))

def mag_twovertex(x):
    return (-6*np.exp(-2*x) + 2*np.exp(2*x)) / (5 + 3*np.exp(-2*x) + np.exp(2*x)) 
#
#def plot_integrand(x):
#    return np.log((np.power((np.exp(2*x) + 1), 2)*np.power((3*np.exp(2*x) + 2), 3))/np.exp(8*x))
#

def plot_integrand_onevertex(x):
    return np.log((2*np.exp(x) + np.exp(-x))/(2*np.exp(-x) + np.exp(x)))

def plot_integrand_twovertex(x):
    return np.log((3*np.exp(-2*x) + np.exp(2*x) + 5)/(3*np.exp(2*x) + np.exp(-2*x) + 5))


mag_data_ny1 = mag_onevertex(data) 

fig, ax = plt.subplots()

ax.set_xlabel(r'$\beta J_R$')
ax.set_ylabel(r'$\langle m \rangle$')
ax.hlines(y=0.0, xmin=-5.0, xmax=5.0, color='black', linestyles='dashed')
ax.vlines(x=0.0, ymin=-2.0, ymax=2.0, color='black', linestyles='dashed')

ax.plot(data, mag_data_ny1, c='green', label=r'$\langle m \rangle$')

fig.legend()
fig.savefig('free_energy_analytic_onevertex.pdf', dpi=300.0)

data = np.arange(0.0, 10.01, 0.01)
mag_int_data = plot_integrand_onevertex(data)
print(np.abs(mag_int_data - np.log(2)))

fig, ax = plt.subplots()
ax.set_xlabel(r'$\beta J$')
ax.set_ylabel(r'$\int \langle m \rangle$')

ax.plot(data, mag_int_data, c='red', label=r'$\int \langle m(x)\rangle dx$')
ax.hlines(y=np.log(3/2), xmin=0.0, xmax=10.0, color='black', linestyles='dashed')

fig.savefig('mag_integrated_onevertex.pdf', dpi=300.0)

data = np.arange(-5.00, 5.01, 0.01)
data2 = np.arange(-0.5, 0.51, 0.01)

fig, ax = plt.subplots()
mag_data_ny2 = mag_twovertex(data) 
mag_data_ny2_inset = mag_twovertex(data2)

left, bottom, width, height = [0.25, 0.55, 0.15, 0.3]
ax2 = fig.add_axes([left, bottom, width, height])
   
ax.set_xlabel(r'$\beta J_R$', fontsize=14)
ax.set_ylabel(r'$\langle m \rangle$', fontsize=14)
ax.hlines(y=0.0, xmin=-np.max(data), xmax=np.max(data), color='black', linestyles='dashed')
ax.vlines(x=0.0, ymin=-2.0, ymax=2.0, color='black', linestyles='dashed')
ax2.hlines(y=0.0, xmin=-np.max(data2), xmax=np.max(data2), linestyles='--', color='black')
ax2.vlines(x=0.0, ymin=np.min(mag_data_ny2_inset), ymax=np.max(mag_data_ny2_inset), linestyles='--', color='black')
ax2.fill([0.0, 0.0, np.log(3)/4], [0.0, -2/5, 0.0], color='green', alpha=0.4)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_xticklabels([])
ax2.set_yticklabels([])

ax.plot(data, mag_data_ny2, c='red')
ax2.plot(data2, mag_data_ny2_inset, c='red')

fig.legend()
fig.savefig('free_energy_analytic_twovertex.pdf', dpi=300.0)
   

data = np.arange(0.0, 10.01, 0.01)
mag_int_data = plot_integrand_twovertex(data)

fig, ax = plt.subplots()
ax.set_xlabel(r'$\delta$', fontsize=14)
ax.set_ylabel(r'$\Delta S / k_B$', fontsize=14)

ax.plot(data, mag_int_data, c='red', label=r'$\int \langle m(x)\rangle dx$')
ax.hlines(y=-np.log(3), xmin=0.0, xmax=10.0, color='black', linestyles='dashed', label=r'$\ln(3)$')
ax2.plot()

fig.savefig('mag_integrated_twovertex.pdf', dpi=300.0)
