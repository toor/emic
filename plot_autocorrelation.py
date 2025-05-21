import numpy as np 
import scipy as sp 
import pandas as pd 
import matplotlib.pyplot as plt
from statsmodels.graphics.tsaplots import plot_acf 
from statsmodels.tsa.stattools import acf
from scipy.optimize import curve_fit
import sys
import os

plt.rcParams['text.usetex'] = True

def rounded_exponential(x):
    return int(np.sign(x)*np.exp(np.abs(x)))

def damped_exp(t, tau):
    return np.exp(-t/tau)

def next_pow_two(n):
    i = 1 
    while i < n:
        i = i << 1 
    return i

# Automated windowing procedure following Sokal (1989)
def auto_window(taus, c):
    m = np.arange(len(taus)) < c * taus
    if np.any(m):
        return np.argmin(m)
    return len(taus) - 1


# Following the suggestion from Goodman & Weare (2010)
def autocorr_gw2010(y, c=5.0):
    f = autocorr_func_1d(np.mean(y, axis=0), False)
    taus = 2.0 * np.cumsum(f) - 1.0
    window = auto_window(taus, c)
    return taus[window]


def autocorr_new(y, c=5.0):
    f = np.zeros_like(y)
    for yy in y:
        f += autocorr_func_1d(yy, False)
    f /= len(y)
    taus = 2.0 * np.cumsum(f) - 1.0
    window = auto_window(taus, c)
    return taus[window]

def autocorr_func_1d(x, norm=True):
    x = np.atleast_1d(x)
    if len(x.shape) != 1:
        raise ValueError("invalid dimensions for 1d correlation function")
    n = next_pow_two(len(x))
    
    f = np.fft.fft(x - np.mean(x), n = 2*n)
    acf = np.fft.ifft(f * np.conjugate(f))[:len(x)].real
    acf /= 4 * n 
    if norm:
        acf /= acf[0] 
    return acf

def calc_autocorrelation_times(n_y, n_x, betaJ):
    taus = np.zeros_like(n_x, dtype=np.float64)
    for i, n in enumerate(n_x):
        filename = f'results/SIM_241018_0002/N_y={n_y}/N_x={n}/betaJ={betaJ}/rep1/data.csv'
        df = pd.read_csv(filename, sep='\s+')
        t = df['iter'].to_numpy()
        rho = df['rho'].to_numpy() 
        
        ac = autocorr_func_1d(rho)
        
        popt, pcov = curve_fit(damped_exp, t, ac)
        taus[i] = popt[0]

    return taus

def calc_autocorrelations(n_x, n_y, betaJ):
    filename = f'results/SIM_241018_0001/N_y={n_y}/N_x={n_x}/betaJ={betaJ:.1f}/rep1/data.csv'
    df = pd.read_csv(filename, sep='\s+')

    t = df['iter'].to_numpy()
    rho = df['rho'].to_numpy()

    ac = autocorr_func_1d(rho)
    print(ac)
    return t, ac

def generate_x_positions(y, max_x, step):
    return np.arange(y, max_x + step, step)


n_y = 2 
n_x_vals = np.array([2])
betaJ_vals = np.array([-10.0, -1.0, 0.0, 1.0, 10.0])

taus = np.zeros((n_x_vals.size, len(betaJ_vals)))
file_path = f'figs/SIM_2410018_0001/autocorrelations/N_y={n_y}/'
os.makedirs(file_path, exist_ok=True)

for j, betaJ in enumerate(betaJ_vals):
    print(f'betaJ = {betaJ}')
    for i, n_x in enumerate(n_x_vals):
        fig, axs  = plt.subplots()
        print(f'N_x = {n_x}')
        t, ac = calc_autocorrelations(n_x, n_y, betaJ)
        
       # popt, pcov = curve_fit(damped_exp, t, ac)
       # tau = popt[0]
       # std_tau = np.sqrt(np.diag(pcov)[0])
       # taus[i,j] = tau 
    
        axs.scatter(t, ac)
        axs.set_xlabel("Time")
        axs.set_ylabel("Autocorrelation")
        fig.legend()
    
        path = f'figs/SIM_241018_0002/autocorrelations/N_y={n_y}'
        os.makedirs(path, exist_ok=True)
    
        fig.savefig(f'figs/SIM_2410018_0001/autocorrelations/N_y={n_y}/autocorrelation_N_x={n_x}_betaJ={betaJ}.pdf', dpi=300.0)
        plt.close()
    # fig, axs = plt.subplots()
    # axs.scatter(n_x_vals, taus[:,j])
    # axs.set_title(rf'$N_y = {n_y}$; $\beta J = {betaJ}$')
    # fig.savefig(f'figs/SIM_241003_0002/autocorrelations/N_y={n_y}/autocorrelation_betaJ={betaJ}.pdf', dpi=300.0) 
