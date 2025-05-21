import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt 
import pandas as pd 
from numba import jit
import matplotlib as mpl
import matplotlib.colors as colours
import sys
import os
from scipy.integrate import simpson

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb}')
#mpl.verbose.level = 'debug-annoying'

def gen_colour(cmap, n):
    c_map = mpl.colormaps[str(cmap)]
    arr = np.linspace(0, 1, n)
    colourlist = list()
    for c in arr:
        rgba=c_map(c)
        clr = colours.rgb2hex(rgba)
        colourlist.append(str(clr))

    return colourlist
class Data:
    def __init__(self, prefix, n_y_vals, max_x, betaJ_vals, step, reps):
        self.prefix = prefix
        self.n_y_vals = n_y_vals
        self.max_x = max_x
        self.betaJ_vals = betaJ_vals
        self.Jmax = np.max(betaJ_vals)
        self.step = step
        self.reps = reps

    def generate_x_positions(self, n_y):
        return np.arange(n_y, self.max_x + self.step, self.step)
    
    # This is a shifted sigmoid function, which is indeed the shape we obtain for the magnetisation
    # as a function of the coupling. Avoids complications involved with computing the integral
    # numerically
    def magnetisation_fit(self, x, A, B, C, D):
        return A*np.tanh(B*(x - C)) + D
    
    # This is the definite integral of the sigmoid function defined in `magnetisation_fit`,
    # integrated over the symmetric interval [-Jmax, Jmax]
    def integrated_magnetisation(self, A, B, C, D):
        x = self.Jmax
        return A * np.log(np.cosh(B*(x - C))/np.cosh(B*(x + C))) / B + 2*D*x
    
    # Propagate the error, assuming no correlation between fitted parameters:
    # see https://doi.org/10.6028/jres.070c.025
    def error_in_integral(self, A, B, C, D, sigma_A, sigma_B, sigma_C, sigma_D):
        delta = self.Jmax
        # These are the partial derivatives of the integrated magnetisation 
        # with respect to each parameter A, B, C, D
        d_dA = np.log(np.cosh(B*(delta - C))/np.cosh(B*(delta + C))) / B 
        d_dB = -(A/B**2)*np.log(np.cosh(B*(delta - C))/np.cosh(B*(delta + C))) + (A / B) * (np.tanh(B*(delta - C)) - np.tanh(B*(delta + C)))
        d_dC = -A*(np.tanh(B*(delta - C)) + np.tanh(B*(delta + C)))
        d_dD = 2*delta 

        
        return np.sqrt((d_dA * sigma_A)**2 + (d_dB * sigma_B)**2 + (d_dC * sigma_C)**2 + (d_dD * sigma_D)**2)
    
    # Compute the entropy difference between parallel
    # and antiparallel configurations of the macrospins, 
    # using \delta S / k_B = \beta\delta_E - \delta F,
    # where delta_F corresponds to the magnetisation of 
    # the right interface, integrated over a sufficiently large range of values of the coupling.
    def calculate_entropy(self, n_y):
        n_x_vals = self.generate_x_positions(n_y)
        delta_S = np.zeros((n_x_vals.size, 3), dtype=np.float64)

        for i, n_x in enumerate(n_x_vals):
            #print(f"N_x = {n_x}")
            total_spins = n_x * n_y 
            s_vals = np.zeros((betaJ_vals.size,2), dtype=np.float64)
            for j, betaJ in enumerate(betaJ_vals):
                avg_mags = np.zeros((self.reps.size, 2))
                for k, rep in enumerate(reps):
                    filename = 'results/' + self.prefix + f'/N_y={n_y}/N_x={n_x}/betaJ={betaJ:.1f}/rep{rep}/data.csv'
                    df = pd.read_csv(filename, sep='\s+')
                    avg_mags[k][0] = df['m_int'].mean()
                    avg_mags[k][1] = df['m_int'].std()
                    #print(f'std for betaJ = {betaJ}, rep {rep} is {avg_mags[k][1]}')
                s_vals[j][0] = np.mean(avg_mags[:,0])
                s_vals[j][1] = np.sqrt(np.sum(np.power(avg_mags[:,1], 2))) / (avg_mags.shape[0])
            
            # Replace any standard deviations which are extremely close to 0.0 with 0.1 
            s_vals[:,1] = np.where(np.isclose(s_vals[:,1], 0.0, atol=1e-9), 0.1, s_vals[:,1])
            # First fit the correct parameters to the magnetisation curve itself. Extract fitted parameters and their standard errors.
            popt, pcov = sp.optimize.curve_fit(self.magnetisation_fit,
                                               xdata=self.betaJ_vals,
                                               ydata=s_vals[:,0],
                                               sigma=s_vals[:,1])
            [A, B, C, D] = popt
            [sigma_A, sigma_B, sigma_C, sigma_D] = np.sqrt(np.diag(pcov))
            
            # Now, use these extracted parameters to integrate the magnetisation 
            # over the range [-Jmax, Jmax]
            delta_S[i][0] = total_spins * self.integrated_magnetisation(A, B, C, D)
            delta_S[i][1] = simpson(self.betaJ_vals, s_vals[:,0])
            delta_S[i][2] = total_spins * self.error_in_integral(A, B, C, D, sigma_A, sigma_B, sigma_C, sigma_D)

        return delta_S

    def calculate_free_energy(self, n_y, n_x_vals, J_vals, beta, reps):
        delta_F = np.zeros_like((n_x_vals.size, 2), dtype=np.float64)
    
        for i, n_x in enumerate(n_x_vals):
            total_spins = n_x * n_y
            s_vals = np.zeros((J_vals.size, 2), dtype=np.float64)
            for j, J in enumerate(J_vals):
                avg_mags = np.zeros((reps.size, 2))
                # Average over repeat measurements
                for k, rep in enumerate(reps):
                    filename = 'results/' + self.prefix + f'/N_y={n_y}/N_x={n_x}/betaJ={betaJ:.1f}/rep{rep}/data.csv'
                    df = pd.read_csv(filename, sep='\s+')
                    avg_mags[k][0] = df['m_int'].mean()
                    avg_mags[k][1] = df['m_int'].var()
                s_vals[j][0] = np.sum(avg_mags[:,0])
                s_vals[j][1] = np.sum(avg_mags[:,1]) # sum variances
            delta_F[i][0] = total_spins * sp.integrate.simpson(s_vals, J_vals)
            delta_F[i][1] = total_spins * error_in_integral(J_vals, s_vals[:,1])
    
        return delta_F

    def calculate_extreme_energies(self, n_y):
        n_x_vals = self.generate_x_positions(n_y)
        deltaE_vals = np.zeros((n_x_vals.size, 2), dtype=np.float64)
        for i, n_x in enumerate(n_x_vals):
            total_spins = n_x * n_y
            E_min = self.process_energy(self.prefix, J_min, num_reps, n_x, n_y)
            E_max = self.process_energy(self.prefix, J_max, num_reps, n_x, n_y)
            
            deltaE_vals[i][0] = (E_max[0] - E_min[0])
            deltaE_vals[i][1] = (E_max[1] + E_max[1])
    
        return deltaE_vals

    def calc_energies(self, J_vals, n_x, n_y, reps):
        energies = np.zeros_like(J_vals)
        for i, J in enumerate(J_vals):
            rep_energies = np.zeros((reps,))
            for r in range(reps):
                filename = 'results/' + self.prefix + f'/N_y={n_y}/N_x={n_x}/betaJ={J:.2f}/rep{r+1}/data.csv'
                rep_energies[r] = pd.read_csv(filename, sep='\s+')['energy'].mean()
            energies[i] = np.mean(rep_energies)
        return energies
    
    def calc_magnetisations(self, n_x, n_y):
        N = n_x * n_y
        mags = np.zeros((self.betaJ_vals.size, 2), dtype=np.float64)
        for i, J in enumerate(self.betaJ_vals):
            rep_mags = np.zeros_like(self.reps, dtype=np.float64)
            for r in self.reps:
                filename = 'results/' + self.prefix + f'/N_y={n_y}/N_x={n_x}/betaJ={J:.1f}/rep{r}/data.csv'
                rep_mags[r-1] = pd.read_csv(filename, sep='\s+')['m_int'].mean() 
            mags[i][0] = np.mean(rep_mags)
            mags[i][1] = np.std(rep_mags)
        return N * mags

    def process_energy(self, betaJ, n_x, n_y):
        repeat_energies = np.zeros((self.reps.size,2), dtype=np.float64)
    
        for r in self.reps:
            # account for zero-based indexing.
            filename = 'results/' + self.prefix + f'/N_y={n_y}/N_x={n_x}/betaJ={betaJ:.1f}/rep{rep}/data.csv'
            df = pd.read_csv(filename, sep='\s+')
            repeat_energies[r-1][0] = df['energy'].mean()
            repeat_energies[r-1][1] = df['energy'].var()
        
        # Sum the variances of the individual (repeat) measurements
        return (np.mean(repeat_energies[:,0]), np.sum(repeat_energies[:,1]))
    
    def linear_fit(self, x, a, b):
        return a*x + b
    
    # This function does all the heavy lifting,
    # calculating the entropy change for all combinations of values of N_y and N_x
    # and generating a plot of entropy vs N_x for each value of N_y.
    def calc_thermodynamic_quantities(self):
        file_path = 'figs/' + self.prefix
        os.makedirs(file_path, exist_ok=True)
        
        fig, ax = plt.subplots()
        ax.set_xlabel(r'$N_x - N_y$', fontsize=14)
        ax.set_ylabel(r'$\Delta S / k_B$', fontsize=14)
        ax.hlines(y=0.0, xmin=0, xmax=self.max_x - np.min(self.n_y_vals), linestyles='--', colors='black')        
        x_ticks = self.generate_x_positions(np.min(self.n_y_vals)) - np.min(self.n_y_vals)
        ax.set_xticks(ticks=x_ticks, labels=x_ticks)
        
        fig1, ax1 = plt.subplots()
        
        max_delta_S = np.zeros((self.n_y_vals.size, 2), dtype=np.float64)
        cmap = colours.ListedColormap(plt.cm.viridis(np.linspace(0,1,6)))
        norm = colours.BoundaryNorm(boundaries=np.arange(1.5, 8.5), ncolors=6)

        for i, n_y in enumerate(self.n_y_vals):
            #print(f"Generating entropy vs N_x for N_y = {n_y}")
            n_x_vals = self.generate_x_positions(n_y)
            delta_S = self.calculate_entropy(n_y)
            max_delta_S[i][0] = np.min(delta_S[:,0])
            max_delta_S[i][1] = delta_S[np.argmin(delta_S[:,0])][2] 

            ax.scatter(n_x_vals - n_y, delta_S[:,0], color=cmap(i), marker='x', label=rf'Analytic $N_y = {n_y}$')
            ax1.scatter(n_x_vals, np.exp(delta_S[:,0]), color=cmap(i), marker='x', label=rf'$N_y = {n_y}$')
        
        # Create inset plot, and add it to the figure.
        rect = (0.35, 0.25, 0.4, 0.3)
        ax2 = fig.add_axes(rect)
        
        # Attempt a linear fit of the maximum value of entropy change per value of N_y to N_y. 
        popt, pcov = sp.optimize.curve_fit(self.linear_fit, xdata=self.n_y_vals, ydata=max_delta_S[:,0])#, sigma=max_delta_S[:,1])
        [a, b] = popt 
        
        # Add a colour bar
        cax = fig.add_axes([0.5, 0.65, 0.2, 0.1])
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = plt.colorbar(sm, cax=cax, orientation='horizontal', ticks=self.n_y_vals)
        cbar.set_label(r'$N_y$')
        cbar.set_ticks(self.n_y_vals)
        cbar.set_ticklabels(self.n_y_vals)


        ax2.set_xlabel(r'$N_y$')
        ax2.set_ylabel(r'Max. $\Delta S/k_B$')
        ax2.scatter(self.n_y_vals, max_delta_S[:,0], marker='x', c='red')
        ax2.errorbar(self.n_y_vals, max_delta_S[:,0], yerr=max_delta_S[:,1], c='black', capsize=2.0, ls='none')
        ax2.plot(self.n_y_vals, self.linear_fit(self.n_y_vals, a, b), c='blue', ls='--')

        fig.savefig('figs/' + self.prefix + '/entropy_vs_Nx.pdf', dpi=300.0)
        
        fig1.savefig('figs/' + self.prefix + '/exp_linear_fit_entropy.pdf', dpi=300.0)

        plt.close()

    def make_magnetisation_plots(self, n_y):
        file_path = 'figs/' + self.prefix + f'/magnetisation_fits/N_y={n_y}'
        os.makedirs(file_path, exist_ok=True)
        n_x_vals = self.generate_x_positions(n_y)
        #print(n_x_vals)

        for n_x in n_x_vals:
            fig, ax = plt.subplots()
            ax.set_xlabel(r'$\beta J$')
            ax.set_ylabel(r'$\langle m \rangle$')
            ax.hlines(y=0.0, xmin=-self.Jmax, xmax=self.Jmax, linestyles='--', colors='black')
            ax.vlines(x=0.0, ymin=-n_y, ymax=n_y, linestyles='--', colors='black')
            mags = self.calc_magnetisations(n_x, n_y)
            popt, pcov = sp.optimize.curve_fit(self.magnetisation_fit, self.betaJ_vals, mags[:,0])
            [A, B, C, D] = popt
            ax.scatter(self.betaJ_vals, mags[:,0], c='red', marker='x', label='Numerical data')
            ax.plot(self.betaJ_vals, mags[:,0], c='blue', linestyle='--', label='Sigmoid fit')
            fig.legend()
            fig.savefig('figs/' + self.prefix + f'/magnetisation_fits/N_y={n_y}/magnetisation_N_x={n_x}.pdf', dpi=300.0)
            plt.close()

    def make_magnetisation_histograms(self, n_y):
        file_path = 'figs/' + self.prefix + f'/magnetisation_histograms/N_y={n_y}'
        os.makedirs(file_path, exist_ok=True)

        #n_x_vals = self.generate_x_positions(n_y)
        n_x_vals = np.array([n_y])

        for n_x in n_x_vals:
            for betaJ in self.betaJ_vals:
                for r in reps:
                    path = 'results/' + self.prefix + f'/N_y={n_y}/N_x={n_x}/betaJ={betaJ}/rep{r}/data.csv'
                    df = pd.read_csv(path, sep='\s+')
                    t = df['iter'].to_numpy()
                    m_int = df['m_int'].to_numpy()
                    
                    fig, ax = plt.subplots()
                    ax.set_xlabel('Time steps (50 MC steps)')
                    ax.set_ylabel(r'$\langle m \rangle$')
                    ax.scatter(t[::100], m_int[::100])
                    
                    os.makedirs(file_path + f'/N_x={n_x}/betaJ={betaJ}/rep{r}', exist_ok=True)
                    fig.savefig(file_path + f'/N_x={n_x}/betaJ={betaJ}/rep{r}/histogram.pdf')
                    plt.close()

def generate_parameters(p):
    return int(p.min_y), int(p.max_y), int(p.y_step), int(p.max_x), int(p.x_step), str(p.prefix), float(p.Jmax), float(p.Jstep), int(reps)

def generate_couplings(Jmax, Jstep):
    couplings = np.arange(-Jmax, Jmax + Jstep, Jstep)
    couplings = np.where(np.isclose(couplings, -0.0, atol=1e-9), 0.0, couplings)
    return couplings

min_y, max_y, y_step, max_x, x_step, prefix, Jmax, Jstep, reps = generate_parameters(snakemake.params)
couplings = generate_couplings(Jmax, Jstep)
n_y_vals = np.arange(min_y, max_y + y_step, y_step)
reps = np.arange(1, reps + 1)


print(n_y_vals)
print(f"Max N_x is {max_x}\n")
print(f"N_x stepsize is {x_step}\n")
print(f"Prefix is {prefix}\n")
print(f"Max J is {Jmax}\n")
print(f"J step is {Jstep}\n")
print(f"Number of reps is {reps.size}\n")

data = Data(prefix, n_y_vals, max_x, couplings, x_step, reps)

data.calc_thermodynamic_quantities()

# Uncomment this block to also plot magnetisation curves for each value of N_y (square lattice only,
# this can be changed by editing make_magnetisation_plots)
#for n_y in data.n_y_vals:
#    data.make_magnetisation_plots(n_y)
#    plt.clf()
