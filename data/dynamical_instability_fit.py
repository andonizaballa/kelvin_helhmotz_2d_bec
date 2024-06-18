import matplotlib.pyplot as plt
import numpy as np
import re
import os 
import pandas as pd
from scipy.optimize import curve_fit
import scipy.constants as sc
from scipy.constants import hbar
import cmath
from scipy.constants import hbar

from fig_config import (
    add_grid,
    figure_features,
)  # <--- import customized functions

def main():
    # Load the data
    instablity_regime_df = import_csv('instability_regime')
    time_df = import_csv('time_vector')

    # First we will show the graph that we want to fit the exponential growth

    p_x_fit_df = plot_instability_allpx(instablity_regime_df,time_df)

    delta_fit = plot_velocity_profile(0,'vel2-015.dat')

    # Now we will plot the fit

    plot_fit(p_x_fit_df,instablity_regime_df, delta_fit)

    #plot_khi(instablity_regime_df, delta_fit)


def import_csv(name):
    df = pd.read_csv(str(name)+'.csv')
    return df

def exponential(x, a, b):
    return a * np.exp(b * x)

def linear(x, a, b):
    return a * x + b



def fit_exponential_growth(t0,t1,px,instablity_regime_df):

    instablity_regime_df = instablity_regime_df[instablity_regime_df['t'] >= t0]
    instablity_regime_df = instablity_regime_df[instablity_regime_df['t'] <= t1]

    instablity_regime_df = instablity_regime_df[instablity_regime_df['px'] == px]

    t = instablity_regime_df['t']
    n_kx = instablity_regime_df['psi_px']

    popt, pcov = curve_fit(linear, t, np.log(n_kx))

    return popt, pcov

def fit_all_the_px(time,instablity_regime_df):

    px_fit_df = pd.DataFrame(columns = ['px', 'a', 'b', 'psi_t60'])

    for px in p_x:
        popt = fit_exponential_growth(t0,t1,px,instablity_regime_df)
        # Create a vector with px and the fitted values
        # if popt[1] is 1 or 1 by 10-4 then the fit is not good

        #if popt[1] < 1 + 1e-4 and popt[1] > 1 - 1e-4:
            #continue

        # We will the value of psi_px at t=60 from the dataframe

        psi_t60 = instablity_regime_df[instablity_regime_df['px'] == px]['psi_px'].values[-1]

        px_fit_df = px_fit_df.append({'px': px, 'a': popt[0], 'b': popt[1], 'psi_t60': psi_t60}, ignore_index = True)
        print(px_fit_df)
    return px_fit_df

def plot_fit(px_fit_df,instablity_regime_df, delta_fit):

    figure_features()

    # We will take tha a
    #a_max = fit_exponential_growth(0,instablity_regime_df)[0]

    #only positive or cero values of px

    px_fit_df = px_fit_df[px_fit_df['px'] >= 0]

    #We will divide the values in px_fit_df in three categories low, medium and high
    px_low_df = px_fit_df[px_fit_df['psi_t60'] < 1e-5]
    px_medium_df = px_fit_df[(px_fit_df['psi_t60'] >= 1e-5) & (px_fit_df['psi_t60'] < 1e-1)]
    px_high_df = px_fit_df[px_fit_df['psi_t60'] >= 1e-1]

    #Now we plot them. We will make three different plots, one for each category

    #Low values of psi_t60

    fig, axs = plt.subplots(4, 1, figsize=(8, 10), sharey=True)

    #now we will plot the points with the error bars

    axs[0].errorbar(px_low_df['px'], px_low_df['a'], yerr = px_low_df['error_a'], fmt = 'o', label = 'Low values of $\psi_{t=60}$', c = 'steelblue', ecolor='lightgray', elinewidth=3, markersize = 5)
    axs[1].errorbar(px_medium_df['px'], px_medium_df['a'], yerr = px_medium_df['error_a'], fmt = 'o', label = 'Medium values of $\psi_{t=60}$', c = 'orangered', ecolor='lightgray', elinewidth=3, markersize = 5)
    axs[2].errorbar(px_high_df['px'], px_high_df['a'], yerr = px_high_df['error_a'], fmt = 'o', label = 'High values of $\psi_{t=60}$',c = 'olivedrab', ecolor='lightgray', elinewidth=3, markersize = 5)

    axs[3].errorbar(px_low_df['px'], px_low_df['a'], yerr = px_low_df['error_a'], fmt = 'o', c = 'steelblue', ecolor='lightgray', elinewidth=3, markersize = 5)
    axs[3].errorbar(px_medium_df['px'], px_medium_df['a'], yerr = px_medium_df['error_a'], fmt = 'o', c = 'orangered', ecolor='lightgray', elinewidth=3, markersize = 5)
    axs[3].errorbar(px_high_df['px'], px_high_df['a'], yerr = px_high_df['error_a'], fmt = 'o',c = 'olivedrab', ecolor='lightgray', elinewidth=3, markersize = 5)
    # Labels fro the x and y axis

    axs[0].set_ylabel('$\sigma_m$' , rotation = 0, labelpad = 20)
    axs[1].set_ylabel('$\sigma_m$', rotation = 0, labelpad = 20)
    axs[2].set_ylabel('$\sigma_m$', rotation = 0, labelpad = 20)
    axs[3].set_ylabel('$\sigma_m$', rotation = 0, labelpad = 20)

    axs[3].set_xlabel('$k_x$')

    #Legends 

    axs[0].legend(frameon = False, loc = 'upper right')
    axs[1].legend(frameon = False, loc = 'upper right')
    axs[2].legend(frameon = False, loc = 'upper right')

    # Save the figure in this folder

    fig.savefig('fit_classified.png', dpi = 300)

    #Close the plot
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()

    # Now we plot all the data in the same graph

    figure_features()

    fig = plt.figure( figsize=(8, 10))

    # Aspect ratio

    aspect_ratio = 1.6168
    fig.set_size_inches(8, 8 / aspect_ratio)

    #plt.errorbar(px_low_df['px'], px_low_df['a'], yerr = px_low_df['error_a'], fmt = 'o', label = 'Low values of $\psi_{t=60}$', c = 'steelblue', ecolor='lightgray', elinewidth=3, markersize = 5)
    #plt.errorbar(px_medium_df['px'], px_medium_df['a'], yerr = px_medium_df['error_a'], fmt = 'o', label = 'Medium values of $\psi_{t=60}$', c = 'orangered', ecolor='lightgray', elinewidth=3, markersize = 5)
    #plt.errorbar(px_high_df['px'], px_high_df['a'], yerr = px_high_df['error_a'], fmt = 'o', label = 'High values of $\psi_{t=60}$',c = 'olivedrab', ecolor='lightgray', elinewidth=3, markersize = 5)
    plt.errorbar(px_fit_df['px'], px_fit_df['a'], yerr = px_fit_df['error_a'], fmt = 'o', label = 'All values of $\psi_{t=60}$',c = 'crimson', ecolor='lightgray', elinewidth=3, markersize = 5)
    plt.ylabel('$\sigma_m$' , rotation = 0, labelpad = 20)
    plt.xlabel('$k_x$')

    #plt.legend(frameon = False, loc = 'upper right')

    plt.savefig('fit_all.png', dpi = 300)

    #Close the plot
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()


def im_omega_kh(delta, kx):
    return  1/(2*delta) * np.sqrt (np.exp(-4 * kx * delta ) - (2 * kx * delta - 1)**2) 



def extract_number(filename):
    # Regular expression to match the numerical part of the file name
    match = re.search(r'-([0-9]+)', filename)

    
    if match:
        return (-1)*int(match.group())
    else:
        return -1  # Return -1 if no number found

def tanh_vel(x,deltav,delta):
    return deltav / 2.0 * np.tanh(x / delta)

def plot_velocity_profile(xprof,file):

    os.chdir('vel2')

    figure_features()

    x , y , vx , vy  = np.loadtxt(file, unpack=True)

    time = extract_number(file)

    #print(time)

    # We want to plot y versus vx for a given value of x 

    for xnum in x:
        # Get the index of the x value
        
        if xnum == xprof :
            index = np.where(x == xnum)
            xlabel = xnum

    popt , pcov = curve_fit(tanh_vel, y[index], vx[index])
            
    plt.plot(y[index], vx[index], label = 'x = ' + str(xlabel))
    plt.plot(y[index], tanh_vel(y[index],popt[0],popt[1]), '--', label = r'$\delta = $' + str(popt[1] ))


    plt.ylabel('$v_x$')
    plt.xlabel('$ y $')
    plt.title('t = ' + str(time))
    plt.legend()
    plt.xlim(-17.5,17.5)
    plt.ylim(-0.5E-8,0.5E-8)
    #plt.title('Velocity in x direction vs y')
    plt.savefig('delta_fit.png', dpi = 300)
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()

    os.chdir('..')

    return popt

def plot_instability_regime_px(t0,t1,df, px):

    # Read the csv file
    instablity_regime_df = df

    instablity_regime_df = instablity_regime_df.sort_values(by = 't')

    #Plot the graph
    figure_features()

    fig = plt.figure( figsize=(8, 10))

    # Aspect ratio

    aspect_ratio = 1.6168
    fig.set_size_inches(8, 8 / aspect_ratio)

    #For every px value, plot the graph

    # # CWe assign a different color to each px value
    # color_map = plt.cm.get_cmap('tab20')
    # for px in instablity_regime_df['px'].unique():
    #     px_df = instablity_regime_df[instablity_regime_df['px'] == px]
    #     plt.scatter(px_df['t'], abs(px_df['psi_px'])**2, label='px = ' + str(px), s=0.5, c= 'tab20')
    


    popt, pcov = fit_exponential_growth(t0,t1,px,instablity_regime_df)

    instablity_regime_df = instablity_regime_df[instablity_regime_df['px'] == px]
    plt.plot(instablity_regime_df['t'], instablity_regime_df['psi_px'], linewidth = 1,  c = 'orangered' )

    instablity_regime_df = instablity_regime_df[instablity_regime_df['t'] >= t0]
    instablity_regime_df = instablity_regime_df[instablity_regime_df['t'] <= t1]
    a = popt[0]
    b = np.exp(popt[1])
    plt.plot(instablity_regime_df['t'], np.exp(popt[1] + popt[0] * instablity_regime_df['t']), '--', c = 'blue', lw = 1, label='$Be^{\sigma t}$ \n $B$ = %.4e \n $\sigma $ = %.5f ' %  (b,a))
    #Two vertical lineas at t0 and t1

    plt.axvline(x = t0, color = 'black', linestyle = '--', lw = 0.5)
    plt.axvline(x = t1, color = 'black', linestyle = '--', lw = 0.5)
 
    plt.xlabel('Time ($t$)')
    plt.xlim(0, 60)
    plt.ylabel(r'$\tilde{n} (\bf{k_x}) $', rotation = 0, labelpad = 20) 
    plt.yscale('log')
    #plt.xscale('log')
    #plt.legend(loc = 'upper right')
    plt.title(r'$\bf{k_x}$ = ' + str(px))
    plt.legend(loc = 'lower right')
    os.chdir('instability_regime')
    plt.savefig('instability_regime_px_' + str(px) + '.png', dpi = 300)
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()

    os.chdir('..')

    return popt, pcov
    
def plot_instability_allpx(instablity_regime_df, time_df):

    p_x_fit_df = pd.DataFrame(columns = ['px', 'a', 'b', 'psi_t60', 'error_a', 'error_b'])
    for px in instablity_regime_df['px'].unique():
        t0 = time_df[time_df['px'] == px]['t0'].values[0]
        t1 = time_df[time_df['px'] == px]['t1'].values[0]
        popt, pcov = plot_instability_regime_px(t0,t1,instablity_regime_df, px)
        psi_t60 = instablity_regime_df[instablity_regime_df['px'] == px]['psi_px'].values[-1]
        errora = np.sqrt(pcov[0][0])
        errorb = np.sqrt(pcov[1][1])
        p_x_fit_df.loc[len(p_x_fit_df.index)] = [px, popt[0], popt[1], psi_t60, errora, errorb]
    return p_x_fit_df


def plot_khi(instability_regime_df, delta_fit):
    a_mu = sc.physical_constants['atomic mass constant'][0]
    m = 12 * a_mu
    v_lab = delta_fit[0]


    k_h = m * delta_fit[0] / hbar - instability_regime_df['px'].unique()
    plt.plot(k_h, im_omega_kh(delta_fit[1], delta_fit[0],k_h), linewidth = 1.5,  c = 'olive', label = r'$\delta = $' + str(delta_fit[1]) + ' \n  $v_{lab} = $' + str(v_lab))
    plt.legend()
    plt.savefig('khi.png', dpi = 300)
    #plt.show()


    #Close the plot
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()
if __name__ == "__main__":
    main()