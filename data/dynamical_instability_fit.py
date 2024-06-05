import matplotlib.pyplot as plt
import numpy as np
import re
import os 
import pandas as pd
from scipy.optimize import curve_fit
import scipy.constants as sc
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

    p_x_fit = plot_instability_allpx(instablity_regime_df,time_df)

    delta_fit = plot_velocity_profile(0,'vel2-015.dat')

    # Now we will fit the exponential growth for each px value

    #px_fit = fit_all_the_px(time,instablity_regime_df)

    # Now we will plot the fit

    plot_fit(p_x_fit,instablity_regime_df, delta_fit)


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

    return popt

def fit_all_the_px(time,instablity_regime_df):

    px_fit = np.array([])

    p_x = time[0::3]
    t0 = time[1::3]
    t1 = time[2::3]

    for px in p_x:
        popt = fit_exponential_growth(t0,t1,px,instablity_regime_df)
        # Create a vector with px and the fitted values
        # if popt[1] is 1 or 1 by 10-4 then the fit is not good

        #if popt[1] < 1 + 1e-4 and popt[1] > 1 - 1e-4:
            #continue
        px_fit = np.append(px_fit, [px, popt[0], popt[1]])

    return px_fit

def plot_fit(px_fit,instablity_regime_df, delta_fit):

    figure_features()

    # We will take tha a
    #a_max = fit_exponential_growth(0,instablity_regime_df)[0]

    px = px_fit[0::3]
    a = px_fit[1::3]
    b = px_fit[2::3]


    #plt.plot(px, a, 'o', color = 'darkgreen')
    omega_t = np.imag(im_omega_kh(delta_fit[1], delta_fit[0], px))
    plt.plot(px, omega_t, '--')
    plt.xlabel('$k_x$')
    plt.ylabel(r'$\sigma_{m} \over \sigma^*$', rotation = 0, labelpad = 20)
    # Save the plot

    plt.savefig('fit.png', dpi = 300)

    #Close the plot
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()

def im_omega_kh(delta, deltav, k):
    vec = np.array([])
    for kx in k:
        vec = np.append(vec, deltav / (4 * delta) * cmath.sqrt(np.exp(-4 * kx * delta) - (2 * kx * delta - 1)**2))
    return  vec

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
    


    popt = fit_exponential_growth(t0,t1,px,instablity_regime_df)

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

    return popt
    
def plot_instability_allpx(instablity_regime_df, time_df):

    p_x_fit = np.array([])
    for px in instablity_regime_df['px'].unique():
        t0 = time_df[time_df['px'] == px]['t0'].values[0]
        t1 = time_df[time_df['px'] == px]['t1'].values[0]
        popt = plot_instability_regime_px(t0,t1,instablity_regime_df, px)
        p_x_fit = np.append(p_x_fit, [px, popt[0], popt[1]])
    return p_x_fit




if __name__ == "__main__":
    main()