import matplotlib.pyplot as plt
import numpy as np
import re
import os 
import pandas as pd
from scipy.optimize import curve_fit

from fig_config import (
    add_grid,
    figure_features,
)  # <--- import customized functions

def main():
    # Load the data
    instablity_regime_df = import_csv()

    # Fit the data

    px_fit = fit_all_the_px(instablity_regime_df)

    # Plot the data

    plot_fit(px_fit, instablity_regime_df)


def import_csv():
    instablity_regime_df = pd.read_csv('instability_regime.csv')
    return instablity_regime_df

def exponential(x, a, b):
    return a * np.exp(b * x)

def linear(x, a, b):
    return a * x + b

def fit_exponential_growth(px,instablity_regime_df):

    t0 = 15
    t1 = 35

    instablity_regime_df = instablity_regime_df[instablity_regime_df['t'] > t0]
    instablity_regime_df = instablity_regime_df[instablity_regime_df['t'] < t1]

    instablity_regime_df = instablity_regime_df[instablity_regime_df['px'] == px]

    t = instablity_regime_df['t']
    n_kx = instablity_regime_df['psi_px']

    popt, pcov = curve_fit(linear, t, np.log(n_kx))

    return popt

def fit_all_the_px(instablity_regime_df):

    px_fit = np.array([])

    for px in instablity_regime_df['px'].unique():
        popt = fit_exponential_growth(px,instablity_regime_df)
        # Create a vector with px and the fitted values
        # if popt[1] is 1 or 1 by 10-4 then the fit is not good

        #if popt[1] < 1 + 1e-4 and popt[1] > 1 - 1e-4:
            #continue
        px_fit = np.append(px_fit, [px, popt[0], popt[1]])

    return px_fit

def plot_fit(px_fit,instablity_regime_df):

    figure_features()

    # We will take tha a of 'px=0' as the reference value

    popt = fit_exponential_growth(0,instablity_regime_df)

    px = px_fit[0::3]
    a = px_fit[1::3]
    b = px_fit[2::3]

    m = px * 60 / 2 / np.pi

    delta_w = 10 

    a_max = np.max(a)


    plt.plot(px, a/a_max, 'o')
    plt.xlabel('$k_x$')
    plt.ylabel('$\sigma_{k_x}$')
    # Save the plot

    plt.savefig('fit.png', dpi = 300)

    #Close the plot
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()











if __name__ == "__main__":
    main()