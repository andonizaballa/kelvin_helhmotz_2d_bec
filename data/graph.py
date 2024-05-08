import matplotlib.pyplot as plt
import numpy as np
import re
import os 
import pandas as pd
import glob
from PIL import Image
import matplotlib.animation as animation

from fig_config import (
    add_grid,
    figure_features,
)  # <--- import customized functions

def create_images():
    #Get the files
    dens_files = glob.glob('density/dens-*.dat')
    phase_file = glob.glob('phase/phase-*.dat')
    vel1_files = glob.glob('vel1/vel1-*.dat')
    vel2_files = glob.glob('vel2/vel2-*.dat')
    os.makedirs('images', exist_ok=True) 
    for file in dens_files:
        plot_graph(file)
    
    for file in phase_file:
        plot_graph(file)
    
    for file in vel1_files:
       plot_graph(file)

    for file in vel2_files:
        plot_graph(file)




def extract_number(filename):
    # Regular expression to match the numerical part of the file name
    match = re.search(r'\d+', filename)
    if match:
        return int(match.group())
    else:
        return -1  # Return -1 if no number found

def make_gif():
    #Remove the files with no numbers in the file

    # Get the list of filenames matching the pattern sorted numerically
    dens_file_names = sorted(glob.glob("density/dens-*.dat.png"), key=extract_number)
    phase_file_names = sorted(glob.glob("phase/phase-*.dat.png"), key=extract_number)
    vel1_file_names = sorted(glob.glob("vel1/vel1-*.dat.png"), key=extract_number)
    vel2_file_names = sorted(glob.glob("vel2/vel2-*.dat.png"), key=extract_number)
    
    # Open images in sorted order
    dens_frames = [Image.open(image) for image in dens_file_names]
    phase_frames = [Image.open(image) for image in phase_file_names]
    vel1_frames = [Image.open(image) for image in vel_file_names]
    vel2_frames = [Image.open(image) for image in vel_file_names]


    print(vel_frames)

    frame_one = dens_frames[0]
    frame_one.save("density.gif", format="GIF", append_images=dens_frames, save_all=True, duration=150, loop=0)
    
    frame_one = phase_frames[0]
    frame_one.save("phase.gif", format="GIF", append_images=phase_frames, save_all=True, duration=150, loop=0)
    
    frame_one = vel1_frames[0]
    frame_one.save("velocity1.gif", format="GIF", append_images=vel1_frames,
                   save_all=True, duration=150, loop=0)
    
    frame_one = vel2_frames[0]
    frame_one.save("velocity2.gif", format="GIF", append_images=vel2_frames,
                   save_all=True, duration=150, loop=0)
    
    
    # Remove the images
    for file in dens_file_names:
        os.remove(file)

    for file in phase_file_names:
        os.remove(file)

    for file in vel1_file_names:
        os.remove(file)

    for file in vel2_file_names:
        os.remove(file)
    
def plot_graph(file_name):
    figure_features()

    #Open the file
    x , y , density  = np.loadtxt(file_name, unpack=True)
    
    #Plot the graph
    plt.figure()
    plt.scatter(x, y, c = density, vmin = 0, cmap = 'plasma', marker = 'o')
    plt.colorbar()
    plt.title(file_name)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    os.chdir('images')
    plt.savefig(file_name + '.png')

    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf() 

def plot_velocity(file):

    figure_features()

    x , y , vx  = np.loadtxt(file, unpack=True)

    # We want to plot y versus vx for a given value of x 

    for xnum in x:
        # Get the index of the x value
        
        if xnum == 17.8125 :
            index = np.where(x == xnum)
            xlabel = xnum
            
    plt.plot(y[index], vx[index], label = 'x = ' + str(xlabel))


    plt.ylabel('$v_x$')
    plt.xlabel('$ y $')
    plt.legend()
    plt.title('Velocity in x direction vs y')
    plt.savefig('velocity.png')
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()

def instablity_regime_df():

    os.chdir('FT-x')

    #Retrieve the data from FT-x
    ft_file_names = sorted(glob.glob("FT-x*.dat"), key=extract_number)



    # List to store vectors for each column
    column_names = ['t', 'px', 'psi_px']

    instablity_regime_df = pd.DataFrame(columns=column_names)

    for file in ft_file_names:
        input_file = open(file)

        lines = input_file.readlines()

        for line in lines:
            #Search for the values

            values = line.split()
            t = float(values[0])
            px = float(values[1])
            psi_px = float(values[2])
            instablity_regime_df.loc[len(instablity_regime_df.index)] = [t, px, psi_px]

    os.chdir('..')
    #Save the dataframe as a csv file

    instablity_regime_df.to_csv('instability_regime.csv')

def plot_instability_regime():

    # Read the csv file
    instablity_regime_df = pd.read_csv('instability_regime.csv')

    #Plot the graph
    figure_features()




    plt.figure()

    plt.scatter(instablity_regime_df['t'], abs(instablity_regime_df['psi_px'])**2, c = instablity_regime_df['px'] ,label = 'FT-x',s = 0.5,cmap = 'tab20')
 
    plt.xlabel('Time ($t$)')
    plt.ylabel(r'$\tilde{n} (p_x) $') 
    plt.yscale('log')
    #plt.xscale('log')
    cbar = plt.colorbar()
    cbar.set_label('$p_x$')
    #plt.legend(loc = 'upper right')
    plt.title('Dynamical (modulational) instability (FT-x)')
    plt.savefig('instability_regime.png')
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()




if __name__ == '__main__':
    create_images()
    make_gif()

    #instablity_regime_df()
    plot_instability_regime()


    






