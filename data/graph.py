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
    dens_files = glob.glob('dens-*.dat')
    phase_file = glob.glob('phase-*.dat')
    vel_files = glob.glob('vel1-*.dat')
    os.makedirs('images', exist_ok=True) 
    #for file in dens_files:
        #plot_graph(file)
    
    #for file in phase_file:
        #plot_graph(file)
    
    for file in vel_files:
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
    dens_file_names = sorted(glob.glob("dens-*.dat.png"), key=extract_number)
    phase_file_names = sorted(glob.glob("phase-*.dat.png"), key=extract_number)
    vel_file_names = sorted(glob.glob("vel1-*.dat.png"), key=extract_number)
    
    # Open images in sorted order
    dens_frames = [Image.open(image) for image in dens_file_names]
    phase_frames = [Image.open(image) for image in phase_file_names]
    vel_frames = [Image.open(image) for image in vel_file_names]

    print(vel_frames)

    frame_one = dens_frames[0]
    frame_one.save("density.gif", format="GIF", append_images=dens_frames, save_all=True, duration=150, loop=0)
    
    frame_one = phase_frames[0]
    frame_one.save("phase.gif", format="GIF", append_images=phase_frames, save_all=True, duration=150, loop=0)
    
    frame_one = vel_frames[0]
    frame_one.save("velocity.gif", format="GIF", append_images=vel_frames,
                   save_all=True, duration=150, loop=0)
    
    
    # Remove the images
    for file in dens_file_names:
        os.remove(file)

    for file in phase_file_names:
        os.remove(file)

    for file in vel_file_names:
        os.remove(file)
    
def plot_graph(file_name):
    figure_features()

    #Open the file
    input_file = open(file_name)

    lines = input_file.readlines()

    #Create a dataframe
    column_names = ['x', 'y', 'density']

    density_df = pd.DataFrame(columns=column_names)

    for line in lines:
        #Search for the values

        search = re.findall(r'[-+]?\d*\.\d+(?:E[-+]?\d+)?', line)
                            
        if re.search(r'[-+]?\d*\.\d+(?:E[-+]?\d+)?', line) is None:
            continue

        x = float(search[0])
        y = float(search[1])
        density = float(search[2])
        density_df.loc[len(density_df.index)] = [x, y, density]
    
    #Plot the graph
    plt.figure()
    plt.scatter(density_df['x'], density_df['y'], c = density_df['density'], vmin = 0, cmap = 'plasma')
    plt.colorbar()
    plt.title(file_name)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
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



if __name__ == '__main__':
    #create_images()
    #make_gif()

    plot_velocity('vel2-053.dat')


    






