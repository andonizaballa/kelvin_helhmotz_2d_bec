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

def main():

    #Create the images for density, phase, vel1 and vel2
    #create_images('dens')
    #create_images('phase')
    create_images('vel1')
    #create_images('vel2')

    #Create the gif for density, phase, vel1 and vel2
    #make_gif('dens')
    #make_gif('phase')
    #make_gif('vel1')
    #make_gif('vel2')

def create_images(folder_name):

    dict_folder = {'dens': 'density', 'phase': 'phase', 'vel1': 'vel1', 'vel2': 'vel2'}

    dict_max = {'dens': 7E-4, 'phase': 2*np.pi, 'vel1': 0.5E-88, 'vel2': 2.5}
    dict_min = {'dens': 0, 'phase': 0, 'vel1': -0.5E-8, 'vel2': -2.5}

    max = dict_max[folder_name]
    min = dict_min[folder_name]

    #Get the files

    os.chdir(dict_folder[folder_name])

    files = sorted(glob.glob(folder_name + "-*.dat"), key=extract_number)

    print(files)

    os.makedirs('images', exist_ok=True)

    #If folder is density or phase

    if folder_name == 'dens' or folder_name == 'phase':
        for file in files:
            plot_graph_dens(file,max,min) 
    else: 
        for file in files:
            plot_graph_vel(file,max,min)

    os.chdir('..')
    




def extract_number(filename):
    # Regular expression to match the numerical part of the file name
    match = re.search(r'\d+', filename)
    if match:
        return int(match.group())
    else:
        return -1  # Return -1 if no number found

def make_gif(folder_name):
    dict_folder = {'dens': 'density', 'phase': 'phase', 'vel1': 'vel1', 'vel2': 'vel2'}

    os.chdir(dict_folder[folder_name])

    # Get the list of filenames matching the pattern sorted numerically
    file_names = sorted(glob.glob(folder_name + "-*.dat.png"), key=extract_number(folder_name + "-*.dat.png"))

    
    # Open images in sorted order
    frames = [Image.open(image) for image in file_names]




    frame_one = frames[0]
    frame_one.save(folder_name  + ".gif", format="GIF", append_images=frames, save_all=True, duration=150, loop=0)
    
    
    
    # Remove the images
    for file in file_names:
        os.remove(file)
    os.chdir('..')
    
def plot_graph_dens(file_name,max,min):
    figure_features()

    #Open the file
    x , y , density  = np.loadtxt(file_name, unpack=True)
    
    #Plot the graph
    fig = plt.figure()

    aspect_ratio = 2
    fig.set_size_inches(8, 8 / aspect_ratio)
    plt.scatter(x, y, c = density, vmin = min, vmax = max,  cmap = 'plasma', marker = 'o')
    plt.colorbar()
    #plt.title(file_name)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.savefig(file_name + '.png')

    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf() 

def plot_graph_vel(file_name,max,min):
    figure_features()

    fig = plt.figure()

    aspect_ratio = 2
    fig.set_size_inches(8, 8 / aspect_ratio)

    #Open the file
    x , y , vx, vy  = np.loadtxt(file_name, unpack=True)
    
    #Plot the graph

    plt.scatter(x, y, c = vx, cmap = 'viridis', vmin = min ,vmax=  max , marker = 'o')
    cbar = plt.colorbar()
    #plt.title(file_name)
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.xlim(-30,30)
    plt.ylim(-17.5,17.5)
    plt.savefig(file_name + '.png')
    cbar.set_label('$v_x$')

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
    plt.xlim(0,100)
    plt.ylim()
    #plt.title('Velocity in x direction vs y')
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
    #Order by px

    instablity_regime_df = instablity_regime_df.sort_values(by = 'px')

    instablity_regime_df.to_csv('instability_regime.csv')

def plot_instability_regime():

    # Read the csv file
    instablity_regime_df = pd.read_csv('instability_regime.csv')

    #Plot the graph
    figure_features()




    fig = plt.figure()

    # Aspect ratio

    aspect_ratio = 1.618
    fig.set_size_inches(8, 8 / aspect_ratio)

    #For every px value, plot the graph

    # # CWe assign a different color to each px value
    # color_map = plt.cm.get_cmap('tab20')
    # for px in instablity_regime_df['px'].unique():
    #     px_df = instablity_regime_df[instablity_regime_df['px'] == px]
    #     plt.scatter(px_df['t'], abs(px_df['psi_px'])**2, label='px = ' + str(px), s=0.5, c= 'tab20')
    


    plt.scatter(instablity_regime_df['t'], instablity_regime_df['px'] ,c = abs(instablity_regime_df['psi_px'])**2 ,label = 'FT-x',s = 2,cmap = 'plasma', vmin = 1E3, vmax = 1E5)
 
    plt.xlabel('Time ($t$)')
    plt.xlim(0, 100)
    plt.ylabel('$p_x$') 
    #plt.yscale('log')
    #plt.xscale('log')
    cbar = plt.colorbar()
    cbar.set_label(r'$\tilde{n} (p_x) $')
    #plt.legend(loc = 'upper right')
    #plt.title('Dynamical (modulational) instability (FT-x)')
    plt.savefig('instability_regime.png')
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()




if __name__ == '__main__':
    #create_images()
    #make_gif()
    #main()

    #instablity_regime_df()
    plot_instability_regime()


    






