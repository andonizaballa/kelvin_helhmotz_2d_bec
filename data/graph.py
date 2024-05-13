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
    #create_images('vel2')

    #Create the images for the velocity profile
    plot_velocity_profile_allx()

    #Create the gif for density, phase, vel1 and vel2
    #make_gif('dens')
    #make_gif('phase')
    #make_gif('vel2')

    #Create the gif for the velocity profile
    #gif_velocity_profile()

    #Create the dataframe for the instability regime   
    #instablity_regime_df()

    #Plot the instability regime
    #plot_instability_regime()

    #Plot the instability regime for all px values
    #plot_instability_allpx()

    


def create_images(folder_name):

    dict_folder = {'dens': 'density', 'phase': 'phase', 'vel1': 'vel1', 'vel2': 'vel2'}

    dict_max = {'dens': 9E-4, 'phase': 2*np.pi, 'vel1': 0.5E-88, 'vel2': 2.5}
    dict_min = {'dens': 0, 'phase': 0, 'vel1': -0.5E-8, 'vel2': -2.5}

    max = dict_max[folder_name]
    min = dict_min[folder_name]

    #Get the files

    os.chdir(dict_folder[folder_name])

    files = sorted(glob.glob(folder_name + "-*.dat"), key=extract_number)

    #print(files)

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
    match = re.search(r'-([0-9]+)', filename)

    
    if match:
        return (-1)*int(match.group())
    else:
        return -1  # Return -1 if no number found

def make_gif(folder_name):
    dict_folder = {'dens': 'density', 'phase': 'phase', 'vel1': 'vel1', 'vel2': 'vel2'}

    os.chdir(dict_folder[folder_name])

    # Get the list of filenames matching the pattern sorted numerically
    file_names = sorted(glob.glob(folder_name + "-*.dat.png"), key=extract_number)

    
    # Open images in sorted order
    frames = [Image.open(image) for image in file_names]




    frame_one = frames[0]
    frame_one.save(folder_name  + ".gif", format="GIF", append_images=frames, save_all=True, duration=100, loop=0)
    
    
    
    # Remove the images
    #for file in file_names:
     #   os.remove(file)
    os.chdir('..')
    
def plot_graph_dens(file_name,max,min):
    figure_features()

    #Open the file
    x , y , density  = np.loadtxt(file_name, unpack=True)

    time = extract_number(file_name)
    
    #Plot the graph
    fig = plt.figure()

    aspect_ratio = 2
    fig.set_size_inches(8, 8 / aspect_ratio)
    plt.scatter(x, y, c = density, vmin = min, vmax = max,  cmap = 'plasma', marker = 'o')
    plt.colorbar()
    plt.xlim(-30,30)
    plt.ylim(-20,20)
    plt.title('t = ' + str(time))
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

def plot_velocity_profile(xprof,file):

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
            
    plt.plot(y[index], vx[index], label = 'x = ' + str(xlabel))


    plt.ylabel('$v_x$')
    plt.xlabel('$ y $')
    plt.title('t = ' + str(time))
    plt.legend()
    plt.xlim(-17.5,17.5)
    plt.ylim(-2.5,2.5)
    #plt.title('Velocity in x direction vs y')
    plt.savefig('velocity_profile_x_'+ str(xprof)+'_'+str(file)+'.png')
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()

def plot_velocity_profile_allx():

    os.chdir('vel2')

    file_names = sorted(glob.glob('vel2-*.dat'), key=extract_number)
                                  
    for file in file_names:
        plot_velocity_profile(17.8125,file)

    os.chdir('..')




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
    

    y_lim =max(instablity_regime_df['px'])

    # We will delete the px values +- 0.5236

    instablity_regime_df = instablity_regime_df[instablity_regime_df['px'] != 0.5236]

    instablity_regime_df = instablity_regime_df[instablity_regime_df['px'] != -0.5236]

    

    #print(instablity_regime_df['px'])


    plt.scatter(instablity_regime_df['t'], instablity_regime_df['px'] ,c = abs(instablity_regime_df['psi_px'])**2 ,label = 'FT-x',s = 9,cmap = 'plasma', vmin= 1e7, vmax = 5e9)
 
    plt.xlabel('Time ($t$)')
    plt.xlim(0, 200)
    plt.ylim(-y_lim, y_lim)
    plt.ylabel('$p_x$', rotation = 0, y = 0.45) 
    #plt.yscale('log')
    #plt.xscale('log')
    cbar = plt.colorbar()
    cbar.set_label(r'$\tilde{n} (p_x) $', rotation = 0)
    #plt.legend(loc = 'upper right')
    #plt.title('Dynamical (modulational) instability (FT-x)')
    plt.savefig('instability_regime_t_vs_px.png')
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()

def plot_instability_regime_px(df, px):

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
    




    instablity_regime_df = instablity_regime_df[instablity_regime_df['px'] == px]
    plt.plot(instablity_regime_df['t'], instablity_regime_df['psi_px'], linewidth = 1,  c = 'orangered' )
 
    plt.xlabel('Time ($t$)')
    plt.xlim(0, 200)
    plt.ylabel(r'$\tilde{n} (\bf{x}) $', rotation = 0, labelpad = 20) 
    plt.yscale('log')
    #plt.xscale('log')
    #plt.legend(loc = 'upper right')
    #plt.title('Dynamical (modulational) instability (FT-x)')
    plt.savefig('instability_regime_px' + str(px) + '.png', dpi = 300)
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()



def plot_instability_allpx():

    instablity_regime_df = pd.read_csv('instability_regime.csv')

    os.chdir('instability_regime')

    for px in instablity_regime_df['px'].unique():
        plot_instability_regime_px(instablity_regime_df, px)

    os.chdir('..')

def gif_velocity_profile():

    os.chdir('vel2')

    # Get the list of filenames matching the pattern sorted numerically
    file_names = sorted(glob.glob('velocity_profile_x*.png'), key=extract_number)

    # Open images in sorted order
    frames = [Image.open(image) for image in file_names]

    frame_one = frames[0]
    frame_one.save('velocity_profile.gif', format="GIF", append_images=frames, save_all=True, duration=100, loop=0)
    
    os.chdir('..')




if __name__ == '__main__':
    main() 



    






