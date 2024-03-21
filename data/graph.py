import matplotlib.pyplot as plt
import numpy as np
import re
import os 
import pandas as pd
import glob
from PIL import Image
import matplotlib.animation as animation

def create_images():
    #Get the files
    files = glob.glob('dens-*.dat')
    os.makedirs('images', exist_ok=True)
    for file in files:
        print(file)
        plot_graph(file)



def extract_number(filename):
    # Regular expression to match the numerical part of the file name
    match = re.search(r'\d+', filename)
    if match:
        return int(match.group())
    else:
        return -1  # Return -1 if no number found

def make_gif():
    # Get the list of filenames matching the pattern sorted numerically
    file_names = sorted(glob.glob("dens-*.dat.png"), key=extract_number)
    
    # Open images in sorted order
    frames = [Image.open(image) for image in file_names]

    frame_one = frames[0]
    frame_one.save("density.gif", format="GIF", append_images=frames,
                   save_all=True, duration=150, loop=0)
    
    # Remove the images
    for file in file_names:
        os.remove(file)
    
def plot_graph(file_name):

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
    plt.tricontourf(density_df['x'], density_df['y'], density_df['density'])
    plt.colorbar()
    plt.title(file_name)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig(file_name + '.png')

    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf() 

if __name__ == '__main__':
    create_images()
    make_gif()


    






