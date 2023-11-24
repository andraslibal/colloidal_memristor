import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolor
import os
import sys
import math
import matplotlib.patheffects as pe

relax_start = 15000000
stat_time = 1000

def get_parameters_from_filename(file):
    params = file.split("_")
    
    relax_time = int(params[1][1:len(params[1])-4])

    return relax_time

def gather_data(directory):
    V = dict()
    Vx = []

    for file in os.listdir(directory):
        with open(os.path.join(directory, file),"r") as f:
            relax_time = get_parameters_from_filename(file)
            data = np.loadtxt(f, unpack=True)
            Vx = data[2]
            time = data[0]

            start_speed = 0
            i = 0
            while time[i] < relax_start:
                i += 1
            start_index = i + relax_time//stat_time - 1
            end_index = start_index + 2500
            start_speed = Vx[start_index : end_index]

            print(f"start {start_index}, time {time[i]}, end {end_index}")
            print(start_speed)

            V[relax_time] = sum(start_speed)/len(start_speed)
    
    V = dict(sorted(V.items()))

    return V
    
def plot(directory):
    data = gather_data(directory)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,7))
    ax.plot(data.keys(), data.values(), color='blue', linewidth=2, marker='o')
    ax.tick_params(axis='both', labelsize=15)
    #ax.set_xlim(xmin=0, xmax=5*10**7)
    #ax.set_ylim(ymin=-0.1)
    ax.set_xlabel(rf'$\tau$', fontsize=20)
    ax.set_ylabel(rf'$\langle V_x \rangle$', fontsize=20)

    plt.show()

if __name__ == '__main__':
    if not os.path.isdir(sys.argv[1]):
        print(f"{sys.argv[1]} not a directory.")
        sys.exit(-1)
    plot(sys.argv[1])

