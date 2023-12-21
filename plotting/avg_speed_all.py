import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolor
import os
import sys
import math
import matplotlib.patheffects as pe

def get_parameters_from_filename(file):
    params = file.split("_")
    
    relax_time = int(params[1][1:len(params[1])-4])

    return relax_time

def gather_data(directory):
    V = dict()

    for file in os.listdir(directory):
        with open(os.path.join(directory, file),"r") as f:
            relax_duration_time = get_parameters_from_filename(file)
            data = np.loadtxt(f, unpack=True)
            Vx = data[2]
            time = data[0]
            V[relax_duration_time] = Vx

    V = dict(sorted(V.items()))

    return V, time

def plot(file):
    V, time = gather_data(file)

    colors = ['b', 'c', 'g', 'k', 'm', 'r', 'y', 'orange', 'brown', 'navy', 'hotpink', 'lime']
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15,7))
    
    i = 0
    for relax_duration, v in V.items():
        ax.plot(time, v, colors[i], linewidth=1) 
        i += 1

    ax.tick_params(axis='both', labelsize=15)
    ax.set_xlim(xmin=0, xmax=1*10**7)
    ax.set_ylim(ymin=-0.05)
    ax.set_xlabel(rf'$\tau$', fontsize=20)
    ax.set_ylabel(rf'$\langle V_x \rangle$', fontsize=20)

    plt.show()


if __name__ == '__main__':
    if not os.path.isdir(sys.argv[1]):
        print(f"{sys.argv[1]} not a directory.")
        sys.exit(-1)
    plot(sys.argv[1])