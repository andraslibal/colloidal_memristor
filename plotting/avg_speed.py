import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolor
import os
import sys
import math
import matplotlib.patheffects as pe

def get_parameters_from_filename(file):
    params = file.split("_")
    
    af = float(params[0][1:len(params[0])])
    N = int(params[1][1:len(params[1])])
    Np = int(params[2][2:len(params[2])])
    Fp = params[3][2:len(params[3])]

    if Fp[len(Fp)-4:len(Fp)] == ".txt":
        Fp = int(Fp[2:len(Fp)-4])
        Fd = 1
    else:
        Fp = int(Fp)
        Fd = float(params[4][2:len(params[4])-4])

    return af, N, Np, Fp, Fd

def gather_data(file):
    Vx = []

    with open(file,"r") as f:
        data = np.loadtxt(f, unpack=True)
        Vx = data[2]
        time = data[0]
    return time, Vx

def plot(file):
    time, Vx = gather_data(file)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(15,7))
    ax.plot(time, Vx, color='red', linewidth=1)
    ax.tick_params(axis='both', labelsize=15)
    ax.set_xlim(xmin=0, xmax=10**7)
    ax.set_ylim(ymin=-0.1)
    ax.set_xlabel(rf'$\tau$', fontsize=20)
    ax.set_ylabel(rf'$\langle V_x \rangle$', fontsize=20)

    plt.show()

if __name__ == '__main__':
    if not os.path.isfile(sys.argv[1]):
        print(f"{sys.argv[1]} not a file.")
        sys.exit(-1)
    plot(sys.argv[1])

