# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last updated: 22-Oct-2023

import numpy as np
import os
import subprocess as sp
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import pandas as pd
import multiprocessing as mp
from functools import partial
import sys
import argparse


matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True

def execute_process(exe):
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    return stderr.decode("utf-8").split("\n")

def get_segs(place):
    temp2 = execute_process(["./getFacets", place])
    segs = []
    skip = False
    temp2 = list(filter(None, temp2))
    for n1 in range(len(temp2)):
        temp3 = temp2[n1].split(" ")
        if not skip:
            temp4 = temp2[n1+1].split(" ")
            x1, y1 = map(float, [temp3[0], temp3[1]])
            x2, y2 = map(float, [temp4[0], temp4[1]])
            segs.extend([((x1, y1),(x2, y2))])
            skip = True
        else:
            skip = False
    return segs

def get_field_values(place, xmin, xmax, ymin, ymax, ny):
    temp2 = list(filter(None, execute_process(["./getData", place, str(xmin), str(ymin), str(xmax), str(ymax), str(ny)])))
    data = np.array([line.split() for line in temp2], dtype=float)
    nx = data.shape[0] // ny
    X = data[:,0].reshape((nx, ny)).transpose()
    Y = data[:,1].reshape((nx, ny)).transpose()
    T = data[:,2].reshape((nx, ny)).transpose()
    D2 = data[:,3].reshape((nx, ny)).transpose()
    Vel = data[:,4].reshape((nx, ny)).transpose()
    return X, Y, T, D2, Vel, nx

def plot_graphics(t, name, xmin, xmax, ymin, ymax, segs, T, D2, Vel):
    fig, axs = plt.subplots(1, 3, figsize=(19.20, 10.80))
    
    # Common attributes for all subplots
    for ax in axs:
        rect = matplotlib.patches.Rectangle((xmin, ymin), xmax-xmin, ymax-ymin, linewidth=2, edgecolor='k', facecolor='none')
        ax.add_patch(rect)
        line_segments = LineCollection(segs, linewidths=4, colors='green', linestyle='solid')
        ax.add_collection(line_segments)
        ax.set_aspect('equal')
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.axis('off')

    # Individual subplots
    im0 = axs[0].imshow(T, cmap="coolwarm", interpolation='Bilinear', origin='lower', extent=[xmin, xmax, ymin, ymax], vmax=10.0, vmin=0.0)
    axs[0].set_title(r"$\phi, t = %5.4f$" % t, fontsize=20)
    fig.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04)

    im1 = axs[1].imshow(D2, cmap="hot", interpolation='Bilinear', origin='lower', extent=[xmin, xmax, ymin, ymax], vmax=1.0, vmin=-3.0)
    axs[1].set_title(r"$\|\mathcal{D}_{ij}\|, t = %5.4f$" % t, fontsize=20)
    fig.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.04)

    im2 = axs[2].imshow(Vel, cmap="Blues", interpolation='Bilinear', origin='lower', extent=[xmin, xmax, ymin, ymax], vmax=10.0, vmin=0.0)
    axs[2].set_title(r"$\|V_i\|, t = %5.4f$" % t, fontsize=20)  
    fig.colorbar(im2, ax=axs[2], fraction=0.046, pad=0.04) 

    plt.savefig(name, bbox_inches='tight', dpi=300)
    plt.close()

def process_file(ti, params):
    """Process a single timestep with all parameters passed as a dictionary"""
    t = params['tSnap'] * ti
    place = f"{params['caseToProcess']}/intermediate/snapshot-{t:5.4f}"
    name = f"{params['folder']}/{int(t*1000):08d}.png"

    if not os.path.exists(place):
        print(f"{place} File not found!")
        return None
    elif os.path.exists(name):
        print(f"{name} Image present!")
        return None

    segs = get_segs(place)
    X, Y, T, D2, Vel, nz = get_field_values(
        place, 
        params['xmin'], 
        params['xmax'], 
        params['ymin'], 
        params['ymax'], 
        params['ny']
    )
    
    plot_graphics(
        t, name, 
        params['xmin'], params['xmax'], 
        params['ymin'], params['ymax'], 
        segs, T, D2, Vel
    )

    print(f"Processed timestep {ti}")
    return nz

def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description='Make videos of standing waves')
    parser.add_argument('--num_workers', type=int, default=4, help='Number of workers')
    parser.add_argument('--tSnap', type=float, default=0.1, help='Snapshot time interval')
    parser.add_argument('--L0', type=float, default=10.0, help='Length of the domain')
    parser.add_argument('--caseToProcess', type=str, default='../', help='Case to process')
    parser.add_argument('--folderToSave', type=str, default='Video', help='Folder to save the video')

    args = parser.parse_args()

    nGFS = 500
    num_workers = min(args.num_workers, mp.cpu_count())
    folder = args.folderToSave
    os.makedirs(folder, exist_ok=True)

    # Parameters dictionary
    params = {
        'ny': 128,
        'xmin': -args.L0/2.,
        'xmax': args.L0/2.,
        'ymin': -args.L0/2.,
        'ymax': args.L0/2.,
        'lw': 2,
        'tSnap': args.tSnap,
        'folder': folder,
        'caseToProcess': args.caseToProcess
    }

    # Create a partial function with the parameters
    process_func = partial(process_file, params=params)

    # Create a pool of workers and process the files
    with mp.Pool(processes=num_workers) as pool:
        results = pool.map(process_func, range(nGFS))
        
        # Filter out None results and print completion
        completed = [r for r in results if r is not None]
        print(f"Completed {len(completed)} out of {nGFS} files")
        if completed:
            print(f"Last known nx value: {completed[-1]}")

if __name__ == "__main__":
    main()
