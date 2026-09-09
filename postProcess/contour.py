"""
# contour.py

Render per-snapshot panels of the species concentration, the deformation-rate
norm and the speed for one ActiveDrops case, in parallel over snapshots.

Run from `postProcess/` with a case directory that contains
`intermediate/snapshot-<t>` files:

    python3 contour.py --caseToProcess ../simulationCases/c1000 --cpus 4

The `getFacets` and `getData` helpers are compiled once with `qcc` before any
snapshot is processed; workers share the read-only executables. Frames are
named by the snapshot time so serial and parallel runs produce the same set.

## Author
Vatsal Sanjay
Email: vatsal.sanjay@comphy-lab.org
Computational Multiphase Physics (CoMPhy) Lab, Durham University
Last updated: Sep 9, 2026
"""

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
import shutil


matplotlib.rcParams['font.family'] = 'serif'
# Use LaTeX text rendering only when a latex executable is available; otherwise mathtext.
matplotlib.rcParams['text.usetex'] = shutil.which('latex') is not None

def execute_process(exe):
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    return stderr.decode("utf-8").split("\n")


def compile_helpers(helpers=("getFacets", "getData")):
    """Compile the Basilisk snapshot readers once, before any snapshot is processed."""
    if not shutil.which("qcc"):
        sys.exit("qcc not found; source the repository .project_config first.")
    # qcc resolves its intermediate files relative to the working directory, so
    # compile with relative names from the postProcess directory.
    here = os.path.dirname(os.path.abspath(__file__))
    for name in helpers:
        src = os.path.join(here, f"{name}.c")
        exe = os.path.join(here, name)
        if os.path.exists(exe) and os.path.getmtime(exe) >= os.path.getmtime(src):
            continue
        cmd = ["qcc", "-O2", "-Wall", "-disable-dimensions", f"{name}.c", "-o", name, "-lm"]
        print("Compiling:", " ".join(cmd), flush=True)
        sp.check_call(cmd, cwd=here)

def get_segs(place):
    temp2 = execute_process(["./getFacets", place])
    segs = []
    temp2 = list(filter(None, temp2))
    # getFacets prints segment endpoints in pairs of lines; ignore a trailing odd line.
    for n1 in range(0, len(temp2) - 1, 2):
        temp3 = temp2[n1].split()
        temp4 = temp2[n1+1].split()
        if len(temp3) < 2 or len(temp4) < 2:
            continue
        x1, y1 = map(float, temp3[:2])
        x2, y2 = map(float, temp4[:2])
        segs.append(((x1, y1), (x2, y2)))
    return segs

def get_field_values(place, xmin, xmax, ymin, ymax, ny):
    temp2 = list(filter(None, execute_process(["./getData", place, str(xmin), str(ymin), str(xmax), str(ymax), str(ny)])))
    rows = [line.split() for line in temp2]
    if not rows or any(len(r) != 5 for r in rows) or len(rows) % ny != 0:
        return None
    data = np.array(rows, dtype=float)
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
    fields = get_field_values(
        place,
        params['xmin'],
        params['xmax'],
        params['ymin'],
        params['ymax'],
        params['ny']
    )
    if fields is None:
        print(f"{place} incomplete field output; skipped")
        return None
    X, Y, T, D2, Vel, nz = fields

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
    parser = argparse.ArgumentParser(description='Render concentration, deformation-rate and speed panels per snapshot')
    parser.add_argument('--cpus', '--CPUs', '--num_workers', dest='cpus', type=int, default=4,
                        help='Number of parallel workers (default 4)')
    parser.add_argument('--tSnap', type=float, default=0.1, help='Snapshot time interval')
    parser.add_argument('--L0', type=float, default=10.0, help='Length of the domain')
    parser.add_argument('--caseToProcess', type=str, default='../simulationCases/c1000',
                        help='Case directory containing intermediate/ (default ../simulationCases/c1000)')
    parser.add_argument('--folderToSave', type=str, default='Video', help='Folder for the rendered frames')
    parser.add_argument('--max-frames', type=int, default=500,
                        help='Maximum number of snapshots to consider (default 500)')

    args = parser.parse_args()
    if args.cpus <= 0:
        parser.error('--cpus must be a positive integer')
    if args.max_frames <= 0:
        parser.error('--max-frames must be a positive integer')

    os.environ.setdefault('OMP_NUM_THREADS', '1')
    compile_helpers()
    here = os.path.dirname(os.path.abspath(__file__))
    os.chdir(here)

    nGFS = args.max_frames
    num_workers = min(args.cpus, mp.cpu_count())
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
