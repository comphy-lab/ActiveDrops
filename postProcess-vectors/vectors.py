import numpy as np
import os
import subprocess as sp
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import matplotlib.gridspec as gridspec
import sys

matplotlib.rcParams['text.usetex'] = True

def execute_process(exe):
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    return stderr.decode("utf-8").split("\n")

def gettingFacets(filename, tracer):

    print('Getting facets values')
    if tracer == 1:
        exe = ["./getFacets", filename]
    else:
        exe = ["./getFacet2", filename]

    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    segs = []
    skip = False
    if (len(temp2) > 1e1):
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == ['']:
                skip = False
                pass
            else:
                if not skip:
                    temp4 = temp2[n1+1].split(" ")
                    r1, z1 = np.array([float(temp3[1]), float(temp3[0])])
                    r2, z2 = np.array([float(temp4[1]), float(temp4[0])])
                    
                    segs.append(((z1, r1),(z2,r2)))
                    skip = True
    print('Got facets values for tracer')

    return segs

def gettingVcm(filename):
    exe = ["./getVelocity_v2", filename]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split(" ")
    return temp2

def gettingCM(filename):
    exe = ["./getCM", filename]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split(" ")
    return temp2

def gettingfield(filename):
    print('Field values')
    exe = ["./getDataSlice", filename, str(xmin), str(xmax), str(ymin), str(ymax), str(nx), str(ny)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    Xtemp, Ytemp, ftemp, Utemp, Vtemp = [], [], [], [], []

    if (len(temp2) > 1e2):
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == ['']:
                pass
            else:
                Xtemp.append(float(temp3[0]))
                Ytemp.append(float(temp3[1]))
                ftemp.append(float(temp3[2]))
                Utemp.append(float(temp3[3]))
                Vtemp.append(float(temp3[4]))

        X = np.asarray(Xtemp)
        Y = np.asarray(Ytemp)
        f = np.asarray(ftemp)
        U = np.asarray(Utemp)
        V = np.asarray(Vtemp)

        X.resize((ny+1, nx+1))
        Y.resize((ny+1, nx+1))
        f.resize((ny+1, nx+1))
        U.resize((ny+1, nx+1))
        V.resize((ny+1, nx+1))

        print('Got Field values')
        return X, Y, f, U, V
    else:
        return Xtemp, Ytemp, ftemp, Utemp, Vtemp
    
def get_field_values(place, xmin, xmax, ymin, ymax, ny):
    temp2 = list(filter(None, execute_process(["./getData", place, str(xmin), str(ymin), str(xmax), str(ymax), str(ny)])))
    data = np.array([line.split() for line in temp2], dtype=float)
    nx = data.shape[0] // ny
    X = data[:,0].reshape((nx, ny)).transpose()
    Y = data[:,1].reshape((nx, ny)).transpose()
    T = data[:,2].reshape((nx, ny)).transpose()
    D2 = data[:,3].reshape((nx, ny)).transpose()
    Vel = data[:,4].reshape((nx, ny)).transpose()
    return T

# ------------------------------------------------------------------------------

nGFS = 270
lw = 4

folder = 'VelVectors'  # output folder
if not os.path.isdir(folder):
    os.makedirs(folder)

ux = np.array([])
uy = np.array([])
tt = np.array([])

LEVEL = 8
nx = 2**(LEVEL)
ny = 2**(LEVEL)

for ti in range(nGFS):
    t = 0.1 * ti
    place = "../intermediate/snapshot-%5.4f" % t
    name = "%s/%8.8d.png" %(folder, int(t*1e3))
    if not os.path.exists(place):
        print("File %s not found!" % place)
    else:
        if os.path.exists(name):
            print("Image %s found!" % name)
        else:
            CM = gettingCM(place)

            CMx = float(CM[0])
            CMy = float(CM[1])

            xmin, xmax, ymin, ymax = [-5.0 + CMx, 5.0 + CMx, -5 + CMy, 5 + CMy]

            X, Y, f, U, V = gettingfield(place)
            facets = gettingFacets(place, 1)
            vcm = gettingVcm(place)
            T = get_field_values(place,xmin, xmax, ymin, ymax, ny)

            if (len(facets)):

                gs = gridspec.GridSpec(15,20)
                
                fig = plt.figure(figsize=(30,20))
                ax = plt.subplot(gs[:, 0:12])
                ax2 = plt.subplot(gs[4:11, 13:19])

                fig.set_size_inches(25.20, 12.80)

                rmin, rmax, zmin, zmax = [-5.0 + CMx, 5.0 + CMx, -5 + CMy, 5 + CMy]
                
                ax.plot([rmin, rmin], [zmin, zmax],'-',color='black',linewidth=lw/2)
                ax.plot([rmin, rmax], [zmin, zmin],'-',color='black',linewidth=lw/2)
                ax.plot([rmin, rmax], [zmax, zmax],'-',color='black',linewidth=lw/2)
                ax.plot([rmax, rmax], [zmin, zmax],'-',color='black',linewidth=lw/2)

                line_segments1 = LineCollection(facets, linewidths=3, colors='red', linestyle='solid')
                ax.add_collection(line_segments1)
                speed = np.sqrt(U**2 + V**2)
                ndx = 10
                ndy = 10
                vScale = 12.5

                ax.quiver(X[::ndx,::ndy], Y[::ndx,::ndy], U[::ndx,::ndy], V[::ndx,::ndy], scale=vScale, color='black',linewidth=0.5)
                ax.imshow(T, cmap="Blues", interpolation='Bilinear', origin='lower', extent=[xmin, xmax, ymin, ymax], vmax=10.0, vmin=0.0)
                ax.set_aspect('equal')
                ax.set_xlim(rmin, rmax)
                ax.set_ylim(zmin, zmax)

                ax.set_title(r'$t/t_{\gamma} = %3.2f$' % t, fontsize=30)

                ax.axis('off')

                ux = np.append(ux, -float(vcm[1]))
                uy = np.append(uy, float(vcm[0]))
                tt = np.append(tt, float(vcm[2]))

                Ucm = np.sqrt(ux**2 + uy**2)

                plt.rc('text', usetex=True)
                plt.rc('font', family='serif')

                ax2.plot(tt, Ucm, linewidth=5)

                ax2.set_xlabel("$t$", fontsize=28)
                ax2.set_ylabel("$V$", fontsize=28)
                ax2.set_xlim([0, 27])
                ax2.set_ylim([0, 2])

                ax2.tick_params(axis='both', which='major', labelsize=20)
                ax2.tick_params(axis='both', which='minor', labelsize=20)

                plt.savefig(name, bbox_inches="tight")
                plt.close()
            else:
                print("Problem in the available file %s" % place)

    print(("Done %d of %d" % (ti+1, nGFS)))
