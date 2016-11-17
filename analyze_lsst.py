
# coding: utf-8

# In[1]:

# imports
from __future__ import division
import os, sys, time, math
import optparse
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (13, 8) if False else (10, 6)

import pyzdde.zdde as pyz
import pyzdde.arraytrace as at

def parse_commandline():
    """
    Parse the options given on the command-line.
    """
    parser = optparse.OptionParser()

    parser.add_option("-t","--type",default="lsst_nominal")

    opts, args = parser.parse_args()

    return opts

def many_spots(ln, hxs,hys,wavenum=1,spirals=10,rays=600):
    fovxs, fovys, xs, ys = np.array([]), np.array([]), np.array([]), np.array([])
    for hx, hy in zip(hxs,hys):
        (x,y,z,intensity) = ln.zSpiralSpot(hx,hy,wavenum,spirals,rays)
        fovxs, fovys, xs, ys = np.append(fovxs,hx*np.ones((len(x),1))), np.append(fovys,hy*np.ones((len(y),1))), np.append(xs,x), np.append(ys,y)
    return fovxs, fovys, xs, ys

# Spot diagram analysis functions
def zGridSpot(ln, hx, hy, waveNum, dx, dy, mode=0):

    pxs = np.arange(-1.0,1.0,dx)
    pys = np.arange(-1.0,1.0,dy)
    [PXs,PYs] = np.meshgrid(pxs,pys)
    pxs, pys = PXs.flatten(), PYs.flatten()
    
    x = [] # x-coordinate of the image surface
    y = [] # y-coordinate of the image surface
    z = [] # z-coordinate of the image surface
    intensity = [] # the relative transmitted intensity of the ray
    for px, py in zip(pxs,pys):
        rayTraceData = ln.zGetTrace(waveNum, mode, -1, hx, hy, px, py)
        if rayTraceData[0] == 0:
            x.append(rayTraceData[2])
            y.append(rayTraceData[3])
            z.append(rayTraceData[4])
            intensity.append(rayTraceData[11])
        else:
            print("Raytrace Error")
            #exit()
            # !!! FIX raise an error here
    return (x, y, z, intensity)

def many_grid_spots(ln, hxs, hys, dx, dy, wavenum=1):
    fovxs, fovys, xs, ys = np.array([]), np.array([]), np.array([]), np.array([])
    for hx, hy in zip(hxs,hys):
        (x,y,z,intensity) = zGridSpot(ln,hx,hy,wavenum,dx,dy,mode=0)
        fovxs, fovys, xs, ys = np.append(fovxs,hx*np.ones((len(x),1))), np.append(fovys,hy*np.ones((len(y),1))), np.append(xs,x), np.append(ys,y)
    return fovxs, fovys, xs, ys

# Parse command line
opts = parse_commandline()

plotDir = os.path.join('C:\Users\CAD User\Desktop\lsstzemax\data',opts.type)
if not os.path.isdir(plotDir): os.mkdir(plotDir)

ln = pyz.createLink() # create a DDE link object for communication

zfile = os.path.join('C:\Users\CAD User\Desktop\lsstzemax\zemax', '%s.zmx'%opts.type)

ln.zLoadFile(zfile)
# Surfaces in the sequential lens data editor
ln.ipzGetLDE()
# General System properties
ln.zGetSystem()
# Paraxial/ first order properties of the system
ln.zGetFirst()
# duplicate of zGetFirst() for use in the notebook
ln.ipzGetFirst()
# ... another example is the zGetSystemAper() that returns information about the aperture. 
# The aperture type is retuned as a code which we might not remember always ...
ln.zGetSystemAper()
# ...with the duplicate, ipzGetSystemAper(), we can immediately know that
# the aperture type is the Entrance Pupil Diameter (EPD)
ln.ipzGetSystemAper()
# information about the field definition
ln.ipzGetFieldData()

# In[13]:

hx = 0.0
hy = 0.0
spirals = 10 #100
rays = 600   #6000
(xu,yu,zu,intensityu) = ln.zSpiralSpot(hx,hy,1,spirals,rays)
(xg,yg,zg,intensityg) = ln.zSpiralSpot(hx,hy,2,spirals,rays)
(xr,yr,zr,intensityr) = ln.zSpiralSpot(hx,hy,3,spirals,rays)
(xi,yi,zi,intensityi) = ln.zSpiralSpot(hx,hy,4,spirals,rays)
(xz,yz,zz,intensityz) = ln.zSpiralSpot(hx,hy,5,spirals,rays)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111)
ax.set_aspect('equal')
ax.scatter(xu,yu,s=5,c='blue',linewidth=0.35,zorder=20)
ax.scatter(xg,yg,s=5,c='green',linewidth=0.35,zorder=21)
ax.scatter(xr,yr,s=5,c='red',linewidth=0.35,zorder=22)
ax.scatter(xi,yi,s=5,c='magenta',linewidth=0.35,zorder=23)
ax.scatter(xz,yz,s=5,c='yellow',linewidth=0.35,zorder=24)
ax.set_xlabel('x');ax.set_ylabel('y')
fig.suptitle('Spiral Spot')
ax.grid(color='lightgray', linestyle='-', linewidth=1)
ax.ticklabel_format(scilimits=(-2,2))

plt.show()
plotName = os.path.join(plotDir,'spiral.png')
plt.savefig(plotName)
plotName = os.path.join(plotDir,'spiral.eps')
plt.savefig(plotName)
plotName = os.path.join(plotDir,'spiral.pdf')
plt.savefig(plotName)
plt.close()

hxs, hys = np.linspace(-1.75,1.75,10), np.linspace(-1.75,1.75,10)
[HXs,HYs] = np.meshgrid(hxs,hys)
hxs, hys = HXs.flatten(), HYs.flatten()
dx = 0.1
dy = 0.1
fovx, fovy, x, y = many_grid_spots(ln, hxs, hys, dx, dy, wavenum = 1)
fov = np.sqrt(fovx**2 + fovy**2)

fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111)
ax.set_aspect('equal')
plt.scatter(x,y,s=5,c=fov,linewidth=0.35,zorder=21)
plt.xlabel('x [mm]')
plt.ylabel('y [mm]')
cbar = plt.colorbar()
cbar.set_label('FOV [degrees]')
plt.grid(color='lightgray', linestyle='-', linewidth=1)
#ax.ticklabel_format(scilimits=(-2,2))
plt.axis([-600,600,-600,600])

plt.show()

plotName = os.path.join(plotDir,'many_spots.png')
plt.savefig(plotName)
plotName = os.path.join(plotDir,'many_spots.eps')
plt.savefig(plotName)
plotName = os.path.join(plotDir,'many_spots.pdf')
plt.savefig(plotName)
plt.close()

filename = os.path.join(plotDir,'spots.dat')
fid = open(filename,'w')
for a,b,c,d in zip(fovx,fovy,x,y):
    fid.write('%.5f %.5f %.5f %.5f\n'%(a,b,c,d))
fid.close()

