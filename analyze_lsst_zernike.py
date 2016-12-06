
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

    parser.add_option("-t","--type",default="lsst_zernike")
    parser.add_option("-o","--output",default="noabberation")
    parser.add_option("-a","--astigmatism",default=0.0,type=float)
    parser.add_option("-d","--defocus",default=0.0,type=float)
    parser.add_option("-c","--coma",default=0.0,type=float)
    parser.add_option("-s","--spherical",default=0.0,type=float)
    parser.add_option("--doPlots", action ="store_true", default=False)

    opts, args = parser.parse_args()

    return opts

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
            x.append(np.nan)
            y.append(np.nan)
            z.append(np.nan)
            intensity.append(np.nan)
            #exit()
            # !!! FIX raise an error here
    return (x, y, z, intensity)

def many_grid_spots(ln, hxs, hys, dx, dy, wavenum=1):
    fovxs, fovys, xs, ys = np.array([]), np.array([]), np.array([]), np.array([])
    for hx, hy in zip(hxs,hys):
        (x,y,z,intensity) = zGridSpot(ln,hx,hy,wavenum,dx,dy,mode=0)
        fovxs, fovys, xs, ys = np.append(fovxs,hx*np.ones((len(x),1))), np.append(fovys,hy*np.ones((len(y),1))), np.append(xs,x), np.append(ys,y)
    return fovxs, fovys, xs, ys

def get_spots(ln):
    hxs, hys = np.linspace(-1.75,1.75,10), np.linspace(-1.75,1.75,10)
    [HXs,HYs] = np.meshgrid(hxs,hys)
    hxs, hys = HXs.flatten(), HYs.flatten()
    dx = 0.1 
    dy = 0.1 
    fovx, fovy, x, y = many_grid_spots(ln, hxs, hys, dx, dy, wavenum = 1)
    fov = np.sqrt(fovx**2 + fovy**2)

    return fov, fovx, fovy, x, y

# Parse command line
opts = parse_commandline()

plotDir = os.path.join('C:\Users\CAD User\Desktop\lsstzemax\data',opts.output)
if not os.path.isdir(plotDir): os.mkdir(plotDir)

ln = pyz.createLink() # create a DDE link object for communication

zfile = os.path.join('C:\Users\CAD User\Desktop\lsstzemax\zemax', '%s.zmx'%opts.type)

ln.zLoadFile(zfile)

fov_original, fovx_original, fovy_original, x_original, y_original = get_spots(ln)

# Set astigmatism, coma. defocus
ln.zSetSurfaceParameter(28,18,opts.defocus)
ln.zSetSurfaceParameter(28,19,opts.astigmatism)
ln.zSetSurfaceParameter(28,20,opts.astigmatism)
ln.zSetSurfaceParameter(28,21,opts.coma)
ln.zSetSurfaceParameter(28,22,opts.coma)
ln.zSetSurfaceParameter(28,25,opts.spherical)

fov, fovx, fovy, x, y = get_spots(ln)

if opts.doPlots:
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
for a,b,c,d,e,f in zip(fovx,fovy,x,y,x_original,y_original):
    fid.write('%.5f %.5f %.5f %.5f %.5f %.5f\n'%(a,b,c,d,e,f))
fid.close()

fovx_unique = np.unique(fovx)
fovy_unique = np.unique(fovy)

filename = os.path.join(plotDir,'spots_mean.dat')
fid = open(filename,'w')
for fox in fovx_unique:
    for foy in fovy_unique:
        idx = np.intersect1d(np.where(fovx==fox)[0],np.where(fovy==foy)[0])
        x_mean = np.nanmean(x[idx])
        y_mean = np.nanmean(y[idx])
        x_original_mean = np.nanmean(x_original[idx])
        y_original_mean = np.nanmean(y_original[idx])
        fid.write("%.5f %.5f %.5f %.5f %.5f %.5f\n"%(fox,foy,x_mean,y_mean,x_original_mean,y_original_mean))
fid.close()


