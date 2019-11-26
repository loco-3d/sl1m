from sl1m.problem_definition import *

from pickle import load
import pickle

import sys

def readFromFile (fileName):
  data = []
  #~ print "filename ", fileName
  with open(fileName,'rb') as f:
    #~ print "ok ", f
    while True:
      try:
        data.append(pickle.load(f))
      except EOFError:
        break
  return data[0]


import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import scipy as sp

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def draw_rectangle(l, ax):
    #~ plotPoints(ax,l)
    lr = l + [l[0]]
    cx = [c[0] for c in lr]
    cy = [c[1] for c in lr]
    cz = [c[2] for c in lr]
    ax.plot(cx, cy, cz)
    
all_surfaces = None
    
def draw_scene(surfaces, ax = None, color = "p"):
    if ax is None:        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    [draw_rectangle(l,ax) for l in all_surfaces]
    return ax


if __name__ == '__main__':
    
    
    from sl1m.fix_sparsity import solveL1, solveMIP
    
    
    fileName = sys.argv[1]    
    pb = readFromFile(fileName+'_problem')    
    surfaces = readFromFile(fileName+'_surf')    
    all_surfaces = readFromFile(fileName+'_scene')  
    pb['p0'] = None
    solveMIP(pb, surfaces, MIP = True, draw_scene = draw_scene, plot = True, convert = False)
    
    pb = readFromFile(fileName+'_problem') 
    pb['p0'] = None
    solveL1(pb, surfaces, draw_scene, plot = True, convert = False)
