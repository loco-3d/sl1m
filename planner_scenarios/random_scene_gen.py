from numpy import arange, array, zeros
from random import *


import matplotlib.pyplot as plt    
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import scipy as sp

NUM_SURF = randint(5,10)
RAND_W = False # default 0.5
RAND_H = True # default 0.15
RAND_HEIGHT = False # default 0.5
LARGE_SURF = False 
# RAND_X_ORI
# RAND_Y_ORI
# RAND_Z_ORI


def draw_rectangle(l, ax):
    #~ plotPoints(ax,l)
    lr = l + [l[0]]
    cx = [c[0] for c in lr]
    cy = [c[1] for c in lr]
    cz = [c[2] for c in lr]
    ax.plot(cx, cy, cz)
    
def draw_scene(ax = None, color = "p"):
    if ax is None:        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    [draw_rectangle(l,ax) for l in all_surfaces]
    return ax

def random_rectangle_gen (p, w, h):
  p1 = (array(p) + array([-h, w, 0])).tolist()
  p2 = (array(p) + array([h, w, 0])).tolist()
  p3 = (array(p) + array([h, -w, 0])).tolist()
  p4 = (array(p) + array([-h, -w, 0])).tolist()
  return [p1, p2, p3, p4]

def random_scene_gen ():
  surfaces = []
  p = zeros(3)
  
  w = round(uniform(0.5, 1.0),2) if RAND_W else 0.5
  h = round(uniform(0.15, 0.3),2) if RAND_H else 0.15
  height = sample([0., 0.5, 1.0, 1.5], 1) if RAND_HEIGHT else 0.
  
  for i in range(0, NUM_SURF):
    prev_w = w
    prev_h = h
  
    surfaces += [random_rectangle_gen(p, w, h)]
    
    w = round(uniform(0.5, 1.0),2) if RAND_W else 0.5
    h = round(uniform(0.15, 0.3),2) if RAND_H else 0.15
    height = sample([0., 0.5, 1.0, 1.5], 1) if RAND_HEIGHT else 0.
    p = (array(p) + array([-prev_h-h - 0.03, 0., height])).tolist()
  
  if not LARGE_SURF:
    p = [0, w*2 + 0.03, 0]
    for i in range(0, NUM_SURF):
      prev_w = w
      prev_h = h
      
      surfaces += [random_rectangle_gen(p, w, h)]
         
      w = round(uniform(0.5, 1.0),2) if RAND_W else 0.5
      h = round(uniform(0.15, 0.3),2) if RAND_H else 0.15
      height = sample([0., 0.5, 1.0, 1.5], 1) if RAND_HEIGHT else 0.
      p = (array(p) + array([-prev_h-h - 0.03, 0., height])).tolist() 
      
    
  return surfaces
    
  



############# main ###################    

if __name__ == '__main__':
  all_surfaces = random_scene_gen()
  ax = draw_scene()
  # ax.set_xlabel('X axis')
  # ax.set_ylabel('Y axis')
  # ax.set_zlabel('Z axis')
  plt.show()
  print "hello world"



"""
floor = [[-0.30, 0.54  , 0.  ], [0.01 ,  0.54, 0. ], [0.01 , -0.46, 0.  ], [-0.30, -0.46, 0.  ], ]
floor2  = [[-3., 0.4  , 0.  ], [-2.7 ,  0.4, 0. ], [-2.7 , 0.1, 0.  ], [-03., 0.1, 0.  ], ]
rub2  = [[ -2.11, 0.30 , 0.05 ], [-2.45 , 0.30, 0.05 ],  [ -2.45, 0.53, 0.05 ], [-2.11, 0.53, 0.05  ], ]
rub3  = [[ -1.91, -0.15 , 0.1 ], [-2.25 , -0.15, 0.1 ],  [ -2.25, 0.18, 0.1 ], [-1.91, 0.18, 0.1  ], ]
rub4  = [[ -1.69, 0.19 , 0.15 ], [-2.03 , 0.19, 0.15 ],  [ -2.03, 0.53, 0.15 ], [-1.69, 0.53, 0.15  ], ]
rub5  = [[ -1.49, -0.15 , 0.2 ], [-1.83 , -0.15, 0.2 ],  [ -1.83, 0.18, 0.2 ], [-1.49, 0.18, 0.2  ], ]
rub6  = [[ -1.29, 0.19 , 0.2 ], [-1.63 , 0.19, 0.2 ],  [ -1.63, 0.53, 0.2 ], [-1.29, 0.53, 0.2  ], ]
rub7  = [[ -1.09, -0.15 , 0.15 ], [-1.43 , -0.15, 0.15],  [ -1.43, 0.18, 0.15], [-1.09, 0.18, 0.15 ], ]
rub75  = [[ -0.89, 0.19 , 0.1 ], [-1.23 , 0.19, 0.1],  [ -1.23, 0.53, 0.1], [-0.89, 0.53, 0.1 ], ]
rub8  = [[ -0.89, -0.15 , 0.025 ], [-1.02 , -0.15, 0.025],  [ -1.02, 0.18, 0.025], [-0.89, 0.18, 0.025 ], ]
# rub9  = [[ -0.35, -0.15 , 0.025 ], [-0.86 , -0.15, 0.025], [-0.86, 0.52, 0.025 ], [ -0.35, 0.52, 0.025], ]
rub8  = [[ -0.89, -0.15 , 0.05 ], [-1.02 , -0.15, 0.05],  [ -1.02, 0.18, 0.05], [-0.89, 0.18, 0.05 ], ]
rub9  = [[ -0.35, -0.15 , 0.05 ], [-0.86 , -0.15, 0.05], [-0.86, 0.52, 0.05 ], [ -0.35, 0.52, 0.05], ]

all_surfaces = [floor, rub2, rub3, rub4, rub5, rub6, rub7, rub75, rub8, rub9, floor2]
# all_surfaces = [floor]

"""
