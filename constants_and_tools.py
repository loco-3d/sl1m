import numpy as np
from sl1m.tools.obj_to_constraints import load_obj, as_inequalities, rotate_inequalities, inequalities_to_Inequalities_object

from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm
import numpy as np

from scipy.spatial import ConvexHull
from qp import solve_lp

#~ import eigenpy
#from curves import bezier3
from random import random as rd
from random import randint as rdi
from numpy import squeeze, asarray


Id = array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
g = array([0.,0.,-9.81])
g6 = array([0.,0.,-9.81,0.,0.,0.])

mu = 0.45

x = array([1.,0.,0.])
z = array([0.,0.,1.])
zero3 = zeros(3) 


eps =0.000001

#### surface to inequalities ####    
def convert_surface_to_inequality(s):
    #TODO does normal orientation matter ?
    #it will for collisions
    n = cross(s[:,1] - s[:,0], s[:,2] - s[:,0])
    if n[2] <= 0.:
        for i in range(3):
            n[i] = -n[i]
    n /= norm(n)
    return surfacePointsToIneq(s, n)
    
def replace_surfaces_with_ineq_in_phaseData(phase):
    phase["S"] = [convert_surface_to_inequality(S) for S in phase["S"]]
    
def replace_surfaces_with_ineq_in_problem(pb):
    [ replace_surfaces_with_ineq_in_phaseData(phase) for phase in pb ["phaseData"]]

def ineqQHull(hull):
    A = hull.equations[:,:-1]
    b = -hull.equations[:,-1]
    return A,b
    
def vectorProjection (v, n):    
    v = v / norm(v)
    n = n / norm(n)
    proj = v - np.dot(v,n)*n
    return proj/norm(proj)
    
def default_transform_from_pos_normal_(transform, pos, normal):
    #return vstack( [hstack([transform,pos.reshape((-1,1))]), [ 0.        ,  0.        ,  0.        ,  1.        ] ] ) # FIXME : temp stuff, only work on flat floor
    # FIXME : there is something wrong the the code above
    #~ print "pos ", pos
    #~ print "normal ", normal
    #align the x-axis of the foot to the root heading direction    
    f = x
    xp = np.dot(transform, x).tolist()  
    t = vectorProjection(xp, normal)
    v = np.cross(f, t)
    c = np.dot(f, t)
    
    if abs(c) > 0.99 :
        return vstack( [hstack([transform,pos.reshape((-1,1))]), [ 0.        ,  0.        ,  0.        ,  1.        ] ] )     
    else:
        u = v/norm(v)
        h = (1. - c)/(1. - c**2)

        vx, vy, vz = v
        rot1 =array([[c + h*vx**2, h*vx*vy - vz, h*vx*vz + vy],
              [h*vx*vy+vz, c+h*vy**2, h*vy*vz-vx],
              [h*vx*vz - vy, h*vy*vz + vx, c+h*vz**2]])
    #align the z-axis of the foot to the surface normal
    f = z
    t = array(normal)
    t = t / norm(t)
    v = np.cross(f, t)
    c = np.dot(f, t)    
    if abs(c) > 0.99 :
        rot2 = identity(3)    
    else:
        u = v/norm(v)
        h = (1. - c)/(1. - c**2)

        vx, vy, vz = v
        rot2 =array([[c + h*vx**2, h*vx*vy - vz, h*vx*vz + vy],
              [h*vx*vy+vz, c+h*vy**2, h*vy*vz-vx],
              [h*vx*vz - vy, h*vy*vz + vx, c+h*vz**2]])
              
    rot = np.dot(rot1,rot2)
    return vstack( [hstack([rot,pos.reshape((-1,1))]), [ 0.        ,  0.        ,  0.        ,  1.        ] ] )
    

def default_transform_from_pos_normal(pos, normal):
    f = array([0.,0.,1.])
    t = array(normal)
    t = t / norm(t)
    v = np.cross(f, t)
    c = np.dot(f, t)    
    if abs(c) > 0.99 :
        rot = identity(3)    
    else:
        u = v/norm(v)
        h = (1. - c)/(1. - c**2)

        vx, vy, vz = v
        rot =array([[c + h*vx**2, h*vx*vy - vz, h*vx*vz + vy],
              [h*vx*vy+vz, c+h*vy**2, h*vy*vz-vx],
              [h*vx*vz - vy, h*vy*vz + vx, c+h*vz**2]])
    return vstack( [hstack([rot,pos.reshape((-1,1))]), [ 0.        ,  0.        ,  0.        ,  1.        ] ] )

    
#last is equality
def surfacePointsToIneq(S, normal):
    n = array(normal)
    tr = default_transform_from_pos_normal(array([0.,0.,0.]),n)
    trinv = tr.copy(); trinv[:3,:3] = tr[:3,:3].T;
        
    trpts = [tr[:3,:3].dot(s)[:2] for s in S.T]
    hull = ConvexHull(array(trpts))
    A,b = ineqQHull(hull)
    A = hstack([A,zeros((A.shape[0],1))])
    ine = inequalities_to_Inequalities_object(A,b)
    ine = rotate_inequalities(ine, trinv)
    #adding plane constraint
    #plane equation given by first point and normal for instance
    d = array([n.dot(S[:,0])])
    A = vstack([ine.A, n , -n ])
    b = concatenate([ine.b, d, -d]).reshape((-1,))
    A = vstack([ine.A, n])
    b = concatenate([ine.b, d]).reshape((-1,))
    return A, b
    
############ BENCHMARKING ###############

from time import  clock
def timMs(t1, t2):
    return (t2-t1) * 1000.
