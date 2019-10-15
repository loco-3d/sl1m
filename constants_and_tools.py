import numpy as np
from hpp.corbaserver.rbprm.tools.obj_to_constraints import load_obj, as_inequalities, rotate_inequalities, inequalities_to_Inequalities_object
from hpp_centroidal_dynamics import *
from hpp_spline import *
from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, stack, vstack, hstack, asmatrix, identity, cross, concatenate
from numpy.linalg import norm
import numpy as np

from scipy.spatial import ConvexHull
from hpp_bezier_com_traj import *
from qp import solve_lp

import eigenpy
import cdd
#from curves import bezier3
from random import random as rd
from random import randint as rdi
from numpy import squeeze, asarray


Id = matrix([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
g = array([0.,0.,-9.81])
g6 = array([0.,0.,-9.81,0.,0.,0.])

mu = 0.45

x = array([1.,0.,0.])
z = array([0.,0.,1.])
zero3 = zeros(3) 


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


def addHeightConstraint(K,k, val):
    K1 = vstack([K, -z])
    k1 = concatenate([k, -ones(1) * val]).reshape((-1,))    
    return K1, k1

def skew(x):
    res = zeros((3,3));
    res[0][1] = - x[2]; res[0][2] =   x[1];
    res[1][0] =   x[2]; res[1][2] = - x[0];
    res[2][0] = - x[1]; res[2][1] =   x[0];
    return res;

gX = skew(g)

def vec3to6(x):
    res = zeros((6,))
    res[:3] = x
    return res

def w(m, c, ddc = zeros(3)):
    res = zeros(6)
    res[:3] = m * (ddc-g )
    res[3:] = cross(c,res[:3])
    return res

def generators(A,b, Aeq = None, beq = None):
    m = np.hstack([b,-A])
    matcdd = cdd.Matrix(m); matcdd.rep_type = cdd.RepType.INEQUALITY
    
    if Aeq is not None:
        meq = np.hstack([beq,-Aeq])
        matcdd.extend(meq.tolist(), True)
    
    H = cdd.Polyhedron(matcdd)
    g = H.get_generators()
    
    return [array(g[el][1:]) for el in range(g.row_size)], H
    
def filter(pts):
    hull = ConvexHull(pts, qhull_options='Q12')
    return [pts[i] for i in hull.vertices.tolist()]
    
def ineq(pts, canonicalize = False):
    apts = array(pts)
    m = np.hstack([ones((apts.shape[0],1)),apts])
    matcdd = cdd.Matrix(m); matcdd.rep_type = cdd.RepType.GENERATOR
    H = cdd.Polyhedron(matcdd)
    bmA = H.get_inequalities()
    if canonicalize:
        bmA.canonicalize()
    Ares = zeros((bmA.row_size,bmA.col_size-1))
    bres = zeros(bmA.row_size )
    for i in range(bmA.row_size):
        l = array(bmA[i])
        Ares[i,:] = -l[1:]
        bres[i]   =  l[0]
    return Ares, bres
    
def ineqQHull(hull):
    A = hull.equations[:,:-1]
    b = -hull.equations[:,-1]
    return A,b
    
    
def canon(A,b):
    m = np.hstack([b,-A])
    matcdd = cdd.Matrix(m); matcdd.rep_type = 1
    H = cdd.Polyhedron(matcdd)
    bmA = H.get_inequalities()
    #~ bmA.canonicalize()
    Ares = zeros((bmA.row_size,bmA.col_size-1))
    bres = zeros((bmA.row_size,1 ))
    for i in range(bmA.row_size):
        #~ print "line ", array(bmA[i])
        #~ print "A ", A[i][:]
        #~ print "b ", b[i]
        l = array(bmA[i])
        Ares[i,:] = -l[1:]
        bres[i]   =  l[0]
        #~ print "Ares ",Ares[i,:]
        #~ print "bres ",bres[i]
    return Ares, bres



def getCrocStateForward(c0, dc0, ddc0, x):
    c0a = array(c0.T.tolist()[0])
    dc0a = array(dc0.T.tolist()[0])
    ddc0a = array(ddc0.T.tolist()[0])
    #~ end_vel = (x - ddc0a) * 3. / 1.
    end_vel = (x - dc0a) * 3. / 1.
    end_acc = (-2.*dc0a  + x + ddc0a) * 3. * 2. / 1.
    return (matrix(x).T, matrix(end_vel).T, matrix(end_acc).T)
    #~ return (x, end_vel, end_acc)
    
def getCrocBackward(x, ddc1, dc1, c1, T):
    n = 3
    start_vel = n / T * (ddc1 - x)
    start_acc = n * (n-1) / (T*T) * (dc1 - 2* ddc1 + x)
    return (x,start_vel,start_acc)
    #~ return (x, end_vel, end_acc)

def genPolytope(A,b):
    pts, H = generators(A,b)
    apts = array(pts)
    if len(apts) > 0:
        hull = ConvexHull(apts)
        return hull, pts, apts, H
    return None, None, None, None
    
def convex_hull_ineq(pts):
    return None
    
    m = cData.contactPhase_.getMass()
    #~ #get 6D polytope
    (H, h) = ineqFromCdata(cData)
    #project to the space where aceleration is 0
    D = zeros((6,3))
    D[3:,:] = m * gX
    
    d = zeros((6,))
    d[:3] = -m * g
    
    A =     H.dot(D)
    b = h.reshape((-1,)) - H.dot(d)    
    #add kinematic polytope
    (K,k) = (cData.Kin_[0], cData.Kin_[1].reshape(-1,))
    
    resA = vstack([A, K])
    resb = concatenate([b, k]).reshape((-1,1))
    
    #DEBUG
    allpts = generators(resA,resb)[0]
    error = False
    for pt in allpts:
        print "pt ", pt
        assert (resA.dot(pt.reshape((-1,1))) - resb).max() <0.001, "antecedent point not in End polytope"  + str((resA.dot(pt.reshape((-1,1))) - resb).max())
        if (H.dot(w(m,pt).reshape((-1,1))) - h).max() > 0.001:
            error = True
            print "antecedent point not in End polytope"  + str((H.dot(w(m,pt).reshape((-1,1))) - h).max())
    assert not error, str (len(allpts))
    
    return (resA, resb)
    #~ return (A, b)
    #~ return (vstack([A, K]), None)

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
    #~ print "pos ", pos
    #~ print "normal ", normal
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

import os
if "INSTALL_HPP_DIR" in os.environ:
    insdir = os.environ['INSTALL_HPP_DIR'] + '/share/'
elif "DEVEL_HPP_DIR" in os.environ:
    insdir = os.environ['DEVEL_HPP_DIR'] + '/install/share/'

# load hrp2 constraints


__ineq_right_foot_hrp2 = None
__ineq_left_foot_hrp2  = None
__ineq_right_foot_hrp2_reduced = None
__ineq_left_foot_hrp2_reduced  = None

# add foot offset
def right_foot_hrp2_constraints(transform):
    global __ineq_right_foot_hrp2
    if __ineq_right_foot_hrp2 is None:
        #~ filekin = rob.rLegKinematicConstraints.replace("package://", insdir)
        filekin = insdir +"hrp2-rbprm/com_inequalities/feet_quasi_flat/hrp2_COM_constraints_in_RF_effector_frame_REDUCED.obj"
        #~ filekin = insdir +"hrp2-rbprm/com_inequalities/feet_quasi_flat/hrp2_COM_constraints_in_RF_effector_frame.obj"            
        #~ filekin = filekin.replace(".obj", "_reduced.obj")
        obj = load_obj(filekin)
        __ineq_right_foot_hrp2 = as_inequalities(obj)
    transform2 = transform.copy()
    transform2[2,3] += 0.105
    #~ transform2[2,3] += 0.105
    #~ transform2[2,3] += 0.105
    #~ transform2[2,3] += 0.01
    #~ transform2[2,3] += 0.105
    ine = rotate_inequalities(__ineq_right_foot_hrp2, transform2)
    return (ine.A, ine.b)
    
def right_foot_hrp2_constraints_reduced(transform):
    global __ineq_right_foot_hrp2_reduced
    if __ineq_right_foot_hrp2_reduced is None:
        #~ filekin = rob.rLegKinematicConstraints.replace("package://", insdir)
        filekin = insdir +"hrp2-rbprm/com_inequalities/feet_quasi_flat/hrp2_COM_constraints_in_RF_effector_frame_REDUCED.obj"
        obj = load_obj(filekin)
        right_foot_hrp2_constraints_reduced = as_inequalities(obj)
    transform2 = transform.copy()
    transform2[2,3] += 0.105
    #~ transform2[2,3] += 0.105
    #~ transform2[2,3] += 0.105
    #~ transform2[2,3] += 0.01
    #~ transform2[2,3] += 0.105
    ine = rotate_inequalities(__ineq_right_foot_hrp2_reduced, transform2)
    return (ine.A, ine.b)
        
def left_foot_hrp2_constraints(transform):
    global __ineq_left_foot_hrp2
    if __ineq_left_foot_hrp2 is None:
        #~ filekin = rob.lLegKinematicConstraints.replace("package://", insdir)
        filekin = insdir +"hrp2-rbprm/com_inequalities/feet_quasi_flat/hrp2_COM_constraints_in_LF_effector_frame_REDUCED.obj"
        #~ filekin = insdir +"hrp2-rbprm/com_inequalities/feet_quasi_flat/hrp2_COM_constraints_in_LF_effector_frame.obj"
        #~ filekin = filekin.replace(".obj", "_reduced.obj")
        obj = load_obj(filekin)
        __ineq_left_foot_hrp2 = as_inequalities(obj)
    transform2 = transform.copy()
    transform2[2,3] += 0.105
    #~ transform2[2,3] += 0.105
    #~ transform2[2,3] += 0.105
    #~ transform2[2,3] += 0.01
    #~ transform2[2,3] += 0.105
    ine = rotate_inequalities(__ineq_left_foot_hrp2, transform2)
    return (ine.A, ine.b)
    
def left_foot_hrp2_constraints_reduced(transform):
    global __ineq_left_foot_hrp2_reduced
    if __ineq_left_foot_hrp2_reduced is None:
        filekin = insdir +"hrp2-rbprm/com_inequalities/feet_quasi_flat/hrp2_COM_constraints_in_LF_effector_frame_REDUCED.obj"
        obj = load_obj(filekin)
        __ineq_left_foot_hrp2_reduced = as_inequalities(obj)
    transform2 = transform.copy()
    transform2[2,3] += 0.105
    #~ transform2[2,3] += 0.105
    #~ transform2[2,3] += 0.105
    #~ transform2[2,3] += 0.01
    #~ transform2[2,3] += 0.105
    ine = rotate_inequalities(__ineq_left_foot_hrp2_reduced, transform2)
    return (ine.A, ine.b)

def right_foot_talos_constraints(transform):
    global __ineq_right_foot_hrp2
    if __ineq_right_foot_hrp2 is None:
        filekin = insdir +"talos-rbprm/com_inequalities/feet_quasi_flat/talos_COM_constraints_in_RF_effector_frame_REDUCED.obj"
        obj = load_obj(filekin)
        __ineq_right_foot_hrp2 = as_inequalities(obj)
    ine = rotate_inequalities(__ineq_right_foot_hrp2, transform.copy())
    return (ine.A, ine.b)
    
def left_foot_talos_constraints(transform):
    global __ineq_left_foot_hrp2
    if __ineq_left_foot_hrp2 is None:
        filekin = insdir +"talos-rbprm/com_inequalities/feet_quasi_flat/talos_COM_constraints_in_LF_effector_frame_REDUCED.obj"
        obj = load_obj(filekin)
        __ineq_left_foot_hrp2 = as_inequalities(obj)
    ine = rotate_inequalities(__ineq_left_foot_hrp2, transform.copy())
    return (ine.A, ine.b)
    

__ineq_rf_in_rl_hrp2 = None
__ineq_lf_in_rf_hrp2  = None
    
# add foot offset
def right_foot_in_lf_frame_hrp2_constraints(transform):
    global __ineq_rf_in_rl_hrp2
    if __ineq_rf_in_rl_hrp2 is None:
        filekin = insdir +"hrp2-rbprm/relative_effector_positions/hrp2_LF_constraints_in_RF_quasi_flat_REDUCED.obj"
        #~ filekin = filekin.replace(".obj", "_reduced.obj")
        obj = load_obj(filekin)
        __ineq_rf_in_rl_hrp2 = as_inequalities(obj)
    transform2 = transform.copy()
    #~ transform2[2,3] += 0.105
    ine = rotate_inequalities(__ineq_rf_in_rl_hrp2, transform2)
    return (ine.A, ine.b)
        
def left_foot_in_rf_frame_hrp2_constraints(transform):
    global __ineq_lf_in_rf_hrp2
    if __ineq_lf_in_rf_hrp2 is None:
        filekin = insdir +"hrp2-rbprm/relative_effector_positions/hrp2_RF_constraints_in_LF_quasi_flat_REDUCED.obj"
        obj = load_obj(filekin)
        __ineq_lf_in_rf_hrp2 = as_inequalities(obj)
    transform2 = transform.copy()
    #~ transform2[2,3] += 0.105
    ine = rotate_inequalities(__ineq_lf_in_rf_hrp2, transform2)
    return (ine.A, ine.b)
    

# add foot offset
def right_foot_in_lf_frame_talos_constraints(transform):
    global __ineq_rf_in_rl_hrp2
    if __ineq_rf_in_rl_hrp2 is None:
        filekin = insdir +"talos-rbprm/relative_effector_positions/talos_RF_constraints_in_LF_quasi_flat_REDUCED.obj"
        obj = load_obj(filekin)
        __ineq_rf_in_rl_hrp2 = as_inequalities(obj)
    ine = rotate_inequalities(__ineq_rf_in_rl_hrp2, transform.copy())
    return (ine.A, ine.b)
        
def left_foot_in_rf_frame_talos_constraints(transform):
    global __ineq_lf_in_rf_hrp2
    if __ineq_lf_in_rf_hrp2 is None:
        filekin = insdir +"talos-rbprm/relative_effector_positions/talos_LF_constraints_in_RF_quasi_flat_REDUCED.obj"
        obj = load_obj(filekin)
        __ineq_lf_in_rf_hrp2 = as_inequalities(obj)
    ine = rotate_inequalities(__ineq_lf_in_rf_hrp2, transform.copy())
    return (ine.A, ine.b)
    
    
#~ def projectAlongAxis(S, n):
    
    
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
    
if __name__ == '__main__':
    from random import random

    #test inequalities computation
    S_flat = [[0.3, 0.6, 0.15],[0.3, -0.16, 0.15],[0.6, -0.16, 0.15],[0.6, 0.6, 0.15]]
    normal = [0.,0.,1.]
    
    #testing surfacePointsToIneq
    A, b =  surfacePointsToIneq(array(S_flat).T , normal)
    
    for _ in range(100):
        l = len(S_flat)
        a = [random() for _ in range(l)]
        total = sum(a)
        a = [el / total for el in a]
        weighted = [array(S_flat[i]) * a[i] for i in range(l)]
        pt = sum(weighted)
        #~ if (A[:-2,:].dot(pt) - b[:-2]).max() > 0.:
        if (A.dot(pt) - b).max() > 0.00001:
            print "fail", (A.dot(pt) - b).max(), " pt ", pt
            break
    
    
    #testing qhull
    a = array
    pts = [a([1,2]),a([4,-1]),a([-2,2]),a([1,3])]
    apts = a(pts)
    
    c = ConvexHull(apts)
    A = c.equations[:,:-1]
    b = -c.equations[:,-1]

    for i in range(100):
        l = len(pts)
        a = [random() for _ in range(l)]
        total = sum(a)
        a = [el / total for el in a]
        weighted = [pts[i] * a[i] for i in range(l)]
        pt = sum(weighted)
        if (A.dot(pt) - b).max() > 0.:
            print "fail", (A.dot(pt) + b).max()
           
           


####################### BENCH ###################"

from time import  clock
def timMs(t1, t2):
    return (t2-t1) * 1000.
    
    
if __name__ == '__main__':
    #testing inequalities
    def default_transform_from_pos_normal(pos = array([0.,0.,0.]), normal = array([0.,0,1.])):
        f = array([0.,0.,1.])
        t = array(normal)
        t = t / norm(t)
        v = np.cross(f, t)
        c = np.dot(f, t)
        if c > 0.99:
            rot = identity(3)    
        else:
            u = v/norm(v)
            h = (1. - c)/(1. - c**2)
            vx, vy, vz = v
            rot =array([[c + h*vx**2, h*vx*vy - vz, h*vx*vz + vy],
                  [h*vx*vy+vz, c+h*vy**2, h*vy*vz-vx],
                  [h*vx*vz - vy, h*vy*vz + vx, c+h*vz**2]])
        return vstack( [hstack([rot,pos.reshape((-1,1))]), [ 0.        ,  0.        ,  0.        ,  1.        ] ] )
    
    tr = default_transform_from_pos_normal()
    A, b = left_foot_hrp2_constraints(tr)
    Ar, br = left_foot_hrp2_constraints_reduced(tr)
    
    from plot_plytopes import plot_polytope_H_rep    
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(131, projection="3d")
    plot_polytope_H_rep(A,b.reshape((-1,1)), color = "r", just_pts = False, ax = ax)
    ax = fig.add_subplot(132, projection="3d")
    plot_polytope_H_rep(Ar,br.reshape((-1,1)), color = "r", just_pts = False, ax = ax)
    ax = fig.add_subplot(133, projection="3d")
    plot_polytope_H_rep(vstack([A,Ar]),concatenate([b, br]).reshape((-1,1)), color = "r", just_pts = False, ax = ax)
    plt.show()
