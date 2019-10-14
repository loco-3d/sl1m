import numpy as np
from hpp_centroidal_dynamics import *
from hpp_spline import *
from numpy import array, asmatrix, matrix, zeros, ones
from numpy import array, dot, vstack, hstack, asmatrix, identity, cross
from numpy.linalg import norm

from scipy.spatial import ConvexHull
from hpp_bezier_com_traj import *
#~ from qp import solve_lp

import eigenpy
import cdd
from curves import bezier3
from random import random as rd
from random import randint as rdi
from numpy import squeeze, asarray


Id = matrix([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
g = array([0.,0.,-9.81])
