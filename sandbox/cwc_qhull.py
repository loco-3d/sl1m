"""
Created on Thurs Aug 4

@author: adelpret, updated by stonneau
"""

import sys
sys.path.insert(0, './tools')

from polytope_conversion_utils import *
from transformations import euler_matrix
from numpy import array, vstack, zeros, sqrt
import numpy as np

from math import cos, sin, tan, atan, pi
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import axes
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull

NUMBER_TYPE = 'float'  # 'float' or 'fraction'

                     

n = 3;      # generator size
cg = 4;     # number of generators per contact
#~ USE_DIAGONAL_GENERATORS = True;
USE_DIAGONAL_GENERATORS = False;
CONTACT_SET = 1;

## 
#  Given a list of contact points
#  as well as a list of associated normals
#  compute the gravito inertial wrench cone
#  \param p array of 3d contact positions
#  \param N array of 3d contact normals
#  \param simplify_cones if true inequality conversion will try to remove 
#  redundancies
#  \param params requires "mass" "g"  and "mu"
#  \return the CWC H, H w <= 0, where w is the wrench
def compute_CWC(p, N,  params):
	#~ print "CWC"
	''' compute generators '''
	mass = params["mass"]
	g = params["g"]
	mu = params["mu"]
	c = p.shape[0];
	#gamma = atan(mu);   # half friction cone angle
	m = c*cg;            # number of generators
	S = np.zeros((n,m));
	T1 = np.zeros((c,n));
	T2 = np.zeros((c,n));
	muu = mu/sqrt(2);
	for i in range(c):
		''' compute tangent directions '''
		N[i,:]  = N[i,:]/np.linalg.norm(N[i,:]);
		T1[i,:] = np.cross(N[i,:], [0,1,0]);
		if(np.linalg.norm(T1[i,:])<1e-5):
			T1[i,:] = np.cross(N[i,:], [1,0,0]);
		T1[i,:] = T1[i,:]/np.linalg.norm(T1[i,:]);
		T2[i,:] = np.cross(N[i,:], T1[i,:]);
		T2[i,:] = T2[i,:]/np.linalg.norm(T2[i,:]);
		
		if(USE_DIAGONAL_GENERATORS):
			S[:,cg*i+0] =  muu*T1[i,:] + muu*T2[i,:] + N[i,:];
			S[:,cg*i+1] =  muu*T1[i,:] - muu*T2[i,:] + N[i,:];
			S[:,cg*i+2] = -muu*T1[i,:] + muu*T2[i,:] + N[i,:];
			S[:,cg*i+3] = -muu*T1[i,:] - muu*T2[i,:] + N[i,:];
		else:
			S[:,cg*i+0] =   mu*T1[i,:] + N[i,:];
			S[:,cg*i+1] =  -mu*T1[i,:] + N[i,:];
			S[:,cg*i+2] =   mu*T2[i,:] + N[i,:];
			S[:,cg*i+3] = - mu*T2[i,:] + N[i,:];
		
		S[:,cg*i+0] = S[:,cg*i+0]/np.linalg.norm(S[:,cg*i+0]);
		S[:,cg*i+1] = S[:,cg*i+1]/np.linalg.norm(S[:,cg*i+1]);
		S[:,cg*i+2] = S[:,cg*i+2]/np.linalg.norm(S[:,cg*i+2]);
		S[:,cg*i+3] = S[:,cg*i+3]/np.linalg.norm(S[:,cg*i+3]);

	''' compute matrix mapping contact forces to gravito-inertial wrench '''
	M = np.zeros((6,3*c));
	for i in range(c):
		M[:3, 3*i:3*i+3] = -np.identity(3);
		M[3:, 3*i:3*i+3] = -crossMatrix(p[i,:]);
		
	''' project generators in 6d centroidal space '''
	S_centr = np.zeros((6,m));
	for i in range(c):
		S_centr[:,cg*i:cg*i+cg] = np.dot(M[:,3*i:3*i+3], S[:,cg*i:cg*i+cg]);
	''' convert generators to inequalities '''
	#~ H = cone_span_to_face(S_centr, simplify_cones); 
    #~ St = S.T
    #~ pts = [S[:,i].reshape((-1,)) for i in range(D)]
	hull = ConvexHull(S_centr.T)

	print "CWC_END", hull.equations[:,:-1].shape
	return hull.equations[:,:-1], -hull.equations[:,-1].reshape((-1))
	
	#~ ''' project inequalities from 6d to 3d com space '''
	#~ g_vector = np.zeros(3);
	#~ g_vector[2] = -g;
	#~ D = np.zeros((6,3));
	#~ d = np.zeros((6,1));
	#~ D[3:,:] = -mass*crossMatrix(g_vector);
	#~ d[:3]   = mass*g_vector.reshape(3,1);
	#~ F_com = np.dot(H, D);
	#~ f_com = np.dot(H, d);
	
	#~ ''' plot inequalities on x-y plane '''
	#~ f, ax = plt.subplots();
	#~ com_x = np.zeros(2);
	#~ com_y = np.zeros(2);
	#~ com = np.zeros(3);
	#~ com[2] = 0;     # com height
	#~ for i in range(F_com.shape[0]):
		#~ if(np.abs(F_com[i,1])>1e-13):
			#~ com_x[0] = -1;   # com x coordinate
			#~ com_x[1] = 4;   # com x coordinate
			#~ com[0] = com_x[0];
			#~ com[1] = 0;
			#~ com_y[0] = (-f_com[i] - np.dot(F_com[i,:],com) )/F_com[i,1];
			#~ 
			#~ com[0] = com_x[1];
			#~ com_y[1] = (-f_com[i] - np.dot(F_com[i,:],com) )/F_com[i,1];
#~ 
			#~ ax.plot(com_x, com_y, 'k-');
		#~ elif(np.abs(F_com[i,0])>1e-13):
			#~ com_y[0] = -4;
			#~ com_y[1] = 4;
			#~ com[0] = 0;
			#~ com[1] = com_y[0];
			#~ com_x[0] = (-f_com[i] - np.dot(F_com[i,:],com) )/F_com[i,0];
#~ 
			#~ com[1] = com_y[1];
			#~ com_x[1] = (-f_com[i] - np.dot(F_com[i,:],com) )/F_com[i,0];
			#~ ax.plot(com_x, com_y, 'k-');
		#~ else:
			#~ print "[WARNING] Could not print one inequality as all coefficients are 0: F_com[%d,:]=[%f,%f]" % (i,F_com[i,0],F_com[i,1]);
			#~ 
	#~ ax.set_xlim([-1, 1]);
	#~ ax.set_ylim([-2, 0]);
#~ 
	#~ ''' test equilibrium for a grid of points '''
	#~ x_range = np.arange(-1, 1, 0.05);
	#~ y_range = np.arange(-2, 0, 0.05);
	#~ equilibrium = np.zeros((len(x_range), len(y_range)));
	#~ w = np.zeros(6);
	#~ com = np.zeros(3);
	#~ w[2] = -mass*g;
	#~ for i in range(len(x_range)):
		#~ for j in range(len(y_range)):
			#~ com[0] = x_range[i];
			#~ com[1] = y_range[j];
			#~ w[3:] = mass*np.cross(com, g_vector);
			#~ equilibrium[i,j] = np.max(np.dot(H, w));
			#~ 
	#~ ''' plot equilibrium points in blue, other points in red '''
	#~ for i in range(len(x_range)):
		#~ for j in range(len(y_range)):
			#~ if(equilibrium[i,j]<=1e-12):
				#~ ax.scatter(x_range[i], y_range[j], c='b');
			#~ else:
				#~ ax.scatter(x_range[i], y_range[j], c='r');
	#~ ''' plot contact points in black '''
	#~ for i in range(c):
		#~ ax.scatter(p[i,0], p[i,1], c='k', s=100);
#~ 
	#~ #fig = plt.figure(); ax = fig.add_subplot(111, projection='3d', aspect='equal');
	#~ fig = plt.figure(figsize=plt.figaspect(0.5)*1.5) #Adjusts the aspect ratio and enlarges the figure (text does not enlarge)
	#~ ax = fig.gca(projection='3d')
	#~ line_styles =["b", "r", "c", "g"];
	#~ ss = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3];
	#~ for i in range(c):
		#~ ax.scatter(p[i,0],p[i,1],p[i,2], c=line_styles[i%len(line_styles)], marker='x');
		#~ for s in ss:
	#~ #        ax.scatter(p[i,0]+s*N[i,0],p[i,1]+s*N[i,1],p[i,2]+s*N[i,2], c=line_styles[i%len(line_styles)], marker='o');
			#~ for j in range(cg):
				#~ ax.scatter(p[i,0]+s*S[0,cg*i+j],p[i,1]+s*S[1,cg*i+j],p[i,2]+s*S[2,cg*i+j], c=line_styles[i%len(line_styles)], marker='o');
	#~ ax.set_xlabel('x');
	#~ ax.set_ylabel('y');
	#~ ax.set_zlabel('z');
	#~ #axisEqual3D(ax);
#~ 
	#~ plt.show();


