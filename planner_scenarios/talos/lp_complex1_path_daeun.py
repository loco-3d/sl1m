from hpp.corbaserver.rbprm.rbprmbuilder import Builder
from hpp.gepetto import Viewer
from hpp.corbaserver import Client
from hpp.corbaserver.rbprm.talos_abstract import Robot

packageName = 'hpp_environments'
meshPackageName = 'hpp_environments'


rbprmBuilder = Robot ()
rbprmBuilder.setJointBounds ("root_joint", [-2.8,1.8,0.19,0.21, 0.95,1.6])
vMax = 0.3# linear velocity bound for the root
aMax = 0.1 # linear acceleration bound for the root
extraDof = 6
mu=0.5# coefficient of friction

rbprmBuilder.setFilter([Robot.rLegId,Robot.lLegId])

rbprmBuilder.setAffordanceFilter(Robot.rLegId, ['Support',])
rbprmBuilder.setAffordanceFilter(Robot.lLegId, ['Support'])
rbprmBuilder.boundSO3([-4.,4.,-0.1,0.1,-0.1,0.1])
# Add 6 extraDOF to the problem, used to store the linear velocity and acceleration of the root
rbprmBuilder.client.robot.setDimensionExtraConfigSpace(extraDof)
# We set the bounds of this extraDof with velocity and acceleration bounds (expect on z axis)
extraDofBounds = [-vMax,vMax,-vMax,vMax,-10.,10.,-aMax,aMax,-aMax,aMax,-10.,10.]
rbprmBuilder.client.robot.setExtraConfigSpaceBounds(extraDofBounds)
indexECS = rbprmBuilder.getConfigSize() - rbprmBuilder.client.robot.getDimensionExtraConfigSpace()


#~ from hpp.corbaserver.rbprm. import ProblemSolver
from hpp.corbaserver.problem_solver import ProblemSolver
ps = ProblemSolver( rbprmBuilder )

v = Viewer (ps)
v.addLandmark(v.sceneName,1)
from hpp.corbaserver.affordance.affordance import AffordanceTool
afftool = AffordanceTool ()
afftool.setAffordanceConfig('Support', [0.5, 0.03, 0.00005])
afftool.loadObstacleModel (packageName, "multicontact/bauzil_ramp_simplified", "planning", v,reduceSizes=[0.,0.,0.])
afftool.visualiseAffordances('Support', v, [0.25, 0.5, 0.5])




ps.setParameter("Kinodynamic/velocityBound",vMax)
ps.setParameter("Kinodynamic/accelerationBound",aMax)
# force the orientation of the trunk to match the direction of the motion
ps.setParameter("Kinodynamic/forceYawOrientation",True)
ps.setParameter("DynamicPlanner/sizeFootX",0.2)
ps.setParameter("DynamicPlanner/sizeFootY",0.12)
ps.setParameter("DynamicPlanner/friction",mu)
# sample only configuration with null velocity and acceleration :
ps.setParameter("ConfigurationShooter/sampleExtraDOF",False)
ps.setParameter("PathOptimization/RandomShortcut/NumberOfLoops",100)

q_init = rbprmBuilder.getCurrentConfig ();
q_init [0:3] = [-1, -0.82, 0.5]; rbprmBuilder.setCurrentConfig (q_init); r (q_init)

q_goal = q_init [::]
q_goal [3:7] = [ 0.,  0.        ,  0.14943813,  0.98877108        ]
q_goal [0:3] = [100.49, -0.65, 1.2]; r (q_goal)
ps.addPathOptimizer("RandomShortcut")
ps.setInitialConfig (q_init)
ps.addGoalConfig (q_goal)
r(q_goal)

from hpp.corbaserver.affordance.affordance import AffordanceTool
afftool = AffordanceTool ()
afftool.setAffordanceConfig('Support', [0.5, 0.03, 0.00005])
afftool.loadObstacleModel (packageName, "multicontact/bauzil_ramp", "planning", r)
afftool.visualiseAffordances('Support', r, [0.25, 0.5, 0.5])

ps.client.problem.selectConfigurationShooter("RbprmShooter")
ps.client.problem.selectPathValidation("RbprmPathValidation",0.05)
#~ t = ps.solve ()



from hpp.gepetto import PathPlayer
pp = PathPlayer (v)
pp.dt=0.1
pp.displayVelocityPath(1)
v.client.gui.setVisibility("path_1_root","ALWAYS_ON_TOP")
pp.dt = 0.01
  

q_far = q_goal[::]
q_far[2] = -5
v(q_far)

######################################################################### 

pId = ps.numberPaths() -1
pathLength = ps.pathLength(pId) #length of the path
discretisationStep = 0.8
configs = []

from numpy import arange
from sl1m.test import *
from tools.display_tools import displaySurfaceFromPoints

# get configuration along the path
for s in arange (0, pathLength, discretisationStep) :
    configs.append(ps.configAtParam(pId, s)) 

# get surface information
surfaces = contactSurfaces(afftool) 
surfaces = sorted(surfaces) #sort the planes in order

del surfaces[0]
floor1  = [[-2.45, 0.54  , 0.  ], [-3.59 ,  0.54, 0. ], [-3.59 , -0.15, 0.  ], [-2.45, -0.15, 0.  ], ]
floor2  = [[-0.35, 0.54  , 0.  ], [0.01 ,  0.54, 0. ], [0.01 , -0.46, 0.  ], [-0.35, -0.46, 0.  ], ]
surfaces.append((floor1,normal(floor1)))
surfaces.append((floor2,normal(floor2)))
#del surfaces[16] #part of the bridge
#bridge = [[ 1.51, -0.46 , 0.6 ], [1.51 , -0.76, 0.6], [-1.49, -0.76, 0.6 ], [-1.49, -0.46, 0.6 ], ]
#surfaces.append((bridge,normal(bridge)))
del surfaces[0] #big bottom plane

# get surface candidate at each discretization step
# suppose that se start with the left leg
seqs = []; indices = []
for i, config in enumerate(configs):    
    seq = [] ; index_list = []
    if i % 2 == 0 : # left leg
        contacts = rbprmBuilder.getContactSurfacesAtConfig(config,'talos_lleg_rom')  
        for contact in contacts:
            if contact != []:
                index = getCollidingAffIndex(contact, surfaces)
                index_list.append(index)
                if index != -1 : seq.append(array(surfaces[index][0]).T)
                #displaySurfaceFromPoints(v,contact,[0,0,1,1])
    else : # right leg
        contacts = rbprmBuilder.getContactSurfacesAtConfig(config,'talos_rleg_rom')  
        for contact in contacts:
            if contact != []:
                index = getCollidingAffIndex(contact, surfaces)
                index_list.append(index)
                if index != -1 :seq.append(array(surfaces[index][0]).T)
                #displaySurfaceFromPoints(v,contact,[0,0,1,1])
    seqs.append(seq)
    indices.append(index_list)


"""
#get surface candidate at each discretization step
seqs = [] ; seqsLF = [] ; seqsRF = []
for config in configs:
    contactsLF = rbprmBuilder.getContactSurfacesAtConfig(config,'hrp2_lleg_rom')  #left leg
    contactsRF = rbprmBuilder.getContactSurfacesAtConfig(config,'hrp2_rleg_rom')  #right leg
    seqLF = [] ; seqRF = []
    
    for contact in contactsLF:
        index = getCollidingAffIndex(contact, surfaces)
        seqLF.append(array(surfaces[index][0]).T)
    
    for contact in contactsRF:
        index = getCollidingAffIndex(contact, surfaces)
        seqRF.append(array(surfaces[index][0]).T)
    
    seqsLF.append(seqLF)
    seqsRF.append(seqRF)

seqs.append(seqsLF)    
seqs.append(seqsRF)    
"""

#get rotation matrix of the root at each discretization step
from pinocchio import XYZQUATToSe3
R = []
for config in configs:
    v = config[3:7]
    R.append(XYZQUATToSe3([0,0,0]+v).rotation)



