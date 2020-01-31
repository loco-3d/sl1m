from numpy import arange, array
from narrow_convex_hull import getSurfaceExtremumPoints, removeDuplicates, normal, area
from tools.display_tools import displaySurfaceFromPoints
from pinocchio import XYZQUATToSe3
import numpy as np

ROBOT_NAME = 'talos'
MAX_SURFACE = 0.3 # if a contact surface is greater than this value, the intersection is used instead of the whole surface
LF = 0
RF = 1  

def removeNull (l):
  ll = []
  for el in l:
    if el != []:
      ll.append(el)
  return ll

# change the format into an array  
def listToArray (seqs):
  nseq = []; nseqs= []
  for seq in seqs:
    nseq = []
    for surface in seq:
      if surface != []:
        nseq.append(array(surface).T)
    nseqs.append(nseq)
  return nseqs

# get configurations along the path  
def getConfigsFromPath (ps, pathId = 0, discretisationStepSize = 1.) :
  configs = []
  pathLength = ps.pathLength(pathId)
  for s in arange (0, pathLength, discretisationStepSize) :
    configs.append(ps.configAtParam(pathId, s))
  return configs

# get all the contact surfaces (pts and normal)
def getAllSurfaces(afftool) :
  l = afftool.getAffordancePoints("Support")
  return [(getSurfaceExtremumPoints(el), normal(el[0])) for el in l]
    
# get surface information
def getAllSurfacesDict (afftool) :
  all_surfaces = getAllSurfaces(afftool) 
  all_names = afftool.getAffRefObstacles("Support") # id in names and surfaces match
  surfaces_dict = dict(zip(all_names, all_surfaces)) # map surface names to surface points
  return surfaces_dict

# get rotation matrices form configs
def getRotationMatrixFromConfigs(configs) :
  R = []
  for config in configs:
    q_rot = config[3:7]
    R.append(XYZQUATToSe3([0,0,0] + q_rot).rotation)
  return R   
    
# get contacted surface names at configuration
def getContactsNames(rbprmBuilder,i,q):
  if i % 2 == LF : # left leg 
    step_contacts = rbprmBuilder.clientRbprm.rbprm.getCollidingObstacleAtConfig(q, ROBOT_NAME + '_lleg_rom') 
  elif i % 2 == RF : # right leg 
    step_contacts = rbprmBuilder.clientRbprm.rbprm.getCollidingObstacleAtConfig(q, ROBOT_NAME + '_rleg_rom')
  return step_contacts

# get intersections with the rom and surface at configuration
def getContactsIntersections(rbprmBuilder,i,q):
  if i % 2 == LF : # left leg
    intersections = rbprmBuilder.getContactSurfacesAtConfig(q, ROBOT_NAME + '_lleg_rom') 
  elif i % 2 == RF : # right leg
    intersections = rbprmBuilder.getContactSurfacesAtConfig(q, ROBOT_NAME + '_rleg_rom')
  # return intersections
  return removeNull(intersections)

# merge phases with the next phase
def getMergedPhases (seqs):
  nseqs = []
  for i, seq in enumerate(seqs):
    nseq = []
    if i == len(seqs)-1: nseq = seqs[i]
    else: nseq = seqs[i]+seqs[i+1]
    nseq = removeDuplicates(nseq)
    nseqs.append(nseq)  
  return nseqs    



def getSurfacesFromPathContinuous(rbprmBuilder, ps, surfaces_dict, pId, viewer = None, phaseStepSize = 1., useIntersection = False):
    pathLength = ps.pathLength(pId) # length of the path
    discretizationStepSize = 0.4 # step at which we check the colliding surfaces
    
    seqs = [] # list of list of surfaces : for each phase contain a list of surfaces. One phase is defined by moving of 'step' along the path
    t = 0.
    current_phase_end = phaseStepSize
    end = False
    i = 0
    while not end: # for all the path
      phase_contacts_names = []
      while t < current_phase_end: # get the names of all the surfaces that the rom collide while moving from current_phase_end-step to current_phase_end
        q = ps.configAtParam(pId, t)
        step_contacts = getContactsNames(rbprmBuilder,i,q)
        for contact_name in step_contacts : 
          if not contact_name in phase_contacts_names:
            phase_contacts_names.append(contact_name)
        t += discretizationStepSize
      # end current phase
        
      # get all the surfaces from the names and add it to seqs: 
      if useIntersection : 
        intersections = getContactsIntersections(rbprmBuilder,i,q)
            
      phase_surfaces = []
      for name in phase_contacts_names:
        surface = surfaces_dict[name][0] # [0] because the last vector contain the normal of the surface
        if useIntersection and area(surface) > MAX_SURFACE : 
          if len(step_contacts) == len(intersections): # in case of the error with contact detection
            if name in step_contacts : 
              intersection = intersections[step_contacts.index(name)]
              phase_surfaces.append(intersection)
          else:
            phase_surfaces.append(surface)
        else :
          phase_surfaces.append(surface) 

      phase_surfaces = sorted(phase_surfaces) # why is this step required ? without out the lp can fail
      seqs.append(phase_surfaces)

      # increase values for next phase
      t = current_phase_end
      i += 1 
      if current_phase_end == pathLength:
        end = True
      current_phase_end += phaseStepSize
      if current_phase_end >= pathLength:
        current_phase_end = pathLength
    # end for all the guide path
    
    seqs = listToArray(seqs) # convert from list to array, we cannot do this before because sorted() require list

    # get rotation matrix of the root at each discretization step
    configs = []
    for t in arange (0, pathLength, phaseStepSize) :
      configs.append(ps.configAtParam(pId, t)) 
        
    R = getRotationMatrixFromConfigs(configs)
    return R,seqs
    


"""
LARGE_STEP_SIZE = 1.3
SMALL_STEP_SIZE = 1.0

def getSurfacesFromPathContinuous(rbprmBuilder, ps, surfaces_dict, pId, viewer = None, phaseStepSize = 1., useIntersection = False):
    pathLength = ps.pathLength(pId) # length of the path
    discretizationStepSize = 0.5 # step at which we check the colliding surfaces
    
    seqs = [] # list of list of surfaces : for each phase contain a list of surfaces. One phase is defined by moving of 'step' along the path
    t = 0.
    current_phase_end = phaseStepSize
    end = False
    i = 0
    phase_contacts_names = []
    configs = []
    
    while not end: # for all the path
      phase_contacts_prev = phase_contacts_names
      phase_contacts_names = []
      
      configs.append(ps.configAtParam(pId, t)) 
      
      while t < current_phase_end: # get the names of all the surfaces that the rom collide while moving from current_phase_end-step to current_phase_end
        q = ps.configAtParam(pId, t)
        step_contacts = getContactsNames(rbprmBuilder,i,q)
        for contact_name in step_contacts : 
          if not contact_name in phase_contacts_names:
            phase_contacts_names.append(contact_name)
        t += discretizationStepSize
      # end current phase
      # phase_contacts_names = sorted(phase_contacts_names) # is this step needed?
      
      # # if the previous phase and the current phase has different surface candidates, decrease the phase step
      # if phase_contacts_prev != phase_contacts_names:
          # # print "different surface candidates"
          # phaseStepSize = SMALL_STEP_SIZE
      # else:
          # # print "same surface candidates"
          # phaseStepSize = LARGE_STEP_SIZE          
        
      # get all the surfaces from the names and add it to seqs: 
      if useIntersection : 
        intersections = getContactsIntersections(rbprmBuilder,i,q)
        if viewer:
            for intersection in intersections:
                displaySurfaceFromPoints(viewer,intersection,[0,0,1,1])
            
      phase_surfaces = []
      for name in phase_contacts_names:
        surface = surfaces_dict[name][0] # [0] because the last vector contain the normal of the surface
        if useIntersection and area(surface) > MAX_SURFACE : 
          if name in step_contacts : 
            intersection = intersections[step_contacts.index(name)]
            phase_surfaces.append(intersection)
        else :
          phase_surfaces.append(surface) 

      phase_surfaces = sorted(phase_surfaces) # why is this step required ? without out the lp can fail
      seqs.append(phase_surfaces)

      # increase values for next phase
      t = current_phase_end
      i += 1 
      if current_phase_end == pathLength:
        end = True
      current_phase_end += phaseStepSize
      if current_phase_end >= pathLength:
        current_phase_end = pathLength
    # end for all the guide path
    
    seqs = listToArray(seqs) # convert from list to array, we cannot do this before because sorted() require list

    # get rotation matrix of the root at each discretization step
    # configs = []
    # for t in arange (0, pathLength, phaseStepSize) :
      # configs.append(ps.configAtParam(pId, t)) 
        
    R = getRotationMatrixFromConfigs(configs)
    return R,seqs
"""

def getSurfacesFromPath(rbprmBuilder, configs, surfaces_dict, viewer = None, useIntersection = False, useMergePhase = False):
  seqs = [] 
  # get sequence of surface candidates at each discretization step
  for i, q in enumerate(configs):    
    seq = [] 
    intersections = getContactsIntersections(rbprmBuilder,i,q) # get intersections at config
    phase_contacts_names = getContactsNames(rbprmBuilder,i,q) # get the list of names of the surface in contact at config        
    for j, intersection in enumerate(intersections):
      if useIntersection and area(intersection) > MAX_SURFACE : # append the intersection
        seq.append(intersection) 
      else:
        if len(intersections) == len(phase_contacts_names): # in case getCollidingObstacleAtConfig does not work (because of the rom?)
          seq.append(surfaces_dict[phase_contacts_names[j]][0]) # append the whole surface
        else: seq.append(intersection) # append the intersection
      # if viewer:
        # displaySurfaceFromPoints(viewer,intersection,[0,0,1,1])
    seqs.append(seq)
    
  # merge candidates with the previous and the next phase
  if useMergePhase: seqs = getMergedPhases (seqs)
    
  seqs = listToArray(seqs) 
  R = getRotationMatrixFromConfigs(configs)
  return R,seqs
    
