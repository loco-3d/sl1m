
from numpy import arange
from mpcroc.test import contactSurfaces,getCollidingAffIndex,removeDuplicates
from tools.display_tools import displaySurfaceFromPoints
from pinocchio import XYZQUATToSe3
import numpy as np
from numpy import array
from mpcroc.problem_definition import LF

MAX_SURFACE = 0.3 # if a contact surface is greater than this value, the intersection is used instead of the whole surface

def rotationMatrixFromConfigs(configs):
  R = []
  for config in configs:
      q_rot = config[3:7]
      R.append(XYZQUATToSe3([0,0,0]+q_rot).rotation)
  return R     
  
def listToArray (seqs):
    nseq = []; nseqs= []
    for seq in seqs:
        nseq = []
        for surface in seq:
            nseq.append(array(surface).T)
        nseqs.append(nseq)
    return nseqs
  return R

def area(s):
    #print "in area, s = ",s
    area = 0
    for i in range(1,len(s)-1):
      #area += x[i]*( y[i+1] - y[i-1] );
      area += abs(s[i][0]*(s[i+1][1] - s[i-1][1]))
    i = len(s) -1 
    area += abs(s[i][0]*(s[0][1] - s[i-1][1]))
    #print "area = ",area*0.5
    return area * 0.5

def getContactsNames(rbprmBuilder,i,q):
  if i % 2 == LF : # left leg moving
    step_contacts = rbprmBuilder.clientRbprm.rbprm.getCollidingObstacleAtConfig(q,'talos_lleg_rom') 
  else :
    step_contacts = rbprmBuilder.clientRbprm.rbprm.getCollidingObstacleAtConfig(q,'talos_rleg_rom')
  return step_contacts

def getContactsIntersections(rbprmBuilder,i,q):
    if i % 2 == LF : # left leg
        intersections = rbprmBuilder.getContactSurfacesAtConfig(q,'talos_lleg_rom')  
    else : # right leg
        intersections = rbprmBuilder.getContactSurfacesAtConfig(q,'talos_rleg_rom')
    return intersections


def getSurfacesFromGuideContinuous(rbprmBuilder,ps,afftool,pId,viewer = None,step = 1.,useIntersection= False):
  pathLength = ps.pathLength(pId) #length of the path
  discretizationStep = 0.5 # step at which we check the colliding surfaces
  #print "path length = ",pathLength
  # get surface information
  all_surfaces = contactSurfaces(afftool) 
  all_names = afftool.getAffRefObstacles("Support") # id in names and surfaces match
  surfaces_dict = dict(zip(all_names, all_surfaces)) # map surface names to surface points
  seqs = [] # list of list of surfaces : for each phase contain a list of surfaces. One phase is defined by moving of 'step' along the path
  t = 0.
  current_phase_end = step
  end = False
  i = 0
  while not end: # for all the path
    #print "Looking for surfaces for phase "+str(len(seqs))+" for t in ["+str(t)+" ; "+str(current_phase_end)+" ] "
    phase_contacts_names = []
    while t < current_phase_end: # get the names of all the surfaces that the rom collide while moving from current_phase_end-step to current_phase_end
      q = ps.configAtParam(pId, t)
      step_contacts = getContactsNames(rbprmBuilder,i,q)
      for contact_name in step_contacts : 
        if not contact_name in phase_contacts_names:
          phase_contacts_names.append(contact_name)
      t += discretizationStep
    # end current phase
    # get all the surfaces from the names and add it to seqs: 
    if useIntersection : 
      intersections = getContactsIntersections(rbprmBuilder,i,q)
    phase_surfaces = []
    for name in phase_contacts_names:
      surface = surfaces_dict[name][0]
      if useIntersection and area(surface) > MAX_SURFACE : 
        if name in step_contacts : 
          intersection = intersections[step_contacts.index(name)]
          phase_surfaces.append(intersection)
      else :
        phase_surfaces.append(surface) # [0] because the last vector contain the normal of the surface
    #print "There was "+str(len(phase_contacts_names))+" surfaces in contact during this phase."
    phase_surfaces = sorted(phase_surfaces) # why is this step required ? without out the lp can fail
    phase_surfaces_array = [] # convert from list to array, we cannot do this before because sorted() require list
    for surface in phase_surfaces:
      phase_surfaces_array.append(array(surface).T)
    #print "phase_surfaces_array = ",phase_surfaces_array    
    seqs.append(phase_surfaces_array)

    # increase values for next phase
    t = current_phase_end
    i += 1 
    if current_phase_end == pathLength:
      end = True
    current_phase_end += step
    if current_phase_end >= pathLength:
      current_phase_end = pathLength
  # end for all the guide path

  #get rotation matrix of the root at each discretization step
  configs = []
  for t in arange (0, pathLength, step) :
    configs.append(ps.configAtParam(pId, t)) 
  R = rotationMatrixFromConfigs(configs)
  return R,seqs


def getSurfacesFromGuide(rbprmBuilder,ps,afftool,pId,viewer = None,discretisationStep = 0.6,useIntersection = False, mergeCandidates = False):
  pathLength = ps.pathLength(pId) #length of the path
  configs = []
  # get configuration along the path
  for s in arange (0, pathLength, discretisationStep) :
      configs.append(ps.configAtParam(pId, s)) 

  # get surface information
  surfaces = contactSurfaces(afftool) 
  surfaces = sorted(surfaces) #sort the planes in order
  
  # get surface candidate at each discretization step
  # suppose that we start with the left leg
  seqs = []; #indices = []
  for i, config in enumerate(configs):    
      seq = [] ; #index_list = []
      if i % 2 == 0 : # left leg
          contacts = rbprmBuilder.getContactSurfacesAtConfig(config,'talos_lleg_rom')  
      else : # right leg
          contacts = rbprmBuilder.getContactSurfacesAtConfig(config,'talos_rleg_rom')  
      for contact in contacts:
          if contact != []:
              if useIntersection:
                seq.append(contact) # changed it to list to make it easier to remove duplicates
                #seq.append(array(contact).T)
              else:
                index = getCollidingAffIndex(contact, surfaces)
                #index_list.append(index)
                #if index != -1 : seq.append(array(surfaces[index][0]).T)
                if index != -1 : seq.append(surfaces[index][0]) # changed it to list to make it easier to remove duplicates
              if viewer:
                displaySurfaceFromPoints(viewer,contact,[0,0,1,1])
      seqs.append(seq)
      #indices.append(index_list)
  
  # merge candidates with the previous and the next phase
  if mergeCandidates:
      nseqs = []
      for i, seq in enumerate(seqs):
          nseq = []
          if i == 0: nseq=seqs[i]+seqs[i+1]
          elif i == len(seqs)-1: nseq=seqs[i-1]+seqs[i]
          else: nseq=seqs[i-1]+seqs[i]+seqs[i+1]
          nseqs.append(nseq)  
      seqs = nseqs
  
  # remove duplicates          
  for i, seq in enumerate(seqs): seqs[i] = removeDuplicates(seq)
  
  seqs = listToArray(seqs) # change the format from list to array
  R = rotationMatrixFromConfigs(configs)
  return R,seqs


"""
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
                if viewer:
                  displaySurfaceFromPoints(viewer,contact,[0,0,1,1])
    else : # right leg
        contacts = rbprmBuilder.getContactSurfacesAtConfig(config,'talos_rleg_rom')  
        for contact in contacts:
            if contact != []:
                index = getCollidingAffIndex(contact, surfaces)
                index_list.append(index)
                if index != -1 :seq.append(array(surfaces[index][0]).T)
                if viewer:
                  displaySurfaceFromPoints(viewer,contact,[0,0,1,1])
    seqs.append(seq)
    indices.append(index_list)
"""
