import numpy as np
from numpy import arange, array, append, cross, dot, zeros
from numpy.linalg import norm

from scipy.spatial import ConvexHull

def normal(points):
    p1 = array(points[0])
    p2 = array(points[1])
    p3 = array(points[2])
    normal = cross((p2 - p1),(p3 - p1))
    normal /= norm(normal)
    return normal.tolist()

def cutList2D(l):
    return [el[:2] for el in l]
    
#precision handling     
def roundPoints (points, precision):  
    return [[round(x,precision) for x in p] for p in points]    
    
def removeDuplicates (points):
    pList = []
    for p in points:
        if p not in pList:
            pList.append(p)
    return pList
    
def computeAxisAngleRotation(u, c):
    ux = u[0] ; uy = u[1] ; uz = u[2]
    s = np.sqrt(1- c*c)
    return [[c+ux*ux*(1-c), ux*uy*(1-c)-uz*s, ux*uz*(1-c)+uy*s],
    [uy*ux*(1-c)+uz*s,    c+uy*uy*(1-c), uy*uz*(1-c)-ux*s],
    [uz*ux*(1-c)-uy*s, uz*uy*(1-c)+ux*s,    c+uz*uz*(1-c)]]

def getSurfaceRotation (surface):
    n = surface[1]
    cosx = np.dot(surface[1],[0,0,1])
    axisx = np.cross(surface[1],[0,0,1]) ; 
    n_axisx = norm(axisx)
    if n_axisx > 0:
      axisx /= n_axisx
    return computeAxisAngleRotation(axisx, cosx)
    
def getPtsRotation (points):
    return getSurfaceRotation((points,normal(points)))

def getSurfaceTranslation (surface):
    return [sum(x)/len(x) for x in zip(*surface[0])]     
    
def getPtsTranslation (points):
    return getSurfaceTranslation((points,normal(points)))

#into xy plane, back to its position    
def allignSurface (surface):        
    R = getSurfaceRotation(surface)
    t = getSurfaceTranslation(surface)
    translatedPts = [(array(p)-array(t)).tolist() for p in surface[0]]
    rotatedPts = [np.dot(R,p).tolist() for p in translatedPts]
    return [(array(p)+array(t)).tolist() for p in rotatedPts]  

def allignPoints (points):
    return allignSurface((points, normal(points)))
    
def pointsTransform (points,R,t):
    translatedPts = [(array(pt)-array(t)).tolist() for pt in points]
    rotatedPts = [np.dot(R,pt).tolist() for pt in translatedPts]
    return [(array(pt)+array(t)).tolist() for pt in rotatedPts]  
    
def getSurfaceExtremumPoints(el):
    pts = removeDuplicates(el[0]+el[1])
    apts = allignPoints(pts)
    hull = ConvexHull(cutList2D(apts))#,False,'QbB')
    return [pts[idx] for idx in hull.vertices]
    
        
#get contact surfaces (pts and normal)
def contactSurfaces(afftool):
    l = afftool.getAffordancePoints("Support")
    return [(getSurfaceExtremumPoints(el), normal(el[0])) for el in l]
    
#get collding affordance index 
def getCollidingAffIndex(contactPolygon,surfaces):  
    for index, surface in enumerate(surfaces):
        R = getSurfaceRotation(surface)
        t = getSurfaceTranslation(surface)
        transformedPts = pointsTransform(contactPolygon, R, t)
        pts = allignSurface(surface) + transformedPts
        #pts = cutList2D(pts) 
        pts = cutList2D(roundPoints(pts,3)) 
        hull = ConvexHull(pts)#,False,'QbB')
        if sorted(hull.vertices) == list(range(len(surface[0]))): return index        
    return -1

    
######################################################################### 
"""
pId = tp.ps.numberPaths() -1
pathLength = tp.ps.pathLength(pId) #length of the path
discretisationStep = 0.1
discretisationSample = 0.05
configs = []

# get configuration along the path
for s in arange (0, pathLength, discretisationStep) :
    configs.append(tp.ps.configAtParam(pId, s)) 

#get surface information
surfaces = contactSurfaces(afftool) 
surfaces = sorted(surfaces) #sort the planes in order
del surfaces[6] #remove handRail from the list

#get surface candidate at each discretization step
sequences = []
for config in configs:
    contacts = rbprmBuilder.getContactSurfacesAtConfig(config,'hrp2_rleg_rom')  #right leg
    sequence = []
    for contact in contacts:
        index = getCollidingAffIndex(contact, surfaces)
        sequence.append(surfaces[index])
    sequences.append(sequence)
    
pLFs = [] ; pRFs = [] ; ps =[]
nLFs = [] ; nRFs = [] ; ns = []
for config in configs:
    contactsLF = rbprmBuilder.getContactSurfacesAtConfig(config,'hrp2_lleg_rom')  #left leg
    contactsRF = rbprmBuilder.getContactSurfacesAtConfig(config,'hrp2_rleg_rom')  #right leg
    pLF = [] ; pRF = [];
    nLF = [] ; nRF = [];
    
    for contact in contactsLF:
        index = getCollidingAffIndex(contact, surfaces)
        pLF.append(surfaces[index][0])
        nLF.append(surfaces[index][1])
    pLFs.append(pLF)
    nLFs.append(nLF)
    
    for contact in contactsRF:
        index = getCollidingAffIndex(contact, surfaces)
        pRF.append(surfaces[index][0])
        nRF.append(surfaces[index][1])
    pRFs.append(pRF)
    nRFs.append(nRF)

seqsLF = [] ; seqsRF = [] ; seqs = []
seqsLF.append(pLFs) ; seqsLF.append(nLFs)
seqsRF.append(pRFs) ; seqsRF.append(nRFs)
seqs.append(seqsLF) ; seqs.append(seqsRF)

#get rotation matrix of the root at each discretization step
from pinocchio import XYZQUATToSe3
R = []
for config in configs:
    v = config[3:7]
    R.append(XYZQUATToSe3([0,0,0]+v).rotation)

"""
