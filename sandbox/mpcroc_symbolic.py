from symbolic import *

#problem: start with known position, velocity, acceleration (3 variables); comes to a stop
# at the  end of motion

deg = 5
pts = [p_is[1], p_is[2], p_is[3], y, y, y]
c = Bezier(5, pts)

(wdc,wuc) =  computewUpDownWaypoints(c)

#curve spans 3 contact phases, and is thus split at t1 and t2
#c goes from 0 to T

t1, t2 = symbols('t1 t2')
u1, u2 = symbols('u1 u2') 
# normalized times
#~ us = (TSplit-Tmin)/(Tmax-Tmin);
# u1 = (t1-0 )/(T-0)  = t1 / T;
# u2 = ((t2-t1)-t1)/(T- t1) = (t2 - 2*t1) / (T-t1);

c01, c1T = split_normalized(c,u1)
c12, c2T = split_normalized(c1T,u2)

#~ (wdc01,wuc01) =  computewUpDownWaypoints(c01)
#~ (wdc12,wuc12) =  computewUpDownWaypoints(c12)
#~ (wdc2T,wuc2T) =  computewUpDownWaypoints(c2T)

# now we have three curves of degree 6
#each waypoint is going to be abstract as a linear function of y
currentIndexp_is = deg - 1

equivalence = {}

#~ pi01 = []
#~ for (py, pc) in c01.wps_as_factor(y):
    #~ currentIndexp_is += 1;
    #~ npc = p_is[currentIndexp_is]
    #~ equivalence [npc] = pc
    #~ pi01 += [py * y + npc]
    
#~ c01 = Bezier(deg, pi01)
#~ (wdc01,wuc01) =  computewUpDownWaypoints(c01)
