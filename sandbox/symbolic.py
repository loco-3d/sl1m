from __future__ import division
from sympy import *

x, y, z, t = symbols('x y z t')
k, m, n, T = symbols('k m n T', integer=True)
#~ f, g, h = symbols('f g h', cls=Function)


init_printing(use_unicode=True)

from scipy.special import binom as bino

p_isX =  [symbols('Cg')] + [symbols('Cp'+str(i)+'') for i in range(40)]
p_is =   [symbols('g' )] + [symbols('p'+str(i))     for i in range(40)] 
    
class Bezier:
    
    def __init__(self, deg, pis = None, eq = None, mult = None):
        assert(deg)
        self.deg = deg;
        n_pts = deg +1 
        if(pis == None):
            pis = p_is[1:]
        pis = pis[:n_pts]
        self.pis = pis; 
        self.mult = 1;
        if(mult != None):
            self.mult = mult;
        if(eq == None):
            factors = [bino(deg,i) * t**i * (1-t)**(deg-i)*pis[i] for i in range(n_pts)]
            self.eq = sum(factors);
        else:
            self.eq = eq;

    def elevate(self):
        old_deg = self.deg
        new_deg = old_deg+1
        n_pts = new_deg+1
        old_pis =  self.pis
        mult = self.mult
        pis = [old_pis[0]] + [simplify( i / new_deg * old_pis[i-1] + (1 - i / new_deg) * old_pis[i] ) for i in range(1,old_deg+1)] + [old_pis[-1]] 
        factors = [bino(new_deg,i) * t**i * (1-t)**(new_deg-i)*pis[i] for i in range(n_pts)]
        return Bezier(deg = new_deg, eq = sum(factors), pis = pis, mult = mult)
                
    def diff(self):
        old_deg = self.deg
        n_pts = self.deg
        new_deg = n_pts-1
        old_pis =  self.pis
        old_mult = self.mult
        mult = 1 / T * old_mult
        pis = [simplify(old_deg * 1 / T  *(old_pis[i+1] - old_pis[i])) for i in range(n_pts)]
        factors = [bino(new_deg,i) * t**i * (1-t)**(new_deg-i)*pis[i] for i in range(n_pts)]
        return Bezier(deg = new_deg, eq = sum(factors),  pis = pis, mult =  mult)
        
    def wps_as_factor(self, sym = x):
        res = []
        for _, pi in enumerate(self.pis):
            pis_pol = pi.as_poly(sym).as_dict()
            if(len(pis_pol) == 1):
                if pis_pol.has_key((0,)):
                    res +=[[0, pi]]
                elif pis_pol.has_key((1,)):
                    res +=[[nsimplify(pis_pol[(1,)]), 0]]
            elif len(pis_pol) > 1:
                res +=[[nsimplify(pis_pol[(1,)]), nsimplify(pis_pol[(0,)])]]
        return res
        
    
    def __repr__(self):
        res  = "Degree : "        + str(self.deg) 
        res += "\nEquation : "    + self.eq.__repr__()
        res += "\nTime mult : "   + str(self.mult) 
        res += "\nwaypoints :\n"  + str(self.pis) 
        return res
        
        
    def __str__(self):
        return self.__repr__()

  
    
def coeff_at_j(j,n,pis_pol):
    return simplify(sum([ simplify(bino(j,k) / bino(n,k) * pis_pol[(k,)]) for k in range(j+1) ] ))
    
def getFullDic(p):
    dic = p.as_dict()
    monoms = p.all_monoms()
    for m in monoms:
        if not dic.has_key(m):
            dic[m] = 0
    return dic
    
def polybasisToBernBasis(eq, mult, symb = t):
    # to get monomials  in dict
    pis_pol = getFullDic(eq.as_poly(symb))
    #get degree
    deg = len(pis_pol.keys()) -1
    pis = [nsimplify(coeff_at_j(j,deg,pis_pol),tolerance=1e-10,rational=True) for j in range(deg+1)]    
    return Bezier(deg = deg, pis = pis, mult = mult)

allSymbols = p_is[:] + [x]


  
#normalized De Casteljeau reduction, so 0 <= us <=1
def de_casteljau_reduction(pts, us):
    if len(pts) == 1:
        return pts
    res = []
    for i, pt in enumerate(pts[:-1]):
        pt1 = pts[i+1]
        res += [simplify((1 - us) * pt  + us * pt1)]
    return res
    
     
def split_normalized(b,us):
     pts = b.pis
     wps_first = [pts[0]]
     wps_second = [pts[1]]
     while len(pts) > 1:
         pts = de_casteljau_reduction(pts, us)
         wps_first = wps_first + [pts[0]]
         wps_second = [pts[-1]] + wps_second
     deg = b.deg
     mult = b.mult
     return Bezier(deg, wps_first, mult = mult), Bezier(deg, wps_second, mult = mult)

def crossElement(s0,e0,s1,e1):
    if s0 == s1:
        return 0
    elif s0 == x:
        return simplify((-e1)*e0*p_isX[allSymbols.index(s1)]*s0)
    else:
        # here, take order of apparition into account to allow 
        # future factorization
        #~ # x is thus always last which is what we want in the end...
        id0 = allSymbols.index(s0)
        id1 = allSymbols.index(s1)
        if(id0 < id1):
            return simplify((e1)*e0*p_isX[id0]*s1)
        else:
            return simplify((-e1)*e0*p_isX[id1]*s0)
        
def getAllSymbols(eq):
    res = []
    for el in allSymbols:
        dic = eq.as_poly(el).as_dict()
        if(len(dic) > 1): 
            res += [(el, dic[(1,)])]
    return res;

def simplifycross(eq1,eq2, symb = t):
    decomp1 = getAllSymbols(eq1)
    decomp2 = getAllSymbols(eq2)
    res = 0
    for (s0, e0) in decomp1:
        for (s1, e1) in decomp2:
            res += crossElement(s0,e0,s1,e1)
    return res
    
def crossBezier(b1,b2, symb = t):
    neq = simplifycross(b1.eq, b2.eq, symb)
    return polybasisToBernBasis(neq, b1.mult * b2.mult, symb = symb)

#function generation
off = 0


def gXC(c_of_t, symb = t):
    return polybasisToBernBasis(simplify(p_isX[0] * c_of_t.eq),  c_of_t.mult, symb = symb)
    
def wdown(c_of_t):
    ddc = c_of_t.diff().diff()
    return  polybasisToBernBasis(simplify(crossBezier(c_of_t,ddc).eq + gXC(c_of_t).eq), ddc.mult)
    
def wup(ddc, wd):
    res = ddc
    while(res.deg < wd.deg):
        res = res.elevate()
    return res

def computewUpDownWaypoints(c_of_t, sym = x):
    wd = wdown(c_of_t)
    wu = wup(c_of_t.diff().diff(), wd)
    #~ return (wd.wps_as_factor(sym),wu.wps_as_factor(sym))
    return (wd,wu)

def wps_as_factor(self, sym = x):
        res = []
        for _, pi in enumerate(self.pis):
            pis_pol = pi.as_poly(sym).as_dict()
            if(len(pis_pol) == 1):
                if pis_pol.has_key((0,)):
                    res +=[[0, pi]]
                elif pis_pol.has_key((1,)):
                    res +=[[nsimplify(pis_pol[(1,)]), 0]]
            elif len(pis_pol) > 1:
                res +=[[nsimplify(pis_pol[(1,)]), nsimplify(pis_pol[(0,)])]]
        return res
        
# tests

if __name__ == '__main__':

    b3 = Bezier(3)

    #~ a = getBezierFromConstraints()
    #~ bc = getBezierFromConstraints(init_pos = True, init_vel = True, init_acc = False, end_pos  = False,  end_vel = True,  end_acc = False)
    a = getBezierFromConstraints(init_pos = True, init_vel = True, init_acc = True, end_pos  = False,  end_vel = True,  end_acc = True)
    b = a.diff()
    c = b.diff()

    wp = b.pis
    d = crossBezier(a,c)

    gc = gXC(a)

    (wd,wu) =  computewUpDownWaypoints(a)

    a1 = getBezierFromConstraints(init_pos = True, init_vel = True, init_acc = True, end_pos  = False,  end_vel = False,  end_acc = False)
    b1 = a1.diff()
    c1 = b1.diff()

    wp1 = b1.pis
    d1 = crossBezier(a1,c1)
        
    gc = gXC(a1)

    (wd1,wu1) =  computewUpDownWaypoints(a1, y)

    midpis = [0.5* a.pis[i] + 0.5 * a1.pis[i] for i in range(4)]
    amid = Bezier(deg = a.deg, pis = midpis)

    (wdmid,wumid) =  computewUpDownWaypoints(amid)


    i = 0
    for i in range (4):
        print simplify((wd.pis[i] + wd1.pis[i] - 2. * wdmid.pis[i])/3*T**2)

