from cdd import Matrix, Polyhedron, RepType
from numpy import array, hstack, zeros, ones
import numpy as np

NUMBER_TYPE = 'float'  # 'float' or 'fraction'

# CROSSMATRIX Compute the projection matrix of the cross product
def crossMatrix( v ):
    VP = np.array( [[  0,  -v[2], v[1] ],
                    [ v[2],  0,  -v[0] ],
                    [-v[1], v[0],  0   ]] );
    return VP;

class ConeException(Exception):

    def __init__(self, M):
        self.M = M


class NotConeFace(ConeException):

    def __str__(self):
        return "Matrix is not a cone face"


class NotConeSpan(ConeException):

    def __str__(self):
        return "Matrix is not a cone span"

class NotPolyFace(ConeException):

    def __str__(self):
        return "Matrix is not a polytope face"


def cone_span_to_face(S, eliminate_redundancies=False):
    """

    Returns the face matrix S^F of the span matrix S,
    that is, a matrix such that

        {x = S z, z >= 0} if and only if {S^F x <= 0}.

    """
    V = hstack([ones((S.shape[1], 1)), S.T])
    # V-representation: first column is 0 for rays
    V_cdd = Matrix(V, number_type=NUMBER_TYPE)
    V_cdd.rep_type = RepType.GENERATOR
    P = Polyhedron(V_cdd)
    H_matrix = P.get_inequalities();
    if(eliminate_redundancies):
        H_matrix.canonicalize();
    H = array(H_matrix);
    if(H.shape[1]>1):
        b = H[:, 0]
        A = H[:, 1:]
    else:
        b, A = H[:, 0], zeros((H.shape[0],S.shape[0]));
    #~ for i in xrange(H.shape[0]):
        #~ if b[i] != 0:
            #~ raise NotConeSpan(S)
    return -A, b
    
def poly_span_to_face(S):
    """

    Returns the face matrix S^F of the span matrix S,
    that is, a matrix such that

        {x = S z, z >= 0, sum(z)=1} if and only if {S^F x <= s}.

    """
    V = hstack([ones((S.shape[1], 1)), S.T])
    # V-representation: first column is 0 for rays, 1 for vertices
    V_cdd = Matrix(V, number_type=NUMBER_TYPE)
    V_cdd.rep_type = RepType.GENERATOR
    P = Polyhedron(V_cdd)
    H = array(P.get_inequalities()) # H-representation: A x + b >= 0
    b, A = H[:, 0], H[:, 1:]
    return (-A,b)

def arbitrary_span_to_face(S, rv):
    """

    Returns the face matrix S^F of the span matrix S,
    that is, a matrix such that

        {x = S z, z >= 0} if and only if {S^F x <= f}.

    The vector rv specifies whether the corresponding column in S
    is a vertex (1) or a ray (0).
    """
    V = hstack([rv.reshape((S.shape[1], 1)), S.T])
    # V-representation: first column is 0 for rays
    V_cdd = Matrix(V, number_type=NUMBER_TYPE)
    V_cdd.rep_type = RepType.GENERATOR
    P = Polyhedron(V_cdd)
    H = array(P.get_inequalities())
    b, A = H[:, 0], H[:, 1:]
    return (-A,b)

def cone_face_to_span(F):
    """

    Compute the span matrix F^S of the face matrix F,
    that is, a matrix such that

        {F x <= 0} if and only if {x = F^S z, z >= 0}.

    """
    b, A = zeros((F.shape[0], 1)), -F
    # H-representation: A x + b >= 0
    F_cdd = Matrix(hstack([b, A]), number_type=NUMBER_TYPE)
    F_cdd.rep_type = RepType.INEQUALITY
    P = Polyhedron(F_cdd)
    V = array(P.get_generators())
    for i in xrange(V.shape[0]):
        if V[i, 0] != 0:  # 1 = vertex, 0 = ray
            raise NotConeFace(F)
    return V[:, 1:].T


def poly_face_to_span(F,f):
    """

    Compute the span matrix F^S of the face matrix F,
    that is, a matrix such that

        {F x <= f} if and only if {x = F^S z, z >= 0, sum(z)=1}.

    """
    b, A = f.reshape((F.shape[0],1)), -F
    # H-representation: A x + b >= 0
    F_cdd = Matrix(hstack([b, A]), number_type=NUMBER_TYPE)
    F_cdd.rep_type = RepType.INEQUALITY
    P = Polyhedron(F_cdd)
    V = array(P.get_generators())
    for i in xrange(V.shape[0]):
        if V[i, 0] != 1:  # 1 = vertex, 0 = ray
            raise NotPolyFace(F)
    return V[:, 1:].T
    
def arbitrary_face_to_span(F,f):
    """

    Compute the span matrix F^S of the face matrix F,
    that is, a matrix such that

        {F x <= f} if and only if {x = F^S z, z >= 0, sum(z)=1}.
        
    Works for both polytopes and cones.

    """
    b, A = f.reshape((F.shape[0],1)), -F
    # H-representation: A x + b >= 0
    F_cdd = Matrix(hstack([b, A]), number_type=NUMBER_TYPE)
    F_cdd.rep_type = RepType.INEQUALITY
    P = Polyhedron(F_cdd)
    V = array(P.get_generators())
    return (V[:, 1:].T, V[:, 0]);
    
def eliminate_redundant_inequalities(A,b):
    ''' Input format is  A x + b <= 0'''
    # H-representation: A x + b >= 0
    A_plot = Matrix(hstack([-b.reshape((A.shape[0],1)), -A]), number_type=NUMBER_TYPE)
    A_plot.rep_type = RepType.INEQUALITY
    A_plot.canonicalize()
    A_plot = np.array(A_plot)
    b_plot, A_plot = -A_plot[:, 0], -A_plot[:, 1:]
    return (A_plot, b_plot);
    
if __name__ == "__main__":
    print "Gonna test polytope conversion utils"

    print " *** Test cone span to face ***"
    S = np.array([[1,-1, 1],
                  [1, 1, 2]]);
    F = cone_span_to_face(S);
    print "Span"
    print S
    print "Face"
    print F
    
    print "\n *** Test polytope span to face ***"
    S = np.array([[1, 1,-1, -1, 0.8,-0.3],
                  [1,-1, 1, -1,-0.9, 0.1]]);
    (F,f) = poly_span_to_face(S);
    print "Span"
    print S
    print "Face"
    print F
    print f
    
    
