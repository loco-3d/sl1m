# do the loading of the obj file
import numpy as np
from collections import namedtuple

ObjectData = namedtuple("ObjectData", "V T N F")
Inequalities = namedtuple("Inequality", "A b N V")

__EPS = 0.0


def toFloat(stringArray):
    res = np.zeros(len(stringArray))
    for i in range(0, len(stringArray)):
        res[i] = float(stringArray[i])
    return res


def load_obj(filename):
    V = []  # vertex
    T = []  # texcoords
    N = []  # normals
    F = []  # face indexies

    fh = open(filename)
    for line in fh:
        if line[0] == "#":
            continue

        line = line.strip().split(" ")
        if line[0] == "v":  # vertex
            V.append(toFloat(line[1:]))
        elif line[0] == "vt":  # tex-coord
            T.append(line[1:])
        elif line[0] == "vn":  # normal vector
            N.append(toFloat(line[1:]))
        elif line[0] == "f":  # face
            face = line[1:]
            for i in range(0, len(face)):
                face[i] = face[i].split("/")
                # OBJ indexies are 1 based not 0 based hence the -1
                # convert indexies to integer
                for j in range(0, len(face[i])):
                    if j != 1:
                        face[i][j] = int(face[i][j]) - 1
            F.append(face)
    fh.close()
    return ObjectData(V, T, N, F)


# find a point such that ax + by + cz = d that is the closest to the origin
def find_point_on_plane(a, b, c, d):
    # project 0 to the plane
    m = np.zeros((4, 4))
    m[:3, :3] = np.identity(3)
    m[:, 3] = [-a, -b, -c, 0]
    m[3, :3] = [a, b, c]
    n = np.zeros(4)
    n[-1] = d
    res = np.linalg.inv(m).dot(n)[:-1]

    return res


# Create an Inequalities object given the inequalities matrices A,b, s.t. Ax <=b
def inequalities_to_Inequalities_object(A, b):
    nrows = A.shape[0]
    V = np.ones([nrows, 4])
    N = np.empty([nrows, 3])
    i = 0
    for ai, bi in zip(A, b):
        N[i, :] = ai
        V[i, :3] = find_point_on_plane(ai[0], ai[1], ai[2], bi)
        i += 1
    return Inequalities(A, b, N, V)


def inequality(v, n):
    # the plan has for equation ax + by + cz = d, with a b c coordinates of the normal
    # inequality is then ax + by +cz -d <= 0
    # last var is v because we need it
    return [n[0], n[1], n[2], np.array(v).dot(np.array(n)) + __EPS]


def as_inequalities(obj):
    # for each face, find first three points and deduce plane
    # inequality is given by normal
    A = np.empty([len(obj.F), 3])
    b = np.empty(len(obj.F))
    V = np.ones([len(obj.F), 4])
    N = np.empty([len(obj.F), 3])
    for f in range(0, len(obj.F)):
        face = obj.F[f]
        v = obj.V[face[0][0]]
        # assume normals are in obj
        n = obj.N[face[0][2]]
        ineq = inequality(v, n)
        A[f, :] = ineq[0:3]
        b[f] = ineq[3]
        V[f, 0:3] = v
        N[f, :] = n
    return Inequalities(A, b, N, V)


def is_inside(inequalities, pt):
    return ((inequalities.A.dot(pt) - inequalities.b) < 0).all()


# ~ def rotate_inequalities_q():

# TODO this is naive, should be a way to simply update d
def rotate_inequalities(ineq, transform):
    # for each face, find first three points and deduce plane
    # inequality is given by normal
    A = np.empty([len(ineq.A), 3])
    b = np.empty(len(ineq.b))
    V = np.ones([len(ineq.V), 4])
    N = np.ones([len(ineq.N), 3])
    for i in range(0, len(b)):
        v = transform.dot(ineq.V[i, :])
        n = transform[0:3, 0:3].dot(ineq.N[i, 0:3])
        ine = inequality(v[0:3], n[0:3])
        A[i, :] = ine[0:3]
        b[i] = ine[3]
        V[i, :] = v
        N[i, :] = n
    return Inequalities(A, b, N, V)


from pickle import dump


def ineq_to_file(ineq, filename):
    f1 = open(filename, "w+")
    res = {"A": ineq.A, "b": ineq.b, "N": ineq.N, "V": ineq.V}
    dump(res, f1)
    f1.close()


from pickle import load


def ineq_from_file(filename):
    f1 = open(filename, "r")
    res = load(f1)
    return Inequalities(res["A"], res["b"], res["N"], res["V"])
