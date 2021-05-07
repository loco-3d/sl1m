import numpy as np

start = [[1.5, 1., 0.], [-2.5, 1., 0.], [-2.5, -1., 0.], [1.5, -1., 0.]]
end = [[2., 1., 0.1], [4., 1., 0.1], [4., -1., 0.1], [2., -1., 0.1]]

astart = np.array(start).T
aend = np.array(end).T

scene = [[astart], [aend]]
surfaces = [[[astart], [astart]], [[astart, aend], [astart, aend]], [[astart, aend], [astart, aend]], [
    [astart, aend], [astart, aend]], [[astart, aend], [astart, aend]], [[astart, aend], [astart, aend]], [[aend], [aend]]]
surfaces = [[[astart, aend], [astart, aend]]]
walk_surfaces = [[[astart]], [[astart]], [[astart, aend]], [[astart, aend]], [[astart, aend]], [[astart, aend]], [
    [astart, aend]], [[astart, aend]], [[astart, aend]], [[astart, aend]], [[astart, aend]], [[astart, aend]], [[aend]], [[aend]]]