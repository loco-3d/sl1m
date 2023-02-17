import numpy as np

start = [
    [1.95, 1.0, 0.0],
    [-2.5, 1.0, 0.0],
    [-2.5, -1.0, 0.0],
    [1.95, -1.0, 0.0],
]
end = [[2.0, 1.0, 0.0], [4.0, 1.0, 0.0], [4.0, -1.0, 0.0], [2.0, -1.0, 0.0]]

astart = np.array(start).T
aend = np.array(end).T

scene = [[astart], [aend]]
jumping_trot_surfaces = [
    [[astart], [astart]],
    [],
    [[astart, aend], [astart, aend]],
    [],
    [[astart, aend], [astart, aend]],
    [],
    [[astart, aend], [astart, aend]],
    [],
    [[aend], [aend]],
]
trot_surfaces = [
    [[astart], [astart]],
    [[astart, aend], [astart, aend]],
    [[astart, aend], [astart, aend]],
    [[astart, aend], [astart, aend]],
    [[astart, aend], [astart, aend]],
    [[astart, aend], [astart, aend]],
    [[astart, aend], [astart, aend]],
    [[astart, aend], [astart, aend]],
    [[astart, aend], [astart, aend]],
    [[astart, aend], [astart, aend]],
    [[aend], [aend]],
]
walk_surfaces = [
    [[astart]],
    [[astart]],
    [[astart, aend]],
    [[astart, aend]],
    [[astart, aend]],
    [[astart, aend]],
    [[astart, aend]],
    [[astart, aend]],
    [[astart, aend]],
    [[astart, aend]],
    [[astart, aend]],
    [[astart, aend]],
    [[astart, aend]],
    [[astart, aend]],
    [[astart, aend]],
    [[aend]],
    [[aend]],
]


# trot_surfaces = [[[astart, aend], [astart, aend]], [[astart, aend], [astart, aend]]]
