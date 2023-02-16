import numpy as np

floor = [[0.5, -0.74, 0.], [0.5, 0.74, 0.], [-0.5, 0.74, 0.], [-0.5, -0.74, 0.]]

step1 = [[1.72, -0.48, -0.15], [1.72, 0.3, -0.15], [1.16, 0.3, -0.15], [1.16, -0.48, -0.15]]
step2 = [[0.66, 0.45, -0.126], [0.66, 0.07, -0.126], [0.98, 0.07, -0.011], [0.98, 0.45, -0.01]]
step3 = [[0.8, -0.03, -0.01], [0.58, -0.24, -0.126], [0.84, -0.51, -0.126], [1.06, -0.3, -0.01]]

afloor = np.array(floor).T
astep1 = np.array(step1).T
astep2 = np.array(step2).T
astep3 = np.array(step3).T
scene = [[afloor], [astep1], [astep2], [astep3]]

quadruped_surfaces = [
    [afloor],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
    [afloor, astep1, astep2, astep3],
]
quadruped_surfaces_gait = [[surface] for surface in quadruped_surfaces]