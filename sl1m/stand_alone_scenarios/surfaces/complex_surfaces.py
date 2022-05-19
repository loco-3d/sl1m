import numpy as np

start = [
    [-3.0, 0.4, 0.0],
    [-2.7, 0.4, 0.0],
    [-2.7, 0.1, 0.0],
    [-3.0, 0.1, 0.0],
]
floor = [
    [-0.30, 0.54, 0.0],
    [0.01, 0.54, 0.0],
    [0.01, -0.46, 0.0],
    [-0.30, -0.46, 0.0],
]
step1 = [
    [0.01, 0.54, 0.1],
    [0.31, 0.54, 0.1],
    [0.31, -0.46, 0.1],
    [0.01, -0.46, 0.1],
]
step2 = [
    [0.31, 0.54, 0.2],
    [0.61, 0.54, 0.2],
    [0.61, -0.46, 0.2],
    [0.31, -0.46, 0.2],
]
step3 = [
    [0.61, 0.54, 0.3],
    [0.91, 0.54, 0.3],
    [0.91, -0.46, 0.3],
    [0.61, -0.46, 0.3],
]
step4 = [
    [0.91, 0.54, 0.4],
    [1.21, 0.54, 0.4],
    [1.21, -0.46, 0.4],
    [0.91, -0.46, 0.4],
]
step5 = [
    [1.21, 0.54, 0.5],
    [1.51, 0.54, 0.5],
    [1.51, -0.46, 0.5],
    [1.21, -0.46, 0.5],
]
step6 = [
    [1.51, 0.54, 0.6],
    [1.81, 0.54, 0.6],
    [1.81, -0.46, 0.6],
    [1.51, -0.46, 0.6],
]
step7 = [
    [1.51, -0.46, 0.6],
    [1.81, -0.46, 0.6],
    [1.81, -0.76, 0.6],
    [1.51, -0.76, 0.6],
]
bridge = [
    [1.51, -0.46, 0.6],
    [1.51, -0.76, 0.6],
    [-1.49, -0.76, 0.6],
    [-1.49, -0.46, 0.6],
]
platfo = [
    [-1.49, -0.35, 0.6],
    [-1.49, -1.06, 0.6],
    [-2.49, -1.06, 0.6],
    [-2.49, -0.35, 0.6],
]
slope = [
    [-1.49, -0.06, 0.6],
    [-1.49, 1.5, 0.0],
    [-2.49, 1.5, 0.0],
    [-2.49, -0.06, 0.6],
]
rub1 = [
    [-3.0, -0.15, 0.0],
    [-2.45, -0.15, 0.0],
    [-2.45, 0.53, 0.0],
    [-3.0, 0.53, 0.0],
]
rub2 = [
    [-2.11, 0.30, 0.05],
    [-2.45, 0.30, 0.05],
    [-2.45, 0.53, 0.05],
    [-2.11, 0.53, 0.05],
]
rub3 = [
    [-1.91, -0.15, 0.1],
    [-2.25, -0.15, 0.1],
    [-2.25, 0.18, 0.1],
    [-1.91, 0.18, 0.1],
]
rub4 = [
    [-1.69, 0.19, 0.15],
    [-2.03, 0.19, 0.15],
    [-2.03, 0.53, 0.15],
    [-1.69, 0.53, 0.15],
]
rub5 = [
    [-1.49, -0.15, 0.2],
    [-1.83, -0.15, 0.2],
    [-1.83, 0.18, 0.2],
    [-1.49, 0.18, 0.2],
]
rub6 = [
    [-1.29, 0.19, 0.2],
    [-1.63, 0.19, 0.2],
    [-1.63, 0.53, 0.2],
    [-1.29, 0.53, 0.2],
]
rub7 = [
    [-1.09, -0.15, 0.15],
    [-1.43, -0.15, 0.15],
    [-1.43, 0.18, 0.15],
    [-1.09, 0.18, 0.15],
]
rub75 = [
    [-0.89, 0.19, 0.1],
    [-1.23, 0.19, 0.1],
    [-1.23, 0.53, 0.1],
    [-0.89, 0.53, 0.1],
]
rub8 = [
    [-0.89, -0.15, 0.025],
    [-1.02, -0.15, 0.025],
    [-1.02, 0.18, 0.025],
    [-0.89, 0.18, 0.025],
]
rub9 = [
    [-0.35, -0.15, 0.025],
    [-0.86, -0.15, 0.025],
    [-0.86, 0.52, 0.025],
    [-0.35, 0.52, 0.025],
]
rub8 = [
    [-0.89, -0.15, 0.05],
    [-1.02, -0.15, 0.05],
    [-1.02, 0.18, 0.05],
    [-0.89, 0.18, 0.05],
]
rub9 = [
    [-0.35, -0.15, 0.05],
    [-0.86, -0.15, 0.05],
    [-0.86, 0.52, 0.05],
    [-0.35, 0.52, 0.05],
]

arub9 = np.array(rub9).T
arub8 = np.array(rub8).T
arub75 = np.array(rub75).T
arub7 = np.array(rub7).T
arub6 = np.array(rub6).T
arub5 = np.array(rub5).T
arub4 = np.array(rub4).T
arub3 = np.array(rub3).T
arub2 = np.array(rub2).T
arub1 = np.array(rub1).T
astart = np.array(start).T
afloor = np.array(floor).T
astep1 = np.array(step1).T
astep2 = np.array(step2).T
astep3 = np.array(step3).T
astep4 = np.array(step4).T
astep5 = np.array(step5).T
astep6 = np.array(step6).T
astep7 = np.array(step7).T
abridge = np.array(bridge).T
aplatfo = np.array(platfo).T
aslope = np.array(slope).T

scene = [
    [astart],
    [afloor],
    [astep1],
    [astep2],
    [astep3],
    [astep4],
    [astep5],
    [astep6],
    [astep7],
    [abridge],
    [aplatfo],
    [arub8],
    [arub9],
    [arub7],
    [arub75],
    [arub6],
    [arub5],
    [arub4],
    [arub3],
    [arub2],
    [arub1],
]

allrub = [arub2, arub3, arub5, arub4, arub6, arub7, arub75, arub9]
allsteps = [astep2, astep1, astep3, astep4, astep5, astep6, astep7]
allstepsandfloor = allsteps + [arub9, afloor]
allrubfloorsteps = allrub + allsteps + [afloor]
end = [astep6, astep7, abridge, aplatfo]

surfaces = [
    [arub1, arub2],
    [arub1, arub2],
    [arub1, arub2],
    [arub1, arub3, arub2],
    [arub4],
    [arub5],
    [arub6],
    [arub7],
    [arub75],
    allrub,
    [arub7, arub8, arub9],
    [arub7, arub8, arub9],
    [afloor],
    [afloor, astep1],
    [afloor, astep1],
    [astep1, astep2, astep3],
    [astep4, astep2, astep3],
    [astep4, astep2, astep3],
    [astep4, astep2, astep5],
    [astep6, astep2, astep5],
    [astep6],
    [astep6],
    [astep6],
    [astep6],
    [astep6],
    [astep6],
    [astep6],
    [astep6],
    [astep6],
    [astep6],
    [astep6],
    [astep6, astep7],
    [astep7],
    [abridge],
    [abridge],
    [abridge],
    [abridge],
    [abridge],
    [abridge],
    [abridge],
    [abridge],
    [abridge],
    [abridge],
    [aplatfo],
]
surfaces_gait = [[surface] for surface in surfaces]

rubble_stairs = [
    [arub1, arub2, arub3],
    [arub1, arub2, arub3],
    [arub1, arub2, arub3],
    [arub3, arub2],
    [arub4, arub3, arub5],
    [arub5, arub4, arub3, arub6],
    [arub6],
    [arub7],
    [arub75],
    [arub9, afloor],
    [arub9, afloor],
    [afloor, arub9],
    [astep1],
    [astep2],
    [astep3],
    [astep4],
    [astep5],
    [astep6],
    [astep6],
]
rubble_stairs_gait = [[surface] for surface in rubble_stairs]
