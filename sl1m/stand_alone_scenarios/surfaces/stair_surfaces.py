import numpy as np

floor = [ [0.16, 1., 0.], [-1.8, 1., 0.], [-1.8, -1., 0.], [0.16, -1., 0.]  ]

step1 = [[0.3, 0.6, 0.1],[0.3, -0.16, 0.1],[0.6, -0.16, 0.1],[0.6, 0.6, 0.1]]
step2 = [[0.6, 0.6, 0.2 ],[0.6, -0.16, 0.2],[0.9, -0.16, 0.2 ],[0.9, 0.6, 0.2 ]]
step3 = [[0.9, 0.6, 0.3],[0.9, -0.16, 0.3],[1.2, -0.16, 0.3],[1.2, 0.6, 0.3]]
step4 = [[1.2, 0.6, 0.4 ],[1.2, -0.16, 0.4 ],[1.5, -0.16, 0.4 ],[1.5, 0.6, 0.4 ]]

afloor = np.array(floor).T
astep1 = np.array(step1).T
astep2 = np.array(step2).T
astep3 = np.array(step3).T
astep4 = np.array(step4).T

surfaces = [[afloor], [afloor], [astep1,astep2,astep3],[astep2,astep3,astep1], [astep3,astep2,astep1,astep4], [astep3,astep4], [astep4],[astep4]]