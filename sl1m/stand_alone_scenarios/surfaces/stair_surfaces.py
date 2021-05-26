import numpy as np

floor = [[0.16, 1., 0.], [-1.8, 1., 0.], [-1.8, -1., 0.], [0.16, -1., 0.]]

step1 = [[0.3, 0.6, 0.1], [0.3, -0.16, 0.1], [0.6, -0.16, 0.1], [0.6, 0.6, 0.1]]
step2 = [[0.6, 0.6, 0.2], [0.6, -0.16, 0.2], [0.9, -0.16, 0.2], [0.9, 0.6, 0.2]]
step3 = [[0.9, 0.6, 0.3], [0.9, -0.16, 0.3], [1.2, -0.16, 0.3], [1.2, 0.6, 0.3]]
step4 = [[1.2, 0.6, 0.4], [1.2, -0.16, 0.4], [1.8, -0.16, 0.4], [1.8, 0.6, 0.4]]

afloor = np.array(floor).T
astep1 = np.array(step1).T
astep2 = np.array(step2).T
astep3 = np.array(step3).T
astep4 = np.array(step4).T

scene = [[afloor], [astep1], [astep2], [astep3], [astep4]]

surfaces = [[afloor], [afloor], [astep1, astep2, astep3], [astep2, astep3, astep1],
            [astep3, astep2, astep1, astep4], [astep3, astep4], [astep4], [astep4]]
surfaces_gait = [[surface] for surface in surfaces]

walk_surfaces = [[[afloor]], [[afloor, astep1, astep2]], [[afloor, astep1, astep2]], [[afloor, astep1, astep2]], [[afloor, astep1, astep2]], [[afloor, astep1, astep2]], [[astep1, astep2, astep3]], [[astep1, astep2, astep3]], [[astep1, astep2, astep3]], [[astep1, astep2, astep3]], [[astep1, astep2, astep3]], [
    [astep1, astep2, astep3]],  [[astep1, astep2, astep3, astep4]], [[astep1, astep2, astep3, astep4]], [[astep1, astep2, astep3, astep4]], [[astep1, astep2, astep3, astep4]], [[astep3, astep4]], [[astep3, astep4]], [[astep3, astep4]], [[astep3, astep4]], [[astep4]], [[astep4]], [[astep4]], [[astep4]]]

quadruped_surfaces = [[afloor],
                      [afloor, astep1, astep2],
                      [afloor, astep1, astep2],
                      [afloor, astep1, astep2],
                      [afloor, astep1, astep2],
                      [afloor, astep1, astep2],
                      [astep1, astep2, astep3],
                      [astep1, astep2, astep3],
                      [astep1, astep2, astep3],
                      [astep1, astep2, astep3],
                      [astep1, astep2, astep3],
                      [astep1, astep2, astep3],
                      [astep1, astep2, astep3, astep4],
                      [astep1, astep2, astep3, astep4],
                      [astep1, astep2, astep3, astep4],
                      [astep1, astep2, astep3, astep4],
                      [astep3, astep4],
                      [astep3, astep4],
                      [astep3, astep4],
                      [astep3, astep4],
                      [astep4],
                      [astep4],
                      [astep4],
                      [astep4]]
quadruped_surfaces_gait = [[surface] for surface in quadruped_surfaces]

floor_small = [[0.29, 1., 0.], [-1.8, 1., 0.], [-1.8, -1., 0.], [0.29, -1., 0.]]
step1_small = [[0.3, 0.6, 0.05], [0.3, -0.16, 0.05], [0.59, -0.16, 0.05], [0.59, 0.6, 0.05]]
step2_small = [[0.6, 0.6, 0.1], [0.6, -0.16, 0.1], [0.89, -0.16, 0.1], [0.89, 0.6, 0.1]]
step3_small = [[0.9, 0.6, 0.15], [0.9, -0.16, 0.15], [1.19, -0.16, 0.15], [1.19, 0.6, 0.15]]
step4_small = [[1.2, 0.6, 0.2], [1.2, -0.16, 0.2], [1.8, -0.16, 0.2], [1.8, 0.6, 0.2]]

afloor_small = np.array(floor_small).T
astep1_small = np.array(step1_small).T
astep2_small = np.array(step2_small).T
astep3_small = np.array(step3_small).T
astep4_small = np.array(step4_small).T

solo_scene = [[afloor_small], [astep1_small], [astep2_small], [astep3_small], [astep4_small]]
solo_surfaces = [[afloor_small],
                 [afloor_small],
                 [afloor_small],
                 [afloor_small],
                 [afloor_small, astep1_small, astep2_small],
                 [afloor_small, astep1_small, astep2_small],
                 [afloor_small, astep1_small, astep2_small],
                 [afloor_small, astep1_small, astep2_small],
                 [afloor_small, astep1_small, astep2_small],
                 [afloor_small, astep1_small, astep2_small],
                 [afloor_small, astep1_small, astep2_small],
                 [afloor_small, astep1_small, astep2_small],
                 [astep1_small, astep2_small, astep3_small],
                 [astep1_small, astep2_small, astep3_small],
                 [astep1_small, astep2_small, astep3_small],
                 [astep1_small, astep2_small, astep3_small],
                 [astep1_small, astep2_small, astep3_small],
                 [astep1_small, astep2_small, astep3_small],
                 [astep1_small, astep2_small, astep3_small],
                 [astep1_small, astep2_small, astep3_small],
                 [astep2_small, astep3_small, astep4_small],
                 [astep2_small, astep3_small, astep4_small],
                 [astep2_small, astep3_small, astep4_small],
                 [astep2_small, astep3_small, astep4_small],
                 [astep2_small, astep3_small, astep4_small],
                 [astep2_small, astep3_small, astep4_small],
                 [astep3_small, astep4_small],
                 [astep3_small, astep4_small],
                 [astep3_small, astep4_small],
                 [astep3_small, astep4_small],
                 [astep4_small],
                 [astep4_small],
                 [astep4_small],
                 [astep4_small]]
solo_surfaces_gait = [[surface] for surface in solo_surfaces]