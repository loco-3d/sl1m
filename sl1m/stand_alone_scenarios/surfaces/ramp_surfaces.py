import numpy as np

begin = np.array([[1.75, 1.65, 1.65, 1.75],
                  [0.3, 0.3, 0.1, 0.1],
                  [0.6, 0.6, 0.6, 0.6]])

platform = np.array([[2.5, 1.5, 1.5, 2.5],
                     [0.9, 0.9, -1.1, -1.1],
                     [0.6, 0.6, 0.6, 0.6]])

bridge = np.array([[-1.5, -1.5, 1.5, 1.5],
                   [-0.5, -0.8, -0.8, -0.5],
                   [0.6, 0.6, 0.6, 0.6]])

end = np.array([[-1.5, -2.0, -2.0, -1.5],
                [-0.4, -0.4, -1.1, -1.1],
                [0.6, 0.6, 0.6, 0.6]])

surfaces = []
surfaces += [[begin]]+[[begin]]
for i in range(20):
    surfaces += [[platform, bridge]]
surfaces += [[end]]+[[end]]