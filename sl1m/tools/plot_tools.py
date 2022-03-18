import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import numpy as np


COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
          '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']


def plot_point(ax, point, color="b", linewidth=2):
    """
    Plot a point
    """
    x = np.array(point)[0]
    y = np.array(point)[1]
    z = np.array(point)[2]
    ax.scatter(x, y, z, color=color, marker='o', linewidth=linewidth)


def plot_surface(points, ax, color_id=0, alpha=1.):
    """
    Plot a surface
    """
    xs = np.append(points[0, :], points[0, 0]).tolist()
    ys = np.append(points[1, :], points[1, 0]).tolist()
    zs = np.append(points[2, :], points[2, 0]).tolist()
    if color_id == -1:
        ax.plot(xs, ys, zs)
    else:
        ax.plot(xs, ys, zs, color=COLORS[color_id % len(COLORS)], alpha=alpha)


def draw_potential_surfaces(surfaces, gait, phase, ax=None, alpha=1.):
    """
    Plot all the potential surfaces of one phase of the problem
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    for surface in surfaces:
        plot_surface(surface, ax, gait[phase % (len(gait))])
    return ax

def draw_potential_surfaces_gait(surfaces, phase, foot_index, ax=None, title=None):
    """
    Plot all the potential surfaces of one phase of the problem
    """
    if ax is None:
        fig = plt.figure()
        if title is not None:
            fig.suptitle(title, fontsize=16)
        ax = fig.add_subplot(111, projection="3d")
    for surface in surfaces[phase][foot_index]:
        plot_surface(surface, ax, foot_index)
    return ax


def draw_whole_scene(surface_dict, ax=None, title=None):
    """
    Plot all the potential surfaces
    """
    if ax is None:
        fig = plt.figure()
        if title is not None:
            fig.suptitle(title, fontsize=16)
        ax = fig.add_subplot(111, projection="3d")
    for key in surface_dict.keys():
        plot_surface(np.array(surface_dict[key][0]).T, ax, 5)
    return ax


def draw_scene(surfaces, gait=False, ax=None, alpha=1.):
    """
    Plot all the potential surfaces of the problem
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    for surfaces_phase in surfaces:
        for surface in surfaces_phase:
            if gait:
                for foot_surface in surface:
                    plot_surface(foot_surface, ax, 0, alpha=alpha)
            else:
                plot_surface(surface, ax, 0, alpha=alpha)
    return ax


def draw_surface(surfaces, foot, ax=None):
    """
    Plot all the potential surfaces
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    for surface in surfaces[0]:
        plot_surface(surface, ax, foot)
    return ax


def plot_initial_contacts(initial_feet_pose, ax=None):
    """
    Plot the initial feet positions
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    ax.grid(False)

    for i, pose in enumerate(initial_feet_pose):
        if pose is not None:
            plot_point(ax, pose, color=COLORS[i])


def plot_new_contact(moving_feet, moving_feet_pos, ax=None):
    """
    Plot the initial feet positions
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    ax.grid(False)

    for i in range(len(moving_feet)):
        plot_point(ax, moving_feet_pos[i], color=COLORS[moving_feet[i]])


def plot_first_step(configs, coms, moving_foot_pos, gait, ax=None):
    """
    Plot the moving feet positions
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    ax.grid(False)

    plot_point(ax, configs[0], color=COLORS[len(gait) + 1])
    plot_point(ax, coms[0], color=COLORS[len(gait) + 2])
    plot_point(ax, moving_foot_pos[0], color=COLORS[gait[0]])


def plot_selected_surfaces(surfaces, surface_indices, gait, ax=None):
    """
    Plot the surface with the minimum alpha
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    for i, surfaces_phase in enumerate(surfaces):
        plot_surface(surfaces_phase[surface_indices[i]], ax, gait[i % len(gait)])
    return ax


def plot_heightmap(heightmap, alpha=1., ax=None):
    """
    Plot the heightmap
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    i = 0
    if alpha != 1.:
        i = 1

    xv, yv = np.meshgrid(heightmap.x, heightmap.y, sparse=False, indexing='ij')
    ax.plot_surface(xv, yv, heightmap.z, color=COLORS[i], alpha=alpha)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_zlim([np.min(heightmap.z), np.max(heightmap.z) + 1.])

    return ax


def plot_point_list(ax, wps, color="b", D3=True, linewidth=2):
    """
    Plot a list of points
    """
    x = np.array(wps)[:, 0]
    y = np.array(wps)[:, 1]
    if(D3):
        z = np.array(wps)[:, 2]
        ax.scatter(x, y, z, c=color, marker='o', linewidth=5)
    else:
        ax.scatter(x, y, color=color, linewidth=linewidth)


def plot_planner_result(all_feet_pos, coms=None, step_size=None, effector_positions=None, shoulder_positions_2D=None, ax=None, show=True):
    """
    Plot the feet positions and com positions
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
    ax.grid(False)
    ax.view_init(elev=8.776933438381377, azim=-99.32358055821186)

    for foot, foot_pose in enumerate(all_feet_pos):
        px = [c[0] for c in foot_pose if c is not None]
        py = [c[1] for c in foot_pose if c is not None]
        pz = [c[2] for c in foot_pose if c is not None]
        ax.scatter(px, py, pz, color=COLORS[foot], marker='o', linewidth=5)

    if step_size is not None:
        for foot_pose in all_feet_pos:
            px = [foot_pose[0][0] + step_size[0]]
            py = [foot_pose[0][1] + step_size[1]]
            pz = [foot_pose[0][2]]
            ax.scatter(px, py, pz, color=COLORS[5], marker='o', linewidth=5)

    if effector_positions is not None:
        for foot, foot_pos in enumerate(effector_positions):
            px = [c[0] for c in foot_pos if c is not None]
            py = [c[1] for c in foot_pos if c is not None]
            pz = [0. for c in foot_pos if c is not None]
            ax.scatter(px, py, pz, color=COLORS[foot], marker='*', linewidth=5)

    if shoulder_positions_2D is not None:
        for foot, foot_pos in enumerate(shoulder_positions_2D):
            px = [c[0] for c in foot_pos if c is not None]
            py = [c[1] for c in foot_pos if c is not None]
            pz = [0. for c in foot_pos if c is not None]
            ax.scatter(px, py, pz, color=COLORS[foot], marker='x', linewidth=5)

    if coms is not None:
        plot_point_list(ax, coms, color=COLORS[len(all_feet_pos)+1])
        cx = [c[0] for c in coms]
        cy = [c[1] for c in coms]
        cz = [c[2] for c in coms]
        ax.plot(cx, cy, cz, color=COLORS[len(all_feet_pos)+1])

    if show:
        plt.draw()
        plt.show()