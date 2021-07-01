#!/usr/bin/python

#--------------------------------------------------
# visualize GJK test case shapes
#--------------------------------------------------

import sys, os
import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from mpl_toolkits.mplot3d import Axes3D

def main():

    #--------------------------------------------------
    # did you say dynamic svn?
    #--------------------------------------------------
    interactive = False
    test_indices = []
    if len(sys.argv) == 1:
        print('\n\n\nNOTE:\nFor interactive svn, pass argument -i.\nIf one or more numbers are passed, only those test cases will be plotted.\n\n')
    else:
        if '-i' in sys.argv:
            interactive = True
        for a in sys.argv:
            if a.isdigit():
                test_indices.append(int(a))

    tp = './test_csv_dirs/'
    tests = sorted(os.listdir(tp))
    if not test_indices == []:
        tmp = []
        for i in test_indices:
            tmp.append(tests[i-1])
        tests = tmp

    #--------------------------------------------------
    # visualize
    #--------------------------------------------------
    for t in tests:
        title = visualize_test(tp, t)
        if interactive:
            plt.show()
        else:
            store_plot(title)
        plt.close()


def visualize_test(test_path, test_file_name):
    fn = test_path + test_file_name
    o1 = np.array(pd.read_csv(fn+'/o1.csv',header=None,names=['x','y','z']))
    o2 = np.array(pd.read_csv(fn+'/o2.csv',header=None,names=['x','y','z']))

    title = f'{test_file_name}_EPS_case' if '10' in fn or '12' in fn else test_file_name
    ax = plt.subplot(projection='3d')
    ax.set_title(title.replace('_',' '))
    for i in ["x", "y", "z"]:
        eval(f'ax.set_{i}label("{i}")')

    #-------------------------------------------------------------------
    # draw lines, trying to leverage 3D complex hull algorithm (qhull),
    # or at least handcraft some lines to connect points of same shape.
    #-------------------------------------------------------------------
    try:
        hull3d(o1,o2,ax)
    except:
        print(f'Cannot construct 3D simplex. Object has only {min(len(o1), len(o2))} corners and is probably only 2D.', end=' ')
        handcraft(o1,o2,ax)

    # also draw just the points.
    for i in [1, 2]:
        lo = 'k+' if i==1 else 'rx'
        eval(f'ax.plot(o{i}[:,0], o{i}[:,1], o{i}[:,2], lo , label="shape {i}")')

    plt.tight_layout()
    return title

def hull3d(o1,o2,ax):
    h1 = ConvexHull(o1)
    h2 = ConvexHull(o2)

    # need to cycle back to first element to get all simplices plotted:
    # https://stackoverflow.com/questions/27270477/3d-convex-hull-from-point-cloud
    for s in h1.simplices:
        s = np.append(s, s[0]); ax.plot(o1[s, 0], o1[s, 1], o1[s,2], 'k:')
    for s in h2.simplices:
        s = np.append(s, s[0]); ax.plot(o2[s, 0], o2[s, 1], o2[s,2], 'r--')
    make_legend(isHandCrafted=False)

def handcraft(o1,o2,ax):
    """
    Stupidly draws lines from point to point in the order they given in the CSV file.
    This is more for visually connecting points from the same shape, when the shapes are too degenerate
    for a 3D convex hull algorithm.
    """
    o1 = o1 if len(o1) < 5 else np.append(o1, [o1[0]], axis=0)
    o2 = o2 if len(o2) < 5 else np.append(o2, [o2[0]], axis=0)

    for i in [1, 2]:
        li = 'k:' if i==1 else 'r--'
        if eval(f'len(o{i}) == 4'):
            eval(f'handcraft_4(o{i}, li, "shape {i}", ax)')
        else:
            eval(f'ax.plot(o{i}[:,0], o{i}[:,1], o{i}[:,2], li , label="shape {i}")')
    make_legend(isHandCrafted=True)

def handcraft_4(o, linestyle, label, ax):
    lines = [[0,1], [1,3], [2,3], [2,0]]
    for l in lines:
        ax.plot(o[l,0], o[l,1], o[l,2], linestyle, label=label)

    # ax.plot(o[l2,0], o[l2,1], o[l2,2], linestyle, label=label)
    # ax.plot(o[l3,0], o[l3,1], o[l3,2], linestyle, label=label)
    # ax.plot(o[l4,0], o[l4,1], o[l4,2], linestyle, label=label)


def make_legend(isHandCrafted):
    s1 = 'shape 1' if not isHandCrafted else 'shape 1 (generated lines may not cover complete hull)'
    s2 = 'shape 2' if not isHandCrafted else 'shape 2 (generated lines may not cover complete hull)'
    black = mlines.Line2D([], [], color='black', linestyle=':',  marker='+', label=s1)
    red = mlines.Line2D([], [], color='red', linestyle='--',  marker='x', label=s2)
    plt.legend(handles=[black,red])

def store_plot(title):
    d = './svn'
    if not os.path.isdir(d):
        os.mkdir(d)
    p = f'{d}/{title}'
    print(f'storing {p}.png...')
    plt.savefig(p)

if __name__ == "__main__":
    main()

