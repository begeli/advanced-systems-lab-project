#!/usr/bin/python

import sys, os
import math

EPS = 1e-9

#--------------------------------------------------
# Main
#--------------------------------------------------

def main():
    gen = ShapePointsGenerator()

    gen.set_shapes(test_name='ProceduralOverlap', o1=procedural(), o2=procedural())
    gen.store_shapes_to_csv()

    gen.set_shapes(test_name='ProceduralNonOverlap', o1=procedural(), o2=procedural(z_offset=EPS))
    gen.store_shapes_to_csv()

    gen.set_shapes(test_name='ProceduralSectionOverlap', o1=procedural(), o2=procedural(x_offset=0.5))
    gen.store_shapes_to_csv()

    gen.set_shapes(test_name='ProceduralAwayNonOverlap', o1=procedural(), o2=procedural(x_offset=(2.0 + 2.0*EPS)))
    gen.store_shapes_to_csv()

    # gen.set_shapes(test_name='IcosahedronCongruentOverlap', o1=icosahedron(), o2=icosahedron())
    # gen.store_shapes_to_csv()

#--------------------------------------------------
# LOGIC
#--------------------------------------------------

class ShapePointsGenerator():
    def __init__(self, path_test_csv_dirs='./test_csv_dirs/'):
        self.set_path_test_csv_dirs(path_test_csv_dirs)

    def set_shapes(self, o1, o2, test_name):
        self.test_name = test_name
        self.shapes = []
        self.shapes.append(o1)
        self.shapes.append(o2)
        self.test_name = test_name
        return self

    def set_path_test_csv_dirs(self, path_test_csv_dirs):
        assert len(path_test_csv_dirs)>2
        self.path_test_csv_dirs = path_test_csv_dirs if path_test_csv_dirs[-1] == '/' else path_test_csv_dirs + '/'
        return self

    def store_shapes_to_csv(self):
        assert len(self.shapes)==2

        old_tests = sorted(os.listdir(self.path_test_csv_dirs))
        highest_test_index = int(old_tests[-1].split('_')[0])
        current_test_index = str(max(highest_test_index, len(old_tests)) + 1)
        dirbase = f'_{self.test_name}'
        dirname = f'{self.path_test_csv_dirs}{current_test_index}{dirbase}'

        for t in old_tests:
            if dirbase in t:
                print(f'Skipping existing test: {dirname}')
                return

        os.mkdir(dirname)

        for i in [1, 2]:
            fn = f'{dirname}/o{i}.csv'
            print(f'Writing {fn}')
            points = []
            for point in self.shapes[i-1]:
                coords = ''
                for j, coord in enumerate(point):
                    coords = coords + str(coord)
                    if j < len(point)-1:
                        coords = coords + ', '
                points.append(coords)
                points.append('\n')
            points = points[:-1] # remove last newline

            o = ''
            for p in points:
                o = o + p

            with open(fn, 'w') as f:
                f.write(o)

def error(*args, err_code=1):
    print('Error:',*args)
    sys.exit(err_code)

#--------------------------------------------------
# SHAPES
#--------------------------------------------------
def icosahedron():
    """
    Returns list of point coordinates of an icosahedron.

    How to build an ico sphere was already solved for blender3D:
    http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html

    Basically the points of 3 orthogonal centered rectangles
    """
    t = (1.0 + 5**(0.5)) * 0.5;
    ico = [[-1,  t, 0],
            [ 1,  t, 0],
            [ 1, -t, 0],
            [-1, -t, 0],

            [ 0, -1, t],
            [ 0,  1, t],
            [ 0,  1,-t],
            [ 0, -1,-t],

            [ t, 0, -1],
            [ t, 0,  1],
            [-t, 0,  1],
            [-t, 0, -1]]
    return ico


def procedural(count=100, x_offset=0.0, y_offset=0.0, z_offset=0.0):
    o = []
    for i in range(count):
        radian = i * 2.0 * math.pi / float(count);
        p = [
                math.cos(radian) + x_offset,
                math.sin(radian) + y_offset,
                z_offset
            ]
        o.append(p)
    return o


if __name__ == "__main__":
    main()
