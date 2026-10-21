#!/usr/bin/python3


import os
import numpy as np

def unquote(s):
    return s[1+s.find('"'):s.rfind('"')]

def uncomment(s):
    return s[2+s.find('/*'):s.rfind('*/')]

def col(c):
    color = c.split('/*')
    value = unquote(color[1])
    color = unquote(color[0]).split()
    return color[0], value

def load_distance_matrix(matrix_file):
    xpm = open(matrix_file)
    meta = [xpm.readline()]
    while not meta[-1].startswith("static char *gromacs_xpm[]"):
        meta.append(xpm.readline())
    dim = xpm.readline()
    nx, ny, nc, nb = [int(i) for i in unquote(dim).split()]
    colors = dict([col(xpm.readline()) for i in range(nc)])
    mdmat = []
    for i in xpm:
        if i.startswith("/*"):
            continue
        j = unquote(i)
        z = [float(colors[j[k:k+nb]]) for k in range(0,nx,nb)]
        mdmat.append(z)
    mdmat = np.array(mdmat)
    return mdmat

def smallest_average_row(matrix_file):
    mdmat = load_distance_matrix(matrix_file)
    row_averages = np.mean(mdmat, axis=1)
    smallest_average = np.argmin(row_averages)
    return smallest_average

# Example usage:

print(smallest_average_row("matrix.xpm"))
