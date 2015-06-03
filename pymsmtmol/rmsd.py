#------------------------------------------------------------------------------
""" Calculate RMSD between two XYZ files
by: Jimmy Charnley Kromann <jimmy@charnley.dk> and Lars Bratholm
project: https://github.com/charnley/rmsd
license: https://github.com/charnley/rmsd/blob/master/LICENSE
"""

import numpy
import sys
import re

def fit(P, Q):
    """ Varies the distance between P and Q, and optimizes rotation for each
    step until a minimum is found."""

    step_size = P.max(0)
    threshold = step_size*1e-9
    rmsd_best = kabsch(P, Q)
    while True:
        for i in range(3):
            temp = numpy.zeros(3)
            temp[i] = step_size[i]
            rmsd_new = kabsch(P+temp, Q)
            if rmsd_new < rmsd_best:
                rmsd_best = rmsd_new
                P[:,i] += step_size[i]
            else:
                rmsd_new = kabsch(P-temp, Q)
                if rmsd_new < rmsd_best:
                    rmsd_best = rmsd_new
                    P[:,i] -= step_size[i]
                else:
                    step_size[i] /= 2
        if (step_size<threshold).all():
            break
    return rmsd_best


def kabsch(P, Q, output=False):
    """ The Kabsch algorithm
    http://en.wikipedia.org/wiki/Kabsch_algorithm
    The algorithm starts with two sets of paired points P and Q.
    P and Q should already be centered on top of each other.
    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.
    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U
    The optimal rotation matrix U is then used to
    rotate P unto Q so the RMSD can be caculated
    from a straight forward fashion.
    """

    # Computation of the covariance matrix
    C = numpy.dot(numpy.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = numpy.linalg.svd(C)
    d = (numpy.linalg.det(V) * numpy.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    # Create Rotation matrix U
    U = numpy.dot(V, W)

    # Rotate P
    P = numpy.dot(P, U)

    if output:
        return P, rmsd(P, Q)

    return rmsd(P, Q)


def centroid(X):
    """ Calculate the centroid from a vectorset X """
    C = sum(X)/len(X)
    return C


def rmsd(V, W):
    """ Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(W)

    rmsd = 0.0

    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])

    return numpy.sqrt(rmsd/N)

def cal_rmsd(P, Q):
    """
    Calculate Root-mean-square deviation (RMSD) between two molecules.
    The two sets of atoms must be in the same order.
    The script will return three RMSD values;
    1) Normal: The RMSD calculated the straight-forward way.
    2) Kabsch: The RMSD after the two coordinate sets are translated and rotated onto eachother.
    3) Fitted: The RMSD after a fitting function has optimized the centers of the two coordinat sets.
    """
    # Calculate 'dumb' RMSD
    #val1 = round(rmsd(P, Q), 3)

    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    # http://en.wikipedia.org/wiki/Centroid
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc

    #val2 = round(kabsch(P, Q), 3)
    val3 = round(fit(P, Q), 3)

    return val3

    #print "Normal RMSD:", normal_rmsd
    #print "Kabsch RMSD:", round(kabsch(P, Q), 3), "Fitted RMSD:", round(fit(P, Q), 3)
#------------------------------------------------------------------------------
