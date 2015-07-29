#!/usr/bin/env python
import argparse
import matplotlib as mpl

import matplotlib.pylab as plt
import numpy as np
import os
import sys
import time

import ResolvedStellarPops as rsp
from astropy.io import fits
from ResolvedStellarPops.tpagb_path_config import tpagb_path
from ResolvedStellarPops.match.utils import match_diagnostic
    
def move_on(ok, msg='0 to move on: '):
    ok = int(raw_input(msg))
    time.sleep(1)
    return ok

def exclude_tpagb(phot, param):
    from scipy.interpolate import interp1d
    atab = rsp.angst_tables.angst_data
    target, filter1, filter2 = param.upper().split('.')[0].split('_')
    
    lines = open(param, 'r').readlines()
    dmag, dcol, _, colmin, colmax = map(float, lines[4].split()[:-1])
    mag1min, mag1max = map(float, lines[5].split()[:-1])

    maxcol = colmax - dcol
    brightlim = mag1min + dmag
    mag1, mag2 = np.loadtxt(phot, unpack=True)
    # choose color cut
    mincol = find_match_limits(mag1, mag2, color_only=True)[0]
    mtrgb, dmod, av = atab.get_tab5_trgb_av_dmod(target, filters=[filter1,
                                                                  filter2])
    xarr, yarr = put_a_line_on_it(mtrgb, consty=False, filter1=filter1)
    f = interp1d(xarr, yarr)
    faintlim1 = f(mincol)
    faintlim2 = f(maxcol)
    exclude_gate = '1 {0} {2} {0} {3} {1} {3} {1} {4} 0 \n'.format(mincol,
                                                                   maxcol,
                                                                   faintlim1,
                                                                   brightlim,
                                                                   faintlim2)
    lines[7] = exclude_gate
    # not so simple ... need them to be parallelograms.
    # PASS!
    print(exclude_gate)
    # write new param file with exclude/include gate
    os.system('mv {0} {0}_bkup'.format(param))
    with open(param, 'w') as outp:
        [outp.write(l) for l in lines]
    print('wrote %s' % param)
    match_diagnostic(param, phot)
    
def put_a_line_on_it(val, npts=100, consty=True, ax=None,
                     ls='--', annotate=True, filter1=None,
                     annotate_fmt='$TRGB=%.2f$', xmin=-3, xmax=6,
                     ymin=14, ymax=32, ann_kwargs={}, pltkw={}):
    """
    if consty is True: plots a constant y value across ax.xlims().
    if consty is False: plots a constant x on a plot of y vs x-y
    """
    xarr = np.linspace(xmin, xmax, npts)
    yarr = np.linspace(ymin, ymax, npts)
    if consty:
        # just a contsant y value over the plot range of x.
        new_xarr = xarr
    else:
        # a plot of y vs x-y and we want to mark
        # where a constant value of x is
        # e.g, f814w vs f555-f814; val is f555
        new_xarr = val - yarr
        # e.g, f555w vs f555-f814; val is f814
        if filter1 is not None:
            yarr = xarr + val
            new_xarr = xarr

    if ax is not None:
        ax.plot(new_xarr, yarr, ls, **pltkw)
        if annotate:
            xy = (new_xarr[-1] - 0.1, yarr[-1] - 0.2)
            ax.annotate(annotate_fmt % val, xy=xy, ha='right', fontsize=16,
                        **ann_kwargs)
    return new_xarr, yarr


def find_match_limits(mag1, mag2, comp1=90., comp2=90., color_only=False,
                      xlim=None, ylim=None):
    """
    click color limits on a cmd and mag1 mag2 limits on a plot of mag1 vs mag2
    """
    col = mag1 - mag2

    fig, ax = plt.subplots()
    ax.plot(col, mag2, 'o', color='k', ms=3, alpha=0.3, mec='none')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim(ax.get_ylim()[::-1])
    
    if comp1 < 90.:
        ax.hlines(comp2, *ax.get_xlim())

    ok = 1
    while ok == 1:
        print 'click on color min then color max'
        pts = plt.ginput(2, timeout=-1)
        colmin, colmax = [pts[i][0] for i in range(2)]
        ax.vlines(colmin, *ax.get_ylim())
        ax.vlines(colmax, *ax.get_ylim())
        plt.draw()
        ok = move_on(0)

    plt.close()

    inds, = np.nonzero((col < colmax) & (col > colmin))
    data = (colmin, colmax)
    if not color_only:
        fig, ax = plt.subplots()
        ax.plot(mag1, mag2, '.', color='k')
        ok = 1
        while ok == 1:
            print 'click the bright mag value of mag1 and mag2, click a second time to finish'
            pts = plt.ginput(2, timeout=-1)
            mag1max, mag2max = pts[0]
            ax.plot(mag1max, mag2max, 'o', color='r')
            plt.draw()
            ok = move_on(ok)

        plt.close()

        inds, = np.nonzero((mag1 < comp1) & (mag1 > mag1max) &
                           (mag2 < comp2) & (mag2 > mag2max) &
                           (col < colmax) & (col > colmin))

    fig, ax = plt.subplots()
    ax.plot(col, mag2, '.', color='k')
    ax.plot(col[inds], mag2[inds], '.', color='r')
    ax.set_ylim(ax.get_ylim()[::-1])
    if comp2 < 90.:
        ax.hlines(comp2, *ax.get_xlim(), lw=2)
    ax.vlines(colmin, *ax.get_ylim(), lw=2)
    ax.vlines(colmax, *ax.get_ylim(), lw=2)
    if not color_only:
        ax.hlines(mag2max, *ax.get_xlim(), lw=2)
        data = (colmin, colmax, mag1min, mag1max, mag2max)
    
    plt.draw()
    
    print data
    return data


def find_gates(mag1, mag2, param):
    col = mag1 - mag2

    lines = open(param, 'r').readlines()
    colmin, colmax = map(float, lines[4].split()[3:-1])
    mag1min, mag1max = map(float, lines[5].split()[:-1])
    #mag2min, mag2max = map(float, lines[5].split()[:-1])
    # click around
    fig, ax = plt.subplots()
    ax.plot(col, mag2, ',', color='k', alpha=0.2)
    ax.set_ylim(mag1max, mag1min)
    ax.set_xlim(colmin, colmax)

    ok = 1
    while ok == 1:
        print 'click '
        pts = np.asarray(plt.ginput(n=4, timeout=-1))
        exclude_gate = '1 {} 0 \n'.format(' '.join(['%.4f' % p for p in pts.flatten()]))
        pts = np.append(pts, pts[0]).reshape(5,2)
        ax.plot(pts[:,0], pts[:,1], color='r', lw=3, alpha=0.3)
        plt.draw()
        ok = move_on(0)
    lines[7] = exclude_gate
    # not so simple ... need them to be parallelograms.
    # PASS!
    
    # write new param file with exclude/include gate
    os.system('mv {0} {0}_bkup'.format(param))
    with open(param, 'w') as outp:
        [outp.write(l) for l in lines]
    print('wrote %s' % param)
    
    os.chdir(here)



def match_limits(mag1, mag2, color_only=False, comp1=99., comp2=99.):
    plt.ion()
    ok = 1
    while ok == 1:
        data = find_match_limits(mag1, mag2, comp2=comp2, comp1=comp1,
                                 color_only=color_only)
        ok = move_on(0)

    plt.close()

    if color_only:
        colmin, colmax = data
        data_str = '%.2f %.2f' % (colmin, colmax)
    else:
        colmin, colmax, mag1max, mag2max = data
        data_str = '%.2f %.2f %.2f %.2f' % (colmin, colmax, mag1max, mag2max)
    

    print data_str


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find color mag limits of CMDs interactively")

    parser.add_argument('-c', '--color_only', action='store_true',
                        help='skip the magnitude finding')

    parser.add_argument('-e', '--exgates', action='store_true',
                        help='Find exclude gates instead')
     
    parser.add_argument('-t', '--extpgates', action='store_true',
                        help='Find exclude tp-agb instead')

    parser.add_argument('-p', '--param', type=str, help='match param file')

    parser.add_argument('phot', type=str, help='match phot file')
    
    args = parser.parse_args(sys.argv[1:])
    
    mag1, mag2 = np.loadtxt(args.phot, unpack=True)
    
    if args.exgates:
        find_gates(mag1, mag2)
    #elif args.interactive:
    #    print 'NOT WORKING!'
    #    param = match_limits(mag1, mag2, color_only=args.color_only)
    #    match_diagnostic(param, args.phot)
    elif args.extpgates:
        exclude_tpagb(args.phot, args.param)
    else:
        magcolor([args.phot])
    
        
def magcolor(phots):
    fmt="""0.10 0.05 5 {cmin} {cmax} {filter1},{filter2}
{m1min} {m1max} {filter1}
{m2min} {m2max} {filter2}"""
    d = {}
    for m in phots:
        m1, m2 = np.loadtxt(m, unpack=True)
        d['filter1'], d['filter2'] = rsp.asts.parse_pipeline(m)[1]
        d['m1min'], d['m1max'] = m1.min(), m1.max()
        d['m2min'], d['m2max'] = m2.min(), m2.max()
        color = m1-m2
        d['cmin'], d['cmax'] = color.min(), color.max()
        print m
        print fmt.format(**d)
        print ''
        
