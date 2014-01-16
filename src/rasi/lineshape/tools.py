###############################################################################
#
#  RASI ... Reliability Analysis of Semiconductor Interfaces
#
#  (c) 2013-2014 Franz Schanovsky 
#  
#  This project is dedicated to the loving memory of Margarete and Johann 
#  Mittermayr.
#
###############################################################################
#
#    This software is licensed under the EUPL V 1.1
#
#    This software is provided "as is" without warranty of any kind, see the 
#    respective section in the EUPL. USE AT YOUR OWN RISK.
#
###############################################################################


from constants import *

def gen_grid(xmin,xmax,dx):
    x = [ xmin ]
    xval = xmin
    while xval < xmax:
        xval += dx
        x += [ xval ]
    return x

def group(LS,         # Line shape list consisting of (E,p) pairs
          group_dist, # maximum distance for grouping
          echo=True   # Report status to stdout
          ):
    if echo:
        print 
        print "Grouping near lines"
        print "-------------------"
        print

    ielim = 0
    LS_n = [ LS[0] ]
    for (E,p) in LS[1:]:
        E_old,p_old = LS_n[-1]
        if E - E_old < group_dist:
            LS_n[-1] = (E_old,p+p_old)
            ielim += 1
        else:
            LS_n += [ (E,p) ]
    if echo:
        print "     %d lines eliminated"%ielim
    return LS_n

def left_right(x0,xs,ys):
    x_old = xs[0]
    y_old = ys[0]
    for x,y in zip(xs[1:],ys[1:]):
        if x > x0:
            return (x_old,y_old,x,y)

def interpolate(x,xs,ys):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    (xl,yl,xr,yr) = left_right(x,xs,ys)
    k = (yr-yl)/(xr-xl)
    return yl+k*(x-xl)

def write_data(filename,*columns):
    outfile = open(filename,"w")
    for line in zip(*columns):
        print >> outfile, ("%e  "*len(line))%line
    outfile.close()

def write_matrix(filename,matrix):
    outfile = open(filename,"w")
    iters = [ vec.__iter__() for vec in matrix ]
    try:
        while True:
            for iter in iters:
                outfile.write("%e  "%iter.next())
            outfile.write("\n")
    except StopIteration: pass


def smear(xy,    # (x,y) pairs
          dx,    # Separation of spatial samples
          sigma, # Standard deviation of the smearing gaussian
          echo=True # Report status to stdout
          ):
    from scipy import array,exp
    xvals = [ x for x,y in xy ]
    xmin = min(xvals)
    xmax = max(xvals)
    if echo:
        from units import eV
        print
        print "Applying gaussian smearing"
        print "--------------------------"
        print
        print "   sigma = %20.5e"%(sigma/eV)
        print "      dx = %20.5e"%(dx/eV)
        print "    xmin = %20.5e"%(xmin/eV)
        print "    xmax = %20.5e"%(xmax/eV)
    xvals  = array(gen_grid(xmin,xmax,dx))
    ysmear = array([ 0. for _ in xvals ])
    for x0,y0 in xy:
        ysmear += (y0*(2.*pi*sigma**2)**-.5)*exp(-(xvals-x0)**2/(2*sigma**2)) 
    return xvals,ysmear

def calc_lsf(Ew_i, # Eigenvalues of initial state
             Ew_f, # Eigenvalues of final state
             S,    # Overlap matrix
             T,    # Temperature
             Z=None, # Partition function (calculated if not given)
             e_capture = True, # Electron capture transition (hole capture if false)
             echo=True # Report status to stdout
             ):
    from math import exp

    if echo:
        print
        print "Computing line shape function"
        print "-----------------------------"
        print

    if Z == None:
        Z = sum( exp(-E/(kB*T)) for E in Ew_i)
        if echo:
            print "Partition function for initial state (calculated): ",Z
            print
    else:
        if echo:
            print "Partition function for initial state (given): ",Z
            print
    LS=[]
    for E_init,S_i in zip(Ew_i,S):
        p0 = exp(-E_init/(kB*T))/Z
        for E_final,o in zip(Ew_f,S_i):
            if e_capture:
                Et = E_final-E_init # Electron capture gives a peak at +Delta E
            else:
                Et = E_init-E_final # Hole capture gives a peak at -Delta E
            LS.append((Et,p0*o**2))
    LS.sort(lambda (E1,_1),(E2,_2): 1 if E1>E2 else -1)
    return LS

from math import exp
harmosc_pf = lambda omega,T: exp(-hred*omega/(kB*T)*1./2)/(1.-exp(-hred*omega/(kB*T)))


def compare_results((Ew_i_1,Ew_f_1,S_1),(Ew_i_2,Ew_f_2,S_2)):
    print 
    print "Comparing results"
    print "-----------------"
    diff = lambda mx,my: [ [ x - y for x,y in zip(vx,vy) ] for vx,vy in zip(mx,my) ]
    norm = lambda mx: sum(sum(x**2 for x in vx) for vx in mx)

    len_i_1 = len(Ew_i_1)
    len_i_2 = len(Ew_i_2)
    len_f_1 = len(Ew_f_1)
    len_f_2 = len(Ew_f_2)

    if len_i_1 != len_i_2:
        print "   Number of initial values differ by  %d"%abs(len_i_1-len_i_2)
    if len_f_1 != len_f_2:
        print "   Number of final   values differ by  %d"%abs(len_f_1-len_f_2)

    E_init_err = sum(abs(E1-E2) for E1,E2 in zip(Ew_i_1,Ew_i_2))
    E_finl_err = sum(abs(E1-E2) for E1,E2 in zip(Ew_f_1,Ew_f_2))
    S_err = norm(diff(S_1,S_2))

    print "          %20s"%"Error" 
    print " E init   %20.5e"%E_init_err
    print " E finl   %20.5e"%E_finl_err
    print " Ovl.Mat. %20.5e"%S_err

   
