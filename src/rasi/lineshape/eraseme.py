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


from scipy import array
from constants import *
from tools import write_data


def hermite_table (n,xvals): 
    h = [ array([  1  for x in xvals ]), 2*xvals ]
    for i in range(2,n+1):
        h += [ 2*xvals*h[-1] - 2*(i-1)*h[-2] ]
    return h


def semi_numeric_overlaps( 
                        # Physical parameters
                        omega_i,    # Parabolic const. of initial state
                        omega_f,    # Parabolic const. of final state
                        m,          # Mass of the system
                        x_s,        # Separation of parabolas
                        # Method parameters
                        dx,         # Spatial sample width
                        x_i_range,  # Distance of boundary to initial minimum
                        x_f_range,  # Distance of boundary to final minimum
                        n_eigs_i,   # Number of initial eigenstates to consider
                        n_eigs_f = None,
                                    # Number of final eigenstates to consider
                        x_i_0 = 0., # Position of initial minimum
                        n_write=0,    # Number of wavefunctions to write
                        prefix=None,  # Prefix for output
                        echo=True     # Report status to stdout
                        ):

    from tools import gen_grid
    from scipy import exp

    if echo:
        print
        print "**********************************************************"
        print "Calculating numerical integrals of analytical solutions..."
        print
        print

    x_f_0 = x_i_0+x_s

    if n_eigs_f == None:
        n_eigs_f = n_eigs_i

    boundaries = [ x_i_0 - x_i_range, x_i_0 + x_i_range, 
                   x_f_0 - x_f_range, x_f_0 + x_f_range ]

    x_left  = min(boundaries)
    x_right = max(boundaries)

    v_x = gen_grid(x_left,x_right,dx)
    xvals = array(v_x)

    if echo:
        print "Generating factorial table..."
    facttab = [1]
    for i in range(1,max(n_eigs_i,n_eigs_f)):
        facttab += [ i*facttab[-1] ]
    facttab = array(facttab)

    if echo: 
        print "Generating Hermite polynomial table..."
    prefactor_i = (m*omega_i/hred)**.5
    prefactor_f = (m*omega_f/hred)**.5
    hermite_i = hermite_table(n_eigs_i,prefactor_i*(xvals-x_i_0))
    hermite_f = hermite_table(n_eigs_f,prefactor_f*(xvals-x_f_0))

    if echo:
        print "Generating Gaussians..."
    gauss_i = exp(-(prefactor_i*(xvals-x_i_0))**2/2)
    gauss_f = exp(-(prefactor_f*(xvals-x_f_0))**2/2)

    if echo:
        print "Generating Wavefunctions..."
    w_i = [ (2.**n*f)**(-.5)*prefactor_i**.5/pi**.25*gauss_i*h
                                for n,(f,h) in enumerate(zip(facttab,hermite_i)) ]
    w_f = [ (2.**n*f)**(-.5)*prefactor_f**.5/pi**.25*gauss_f*h 
                                for n,(f,h) in enumerate(zip(facttab,hermite_f)) ]

    Ew_i = [ (i+.5)*hred*omega_i for i in range(n_eigs_i) ]
    Ew_f = [ (i+.5)*hred*omega_f for i in range(n_eigs_f) ]

    if prefix != None:
        if n_write > 0:
            write_data("semi_numeric_wf_i.dat",xvals,*tuple(w_i[:n_write]))
            write_data("semi_numeric_wf_f.dat",xvals,*tuple(w_f[:n_write]))
        write_data("semi_numeric_eigs.dat", Ew_i, Ew_f)

    if echo:
        print "Calculating overlaps..."
    S = []
    for i in range(n_eigs_i):
        S_i = []
        for j in range(n_eigs_f):
            S_i.append(((w_i[i]*w_f[j]).sum()*dx)**2)
        S.append(S_i)
    return Ew_i,Ew_f,S
def hermite_function(
                    n, # Quantum number
                    a, # sqrt(m*omega/hbar)
                    xvals,  # list of normalized x values (x' = a*x)
                    i_log_from = 90, # Specify from which number to switch to logarithmic mode 
                    x_log_from = 1e100 # Specify when the x**m part should be switched to log mode
                    ):
    prf0 = sqrt(a/sqrt(pi)) #The nicer prefactor

    fact = facttable(min(n,i_log_from))

    logfact = lambda i: log(fact[i]) if i < i_log_from else stirling(i)

    psi = []
    for x in xvals:
        psi_x = []
        for m in range(n/2+1):
            # First get the sign right
            psi_x_m = (-1)**m
            if n%2 == 1 and x < 0:
                psi_x_m *= -1
            ax = abs(x)

            # Now get done with all those factorials
            fact_part = 0.
            fact_log = False
            if n > i_log_from or False:
                fact_log = True
                fact_part = stirling(n)/2-(logfact(m)+logfact(n-2*m))
            else:
                fact_part = sqrt(fact[n])/(fact[m]*fact[n-2*m])

            # Care for the exponent of 2 part
            two_part = (float(n)/2-2*m)
            two_log  = False
            if two_part < 300 or True:
                two_part = 2.**two_part
                two_log = False
            else:
                two_part = two_part*log(2)
                two_log = True

            # Finally, the spatial expressions
            x_part = 0.
            x_log  = False
            if n-2*m != 0 and ax > x_log_from**(1./(n-2*m) and False):
                x_log=True
                x_part = (n-2*m)*log(ax)
            else:
                x_part = ax**(n-2*m)

            # Now, let's tuck all this together
            log_part = 0.
            lin_part = 1.
            if fact_log: log_part += fact_part
            else: lin_part *= fact_part

            if two_log: log_part += two_part 
            else: lin_part *= two_part

            if x_log: log_part += x_part 
            else: lin_part *= x_part

            psi_x_m *= lin_part*exp(log_part-ax**2/2)
            psi_x += [ psi_x_m ]
   #     print "Summands: ",psi_x
        psi += [ prf0*sum(psi_x) ]
    
    return psi

def facttable(n):
    ft = [ 1 ]
    for i in range(1,n+1):
        ft += [ ft[-1]*i ]
    return ft

def stirling(n):
    return n*log(n)-n

factorial_simple = lambda n: 1 if n==0 else n*factorial_simple(n-1)

def hermite_poly_simple(n,xvals):
    h = [ [ 1 for x in xvals ] , [ 2*x for x in xvals ]]
    for i in range(2,n+1):
        h += [[ 2*x*hn - 2*(i-1)*hnn for x,hn,hnn in zip(xvals,h[-1],h[-2]) ]]
    return h[-1]

def hermite_function_simple(n,a,xvals):
    return [ a/(sqrt(factorial_simple(n)*2**n*sqrt(pi)))*exp(-x**2/2)*h for x,h in zip(xvals,hermite_poly_simple(n,xvals)) ]

