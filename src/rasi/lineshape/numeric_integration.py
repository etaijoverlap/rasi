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


from constants import hred
from math import sqrt,log,pi,exp
from scipy.special import gammaln
from tools import write_data


#################################################################################
#
# Computes the values of hermite functions up to order n at point x using the
# recurrence relation
#
#################################################################################

def hermite_function_rec(
                    n, # Quantum number
                    x  # normalized x value 
                    ):
    psi = [ 1./sqrt(sqrt(pi))*exp(-x**2/2) ,  2./sqrt(2*sqrt(pi))*exp(-x**2/2)*x ]
    for i in range(2,n+1):
        psi += [ sqrt(2./i)*x*psi[-1] - sqrt(float(i-1)/i)*psi[-2] ]
    
    return psi
 

def trapezoidal_eval(
                    # Physical parameters
                    omega_i,     # Parabolic const. of initial state
                    omega_f,     # Parabolic const. of final state
                    m,           # Mass of the system
                    x_s,         # Separation of parabolas
                    # Method parameters
                    maxchg,      # Max. Ansbacher recurrence error
                    maxref,      # Max. Refinement level
                    x_range,     # Normalized dist. of int. boundary from origin
                    n_eigs_i,    # Number of initial eigenstates to consider
                    n_eigs_f = None,
                                 # Number of final eigenstates to consider
                    minref = 2,  # Min. Refinement steps
                    prefix=None, # Prefix for output files
                    echo=True,   # Report status to stdout
                    hermite=hermite_function_rec,
                                 # Method to calculate hermite functions
                    nx_wf = 1000  # Number of samples for wavefuncion files
                    ):
   
    if n_eigs_f == None:
        n_eigs_f = n_eigs_i

    if echo:
        print
        print "Calculating numerical integrals of analytical solutions..."
        print "=========================================================="
        print
        print

    # Normalization constants for the hermite functions
    a_init = sqrt(m*omega_i/hred)
    a_finl = sqrt(m*omega_f/hred)
    
    # The subsequently used wave functions
    wavefuncs_init = lambda x: [ sqrt(a_init)*h for h in hermite(n_eigs_i-1,a_init*x)]
    wavefuncs_finl = lambda x: [ sqrt(a_finl)*h for h in hermite(n_eigs_f-1,a_finl*(x - abs(x_s)))]

    # The actual integration boundaries
    x_lower = - x_range/a_init
    x_upper =   x_range/a_finl + abs(x_s)

    # Sample the wavefunctions, mainly for debugging reasons
    if prefix != None:
        if echo:
            print;print "Writing wave functions...";print
        dx = (x_upper-x_lower)/nx_wf
        xx = x_lower
        xvals = [x_lower] 
        w_i = [ [ vv ] for vv in wavefuncs_init(xx) ]
        w_f = [ [ vv ] for vv in wavefuncs_finl(xx) ]
        for i in range(nx_wf):
            xx += dx
            w_i = [ ll + [vv] for ll,vv in zip(w_i,wavefuncs_init(xx)) ]
            w_f = [ ll + [vv] for ll,vv in zip(w_f,wavefuncs_finl(xx)) ]
            xvals += [ xx ]
        write_data(prefix+"_wf_i.dat",xvals,[xx*a_init for xx in xvals],*tuple(w_i))
        write_data(prefix+"_wf_f.dat",xvals,[(xx+abs(x_s))*a_finl for xx in xvals],*tuple(w_f))


    if echo:
        from units import Ang
        print "Prefactors: "
        print "----------"
        print "    a_init = %e ; a_finl = %e"%(a_init,a_finl)
        print;print
        print "Integration boundaries: "
        print "----------------------"
        print "    left = %e Angstroms; right = %e Angstroms"%(x_lower/Ang,x_upper/Ang)
        print;print

    # The function that will actually be integrated (returns a matrix of all the
    # possible wavefunction multiplications)
    _f = lambda wfi,wff: [ [ wf_i * wf_j for wf_j in wff ] for wf_i in wfi ]
    f = lambda x: _f(wavefuncs_init(x),wavefuncs_finl(x))

    # Some helper functions for overlap matrix handling
    add = lambda a,b: [ [ x + y for x,y in zip(vx,vy) ] for vx,vy in zip(a,b) ]
    mult = lambda a,mx: [ [ a*x  for x in vx ] for vx in mx ]
    diff = lambda mx,my: [ [ x - y for x,y in zip(vx,vy) ] for vx,vy in zip(mx,my) ]
    norm = lambda mx: sum(sum(x**2 for x in vx) for vx in mx)
    mmax = lambda mx: max(max(vx) for vx in mx)
    mmin = lambda mx: min(min(vx) for vx in mx)
   
    # Initial computation
    S = mult(.5*(x_upper-x_lower),add(f(x_upper),f(x_lower)))

    normS=norm(S)
    print "Initial computation norm (should be as low as possible):"
    print "-------------------------------------------------------"
    print "                   norm(S) = %e"%normS
    if normS > 1e-5:
        print "WARNING: Your initial S matrix norm is rather high, you should"
        print "         consider increasing the boundary distance parameter"
        print "         ( now: %e )"%x_range
    print;print
    if echo:
        print "%10s %20s %20s %20s"%("Ref.Step","Diff.Norm","Max.Ovrlp","Recurrence Err.")
        print "-"*70 
    iref=0
    beta = a_finl/a_init
    gamma = x_s*a_init
    while(True):
        iref += 1
        nsteps = 1 << (iref-1)
        dx = (x_upper-x_lower)/nsteps
        xx = x_lower + .5*dx
        Simprove = f(xx)
        for i in range(1,nsteps): Simprove = add(Simprove,f(xx+i*dx))
        S_new = mult(.5,add(S,mult(dx,Simprove)))
        nrm = norm(diff(S_new,S))
        S = S_new
        err = 0.
        lrgst = 0.
        lrgidx = ()
        for i in range(2,len(S)):
            for j in range(1,len(S[0])):
                Sans = beta**2*gamma/(1+beta**2)*sqrt(2./i)*S[i-1][j]+2*beta/(1+beta**2)*sqrt(float(j)/i)*S[i-1][j-1]+(1-beta**2)/(1+beta**2)*sqrt(float(i-1)/i)*S[i-2][j]
                ierr = (S[i][j] - Sans)**2
                err += ierr

        if echo:
            print "%10d %20e %20f %20e"%(iref,nrm,mmax(S),err)
        
        if iref > minref:
            if iref > maxref:
                if echo and err > maxchg:
                    print "Warning: Could not reach convergence!"
                break
            if err < maxchg:
                if echo:
                    print "Seem to have reached sufficient accuracy..."
                break
    Ew_i = [ (i+.5)*hred*omega_i for i in range(n_eigs_i) ]
    Ew_f = [ (i+.5)*hred*omega_f for i in range(n_eigs_f) ]
    return (Ew_i,Ew_f,S)

if __name__ == "__main__":
    from tools import write_data
    from units import *
    #print "Testing mode!!!"
    mass = 24.6619642123*amu
    #x_s = 0.
    x_s = 0.55779*Ang 
    k_i = 2.56092202778*eV/Ang**2
    #k_f = k_i
    k_f = 6.29591112434*eV/Ang**2
    Es = +(2.31509-5.06072699999997)*eV
    omega_i = (2*k_i/mass)**.5
    omega_f = (2*k_f/mass)**.5
    (_,_,S) = trapezoidal_eval( 
                   # Physical parameters
                   omega_i,    # Parabolic const. of initial state
                   omega_f,    # Parabolic const. of final state
                   mass,          # Mass of the system
                   x_s,        # Separation of parabolas
                   # Method parameters
                   1e-16,# Max. Ansbacher recurrence error
                   14,   # Max. Refinement level
                   25.,  # Boundary distance
                   300,   # Number of initial eigenstates to consider
                   #minref=13,
                   echo=True,     # Report status to stdout
                   prefix="sn"
                   )
    from zapol import ZAPOL_overlaps
    (_,_,Szap) = ZAPOL_overlaps(
                   # Physical parameters
                   omega_i,    # Parabolic const. of initial state
                   omega_f,    # Parabolic const. of final state
                   mass,          # Mass of the system
                   x_s,        # Separation of parabolas
                   # Method parameters
                   20,   # Number of initial eigenstates to consider
                   echo=True     # Report status to stdout
                   )
    diff = lambda mx,my: [ [ x - y for x,y in zip(vx,vy) ] for vx,vy in zip(mx,my) ]
    norm = lambda mx: sum(sum(x**2 for x in vx) for vx in mx)
    print "ZAPOL Diference: %e"%(norm(diff(Szap,S)))
    sfile = open("Smatrix.dat","w")
    for line in S:
        l = ""
        for val in line:
            l += "%15.7e"%float(val)
        print >>sfile, l

    sfile = open("Smatrix_z.dat","w")
    for line in Szap:
        l = ""
        for val in line:
            l += "%15.7e"%float(val)
        print >>sfile, l
