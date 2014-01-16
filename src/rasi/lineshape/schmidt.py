from constants import hred
from math import sqrt,log,pi,exp
from scipy.special import gammaln
from tools import write_data


def schmidt_recurrence_relations(
                    # Physical parameters
                    omega_i,     # Parabolic const. of initial state
                    omega_f,     # Parabolic const. of final state
                    m,           # Mass of the system
                    x_s,         # Separation of parabolas
                    # Method parameters
                    n_eigs_i,    # Number of initial eigenstates to consider
                    n_eigs_f = None,
                                 # Number of final eigenstates to consider
                    echo=True,   # Report status to stdout
                    ):
   
    if n_eigs_f == None:
        n_eigs_f = n_eigs_i

    if echo:
        print
        print "Using recurrence relations by P. P. Schmidt"
        print "==========================================="
        print
        print

    # Just change the symbols to those used in Schmidt's paper
    a_init = m*omega_i/hred
    a_finl = m*omega_f/hred
    delta = x_s
    if echo:
        print " Schmidt's parameters:"
        print " -------------------- "
        print
        print " a     = %17.5e"%a_init
        print " a'    = %17.5e"%a_finl
        print " delta = %17.5e"%delta
        print
        print

    # Allocate the mighty mighty overlap array
    I = [ [ 0. for i in xrange(n_eigs_f) ] for j in xrange(n_eigs_i) ]

    # Calculate the initial values
    I[0][0] = sqrt(2*sqrt(a_init*a_finl)/(a_init+a_finl))*exp(-(a_init*a_finl*delta**2)/(2*(a_init+a_finl)))
    I[1][0] = (sqrt(2*a_init)*a_finl)/(a_init+a_finl)*delta*I[0][0]
    I[0][1] = -(sqrt(2*a_finl)*a_init)/(a_init+a_finl)*delta*I[0][0]

    # Now calculate along the edges
    for m in range(0,n_eigs_i-2):
        I[m+2][0] =  sqrt((float(m)+1)/(float(m)+2)) * (a_init-a_finl)/(a_init+a_finl) * I[m][0] \
                           + sqrt(2./(float(m)+2)) * sqrt(a_init)*a_finl/(a_init+a_finl) * delta * I[m+1][0]
    for n in range(0,n_eigs_f-2):
        I[0][n+2] =  sqrt((float(n)+1)/(float(n)+2)) * (a_finl-a_init)/(a_init+a_finl) * I[0][n] \
                           - sqrt(2./(float(n)+2)) * a_init*sqrt(a_finl)/(a_init+a_finl) * delta * I[0][n+1]

    def do_fill(m,n):
        II = lambda a,b: 0. if a<0 or b<0 else I[a][b]
        I[m+1][n+1] = - a_init*sqrt(a_finl)*delta/(sqrt(2)*(a_init+a_finl))*sqrt(1./(float(n)+1))*II(m+1,n) \
                     - (a_init-a_finl)/(a_init+a_finl)*sqrt(float(n)/(float(n)+1))*II(m+1,n-1) \
                     + sqrt(a_init)*a_finl*delta/(sqrt(2)*(a_init+a_finl))*sqrt(1./float(m+1))*II(m,n+1) \
                     + 2*sqrt(a_init*a_finl)/(a_init+a_finl)*sqrt(1./(float(m+1)*float(n+1)))*II(m,n) \
                     - sqrt(a_init)*a_finl*delta/(sqrt(2)*(a_init+a_finl))*sqrt(float(n)/(float(m+1)*float(n+1)))*II(m,n-1) \
                     + (a_init-a_finl)/(a_init+a_finl)*sqrt(float(m)/float(m+1))*II(m-1,n+1)\
                     + a_init*sqrt(a_finl)*delta/(sqrt(2)*(a_init+a_finl))*sqrt(float(m)/(float(m+1)*float(n+1)))*II(m-1,n) \
                     + sqrt(float(m*n)/(float(m+1)*float(n+1)))*II(m-1,n-1)
    for i in range(n_eigs_i-1):
        for j in range(min(i+1,n_eigs_f-1)):
                do_fill(i-j,j)
            
    for j in range(1,n_eigs_f-1):
        c = 0
        for i in range(n_eigs_i-2,-1,-1):
            m = i
            n = j+c
            if n > n_eigs_f-2:
                break
            do_fill(m,n)
            c += 1

    Ew_i = [ (i+.5)*hred*omega_i for i in range(n_eigs_i) ]
    Ew_f = [ (i+.5)*hred*omega_f for i in range(n_eigs_f) ]
    return (Ew_i,Ew_f,I)
