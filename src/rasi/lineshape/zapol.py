from constants import *

def a_template(binom):
    def a(k,l,j):
        t = j+k+l
        if t%2:
            return 0
        else:
            t /= 2
            pmin = max(0,(k-j-l)/2)
            pmax = min(k,(k-j+l)/2)
            return sum( (-1)**p*binom(k,p)*binom(l,t-j-p) for p in range(pmin,pmax+1))
    return a

def h_table (n,xval): 
    h = [ 1 , 2*xval ]
    for i in range(2,n+1):
        h += [ 2*xval*h[-1] - 2*(i-1)*h[-2] ]
    return h

fast_binom = lambda factab: \
                lambda n,k: factab[n]/(factab[k]*factab[n-k])

def ZAPOL_overlaps( 
                   # Physical parameters
                   omega_i,    # Parabolic const. of initial state
                   omega_f,    # Parabolic const. of final state
                   m,          # Mass of the system
                   x_s,        # Separation of parabolas
                   # Method parameters
                   n_eigs_i,   # Number of initial eigenstates to consider
                   n_eigs_f = None,
                               # Number of final eigenstates to consider
                   echo=True     # Report status to stdout
                   ):
    from math import log,sin,cos,exp,atan
    from scipy.special import gammaln
    from sys import stdout

    if echo:
        print
        print "Calculating analytical integrals after B.P. ZAPOL "
        print "=================================================="
        print
        print

    if n_eigs_f == None:
        n_eigs_f = n_eigs_i
    
    if echo:
        print;print "Building factorial table...";print
    factab = [ 1 ]
    for i in range(1,max(n_eigs_i,n_eigs_f)+1):
        factab += [ factab[-1] * i ]
    binom = fast_binom(factab)
    a = a_template(binom)        
    ld_e = 1/log(2)
    theta = atan((omega_f/omega_i)**.5)
    s = sin(theta)
    c = cos(theta)

    omega = omega_i*omega_f/(omega_i+omega_f)
    rho = (m*omega/(2*hred))**.5 *x_s

    if echo:
        print;print "Building Hermite table...";print
    H = h_table(n_eigs_i+n_eigs_f,rho)

    if echo:
        print;print "Calculating overlaps...";print
    LS = []

    Ew_i = [ hred*omega_i*(.5+m) for m in range(n_eigs_i) ]
    Ew_f = [ hred*omega_f*(.5+m) for m in range(n_eigs_f) ]

    amn = []
    S = []
    pct=(n_eigs_i-1)/10
    nxt=pct
    if echo:
        stdout.write("0%...");stdout.flush()
    for m in range(n_eigs_i):
        if echo and m >= nxt:
            if m/pct < 10:
                stdout.write("%d0%%..."%(m/pct));stdout.flush()
            else:
                stdout.write("100%\n");stdout.flush()
            nxt+=pct
        S_i = []
        for n in range(n_eigs_f):
            if n == 0:
                amn = [ a(m,n,j) for j in range(m+n+1) ]
            else:
                amn = [ a(m,n,0) ] + [ amn[j-1]+amn[j+1] for j in range(1,m+n-1) ]
                if m+n > 1:
                    amn += [ a(m,n,m+n-1) ]
                if m+n > 0:
                    amn += [ a(m,n,m+n) ]
            sgn=1
            if m%2:
                I = 0.
                for j in range(1,m+n+1):
                    u=j+m+n
                    if u%2:
                        continue
                    else:
                        u/=2
                        I += amn[j]*sin(j*theta)*H[u]*H[u-j]
                Ipre = (1-(m+n)/2.)/ld_e
                sgn*=(-1)**((m+1)/2)
            else:
                I= 0.
                for j in range(1,m+n+1):
                    u=j+m+n
                    if u%2:
                        continue
                    else:
                        u/=2
                        I += amn[j]*cos(j*theta)*H[u]*H[u-j]
                I *= 2
                if n%2 == 0:
                    I += a(m,n,0)*H[(m+n)/2]**2
                Ipre = (-(m+n)/2.)/ld_e
                sgn*=(-1)**(m/2)
            integ = sgn*(-1)**(m+n)*sin(2*theta)**.5*I*exp(Ipre-((m+n)/ld_e+gammaln(m+1)+gammaln(n+1))/2-rho**2)
            Ef = hred*omega_f*(n+.5)
            S_i.append(integ)
        S.append(S_i)
    return Ew_i,Ew_f,S


if __name__=="__main__":
    from units import *
    #print "Testing mode!!!"
    mass = 24.6619642123*amu
    x_s = 0.55779*Ang 
    k_i = 2.56092202778*eV/Ang**2
    k_f = 6.29591112434*eV/Ang**2
    Es = +(2.31509-5.06072699999997)*eV
    omega_i = (2*k_i/mass)**.5
    omega_f = (2*k_f/mass)**.5
    (Ew,Ef,S) = ZAPOL_overlaps( 
                   # Physical parameters
                   omega_i,    # Parabolic const. of initial state
                   omega_f,    # Parabolic const. of final state
                   mass,          # Mass of the system
                   x_s,        # Separation of parabolas
                   # Method parameters
                   50,   # Number of initial eigenstates to consider
                   50,
                               # Number of final eigenstates to consider
                   echo=True     # Report status to stdout
                   )
    sfile = open("Smatrix.dat","w")
    for line in S:
        l = ""
        for val in line:
            l += "%15.7e"%float(val)
        print >>sfile, l

