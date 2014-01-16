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


from units import *
from constants import *
from math import log,sin,cos,exp,atan

def a_template(binom):
    def a(k,l,j):
        t = j+k+l
        if t%2:
            return 0
        else:
            t /= 2
            pmin = max(0,(k-j-l)/2)
            pmax = min(k,(k-j+l)/2)
            res=sum( (-1)**p*exp(binom(k,p)+binom(l,t-j-p)) for p in range(pmin,pmax+1))
            if res < 0:
                raise ValueError("My a is negative!!!")
            else:
                return log(res)
    return a

# def h_table (n,xval): 
#     h = [ 0. , log(2*xval) ]
#     for i in range(2,n+1):
#         h += [ h[-1]+log(2*xval - 2*(i-1)*exp(h[-1]-h[-2])) ]
#     return h

def h_table (n,xval): 
    h = [ 1 , 2*xval ]
    for i in range(2,n+1):
        h += [ h[-1]*(2*xval - 2*(i-1)*h[-2]/h[-1]) ]
    return h

fast_binom = lambda factab: \
                lambda n,k: factab[n]-(factab[k]+factab[n-k])

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
    from scipy.special import gammaln

    if echo:
        print
        print "**********************************************************"
        print "Calculating analytical integrals after B.P. ZAPOL ..."
        print
        print

    if n_eigs_f == None:
        n_eigs_f = n_eigs_i
    
    if echo:
        print "Building factorial table..."
    factab = [ 1 ]
    for i in range(1,max(n_eigs_i,n_eigs_f)+1):
        factab += [ factab[-1] * i ]

    # Switch to logarithmic mode
    factab = [ log(x) for x in factab ]
    binom = fast_binom(factab)
    a = a_template(binom)
    ld_e = 1/log(2)
    theta = atan((omega_f/omega_i)**.5)
    s = sin(theta)
    c = cos(theta)

    omega = omega_i*omega_f/(omega_i+omega_f)
    rho = (m*omega/(2*hred))**.5 *x_s

    if echo:
        print "Building Hermite table..."
    H = h_table(n_eigs_i+n_eigs_f,rho)
    print H

    Hsgn = [ (1 if x >= 0 else -1) for x in H ]
    H = [ abs(x) for x in H ]
    return

    if echo:
        print "Calculating overlaps..."
    LS = []

    Ew_i = [ hred*omega_i*(.5+m) for m in range(n_eigs_i) ]
    Ew_f = [ hred*omega_f*(.5+m) for m in range(n_eigs_f) ]

    amn = []
    S = []
    for m in range(n_eigs_i):
        S_i = []
        for n in range(n_eigs_f):
            try:
                if n == 0:
                    amn = [ a(m,n,j) for j in range(m+n+1) ]
                else:
                    amn = [ a(m,n,0) ] + [ log(exp(amn[j-1])+exp(amn[j+1])) for j in range(1,m+n-1) ]
                    if m+n > 1:
                        amn += [ a(m,n,m+n-1) ]
                    if m+n > 0:
                        amn += [ a(m,n,m+n) ]
                sgn=1
                if m%2:
                    I = 0.
                    Ipre = (1-(m+n)/2.)/ld_e
                    for j in range(1,m+n+1):
                        u=j+m+n
                        if u%2:
                            continue
                        else:
                            u/=2
                            I += sin(j*theta)*Hsgn[u]*Hsgn[u-j]*exp(amn[j]+H[u]+H[u-j]+Ipre-((m+n)/ld_e+gammaln(m+1)+gammaln(n+1))/2-rho**2)
                    sgn*=(-1)**((m+1)/2)
                else:
                    I= 0.
                    Ipre = (-(m+n)/2.)/ld_e
                    for j in range(1,m+n+1):
                        u=j+m+n
                        if u%2:
                            continue
                        else:
                            u/=2
                            I += cos(j*theta)*Hsgn[u]*Hsgn[u-j]*exp(amn[j]+H[u]+H[u-j]+Ipre-((m+n)/ld_e+gammaln(m+1)+gammaln(n+1))/2-rho**2)
                    I *= 2
                    if n%2 == 0:
                        I += exp(a(m,n,0)*2*H[(m+n)/2]+Ipre-((m+n)/ld_e+gammaln(m+1)+gammaln(n+1))/2-rho**2)
                    #print "n= %d m=%d  I=%e  Ipre=%e"%(n,m,I,Ipre)
                    sgn*=(-1)**(m/2)
                integ = sgn*(-1)**(m+n)*sin(2*theta)**.5*I
                Ef = hred*omega_f*(n+.5)
                try:
                    S_i.append(integ**2)
                except OverflowError,e:
                    print "Overflow at n=%d ; m=%d ; integ=%e ; I=%e ; Ipre=%e"%(n,m,integ,I,Ipre)
                    return   
            except ValueError,e:
                print e
                print "Overflow at n=%d ; m=%d"%(n,m)
                return
        S.append(S_i)
    return Ew_i,Ew_f,S

if __name__=="__main__":
    print "Testing mode!!!"
    mass = 24.6619642123*amu
    x_s = 0.55779*Ang 
    k_i = 2.56092202778*eV/Ang**2
    k_f = 6.29591112434*eV/Ang**2
    Es = +(2.31509-5.06072699999997)*eV
    omega_i = (2*k_i/mass)**.5
    omega_f = (2*k_f/mass)**.5
    ZAPOL_overlaps( 
                   # Physical parameters
                   omega_i,    # Parabolic const. of initial state
                   omega_f,    # Parabolic const. of final state
                   mass,          # Mass of the system
                   x_s,        # Separation of parabolas
                   # Method parameters
                   150,   # Number of initial eigenstates to consider
                   150,
                               # Number of final eigenstates to consider
                   echo=True     # Report status to stdout
                   )
