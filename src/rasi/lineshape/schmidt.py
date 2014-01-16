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


class SchmidtOverlaps(object):
    def __init__(self, omega_occupied=None, omega_unoccupied=None, mass=None, equilibrium_shift=None, n_states_occupied=None, n_states_unoccupied=None):
        """
                Properties:
        """ 
        self.__omega_occupied      = None
        self.__omega_unoccupied    = None
        self.__mass                = None
        self.__equilibrium_shift   = None
        self.__n_states_occupied   = 200
        self.__n_states_unoccupied = 200

        self.__changed = False

        if omega_occupied      != None: self.omega_occupied      = omega_occupied
        if omega_unoccupied    != None: self.omega_unoccupied    = omega_unoccupied
        if mass                != None: self.mass                = mass
        if equilibrium_shift   != None: self.equilibrium_shift   = equilibrium_shift
        if n_states_occupied   != None: self.n_states_occupied   = n_states_occupied
        if n_states_unoccupied != None: self.n_states_unoccupied = n_states_unoccupied

    def update(self):
        from scipy.constants import hbar
        from math import sqrt,log,pi,exp
        from scipy import arange,array

        changed = self.__changed
        self.__changed = False
        if changed :
            omega_occupied = self.omega_occupied ; omega_unoccupied = self.omega_unoccupied
            mass = self.mass ; equilibrium_shift = self.equilibrium_shift
            n_states_occupied = self.n_states_occupied ; n_states_unoccupied = self.n_states_unoccupied

            # Just change the symbols to those used in Schmidt's paper
            a      = mass*omega_occupied/hbar
            aprime = mass*omega_unoccupied/hbar
            delta  = equilibrium_shift

            # Allocate the mighty mighty overlap array
            I = [ [ 0. for i in xrange(n_states_unoccupied) ] for j in xrange(n_states_occupied) ]

            # Calculate the initial values
            I[0][0] = sqrt(2*sqrt(a*aprime)/(a+aprime))*exp(-(a*aprime*delta**2)/(2*(a+aprime)))
            I[1][0] = (sqrt(2*a)*aprime)/(a+aprime)*delta*I[0][0]
            I[0][1] = -(sqrt(2*aprime)*a)/(a+aprime)*delta*I[0][0]

            # Now calculate along the edges
            for m in range(0,n_states_occupied-2):
                I[m+2][0] =  sqrt((float(m)+1)/(float(m)+2)) * (a-aprime)/(a+aprime) * I[m][0] \
                                   + sqrt(2./(float(m)+2)) * sqrt(a)*aprime/(a+aprime) * delta * I[m+1][0]
            for n in range(0,n_states_unoccupied-2):
                I[0][n+2] =  sqrt((float(n)+1)/(float(n)+2)) * (aprime-a)/(a+aprime) * I[0][n] \
                                   - sqrt(2./(float(n)+2)) * a*sqrt(aprime)/(a+aprime) * delta * I[0][n+1]

            def do_fill(m,n):
                II = lambda a,b: 0. if a<0 or b<0 else I[a][b]
                I[m+1][n+1] = - a*sqrt(aprime)*delta/(sqrt(2)*(a+aprime))*sqrt(1./(float(n)+1))*II(m+1,n) \
                             - (a-aprime)/(a+aprime)*sqrt(float(n)/(float(n)+1))*II(m+1,n-1) \
                             + sqrt(a)*aprime*delta/(sqrt(2)*(a+aprime))*sqrt(1./float(m+1))*II(m,n+1) \
                             + 2*sqrt(a*aprime)/(a+aprime)*sqrt(1./(float(m+1)*float(n+1)))*II(m,n) \
                             - sqrt(a)*aprime*delta/(sqrt(2)*(a+aprime))*sqrt(float(n)/(float(m+1)*float(n+1)))*II(m,n-1) \
                             + (a-aprime)/(a+aprime)*sqrt(float(m)/float(m+1))*II(m-1,n+1)\
                             + a*sqrt(aprime)*delta/(sqrt(2)*(a+aprime))*sqrt(float(m)/(float(m+1)*float(n+1)))*II(m-1,n) \
                             + sqrt(float(m*n)/(float(m+1)*float(n+1)))*II(m-1,n-1)
            for i in range(n_states_occupied-1):
                for j in range(min(i+1,n_states_unoccupied-1)):
                        do_fill(i-j,j)
                    
            for j in range(1,n_states_unoccupied-1):
                c = 0
                for i in range(n_states_occupied-2,-1,-1):
                    m = i
                    n = j+c
                    if n > n_states_unoccupied-2:
                        break
                    do_fill(m,n)
                    c += 1

            self.energies_occupied   =   (arange(float(n_states_occupied))+.5)*hbar*omega_occupied
            self.energies_unoccupied = (arange(float(n_states_unoccupied))+.5)*hbar*omega_unoccupied
            self.overlap_matrix      = array(I)
        return changed
        
    def get_omega_occupied(self):
        return self.__omega_occupied
    def set_omega_occupied(self,omega):
        self.__omega_occupied = omega
        self.__changed = True
    omega_occupied = property(get_omega_occupied,set_omega_occupied)

    def get_omega_unoccupied(self):
        return self.__omega_unoccupied
    def set_omega_unoccupied(self,omega):
        self.__omega_unoccupied = omega
        self.__changed = True
    omega_unoccupied = property(get_omega_unoccupied,set_omega_unoccupied)

    def get_mass(self):
        return self.__mass
    def set_mass(self,m):
        self.__mass = m
        self.__changed = True
    mass = property(get_mass,set_mass)

    def get_equilibrium_shift(self):
        return self.__equilibrium_shift
    def set_equilibrium_shift(self,e):
        self.__equilibrium_shift = e
        self.__changed = True
    equilibrium_shift = property(get_equilibrium_shift,set_equilibrium_shift)

    def get_n_states_occupied(self):
        return self.__n_states_occupied
    def set_n_states_occupied(self,n):
        self.__n_states_occupied = n
        self.__changed = True
    n_states_occupied = property(get_n_states_occupied,set_n_states_occupied)

    def get_n_states_unoccupied(self):
        return self.__n_states_unoccupied
    def set_n_states_unoccupied(self,n):
        self.__n_states_unoccupied = n
        self.__changed = True
    n_states_unoccupied = property(get_n_states_unoccupied,set_n_states_unoccupied)
