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

from rasi.base import BasicCalculator

class SchmidtOverlaps(BasicCalculator):
    """
        SchmidtOverlaps calculates harmonic oscillator overlaps using the
        recursion formulas published in P.P. Schmidt Molecular Physics 108 pp.1513-1529 (2010).

        INPUT PARAMETERS:
        ----------------

        omega_occupied .................. Angular oscillation frequency for the occupied state
        omega_unoccupied ................ Angular oscillation frequency for the unoccupied state
        mass ............................ Modal mass
        equilibrium_shift ............... Shift between the equilibrium coordinates of the two states
        n_states_occupied ............... Number of vibrational states to consider in the occupied state
        n_states_unoccupied ............. Number of vibrational states to consider in the unoccupied state


        OUTPUT PARAMETERS:
        -----------------
        
        energies_occupied ............... Array of vibrational energies in the occupied state.
        energies_unoccupied ............. Array of vibrational energies in the unoccupied state.
        overlap_matrix .................. The overlap matrix between the occupied and unoccupied states
                                          (first index corresponds to occupied, second to unoccupied.
        """


    def __init__(self, **kwargs):
        
        self.init_input_variables(
                                  omega_occupied      = None,
                                  omega_unoccupied    = None,
                                  mass                = None,
                                  equilibrium_shift   = None,
                                  n_states_occupied   = 200,
                                  n_states_unoccupied = 200
                                  )
        self.init_output_variables(
                                  energies_occupied   = None,
                                  energies_unoccupied = None,
                                  overlap_matrix      = None
                                  )
        self.set_variables(kwargs)

    def do_update(self):
        from scipy.constants import hbar
        from math import sqrt,log,pi,exp
        from scipy import arange,array

        if self.changed :
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

            self.internal_energies_occupied   =   (arange(float(n_states_occupied))+.5)*hbar*omega_occupied
            self.internal_energies_unoccupied = (arange(float(n_states_unoccupied))+.5)*hbar*omega_unoccupied
            self.internal_overlap_matrix      = array(I)
            return True
        return False

    def partition_function_occupied(self,temperature):
        from harmonic import harmonic_oscillator_partition_function
        return harmonic_oscillator_partition_function(self.omega_occupied,temperature)

    def partition_function_unoccupied(self,temperature):
        from harmonic import harmonic_oscillator_partition_function
        return harmonic_oscillator_partition_function(self.omega_unoccupied,temperature)

        
