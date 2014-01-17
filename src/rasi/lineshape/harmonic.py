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

def harmonic_oscillator_partition_function(omega,T):
    from scipy.constants import hbar
    from scipy.constants import k as kB
    from math import exp
    return exp(-hbar*omega/(kB*T)*1./2)/(1.-exp(-hbar*omega/(kB*T)))
