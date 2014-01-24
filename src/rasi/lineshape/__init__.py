"""
    Reliability Analysis of Semiconductor Interfaces -- Lineshape collection
    ------------------------------------------------------------------------

    Lineshapes are used to model the influence of the vibrational motion of
    a system on its electronic transitions. There are plenty of different
    approaches for calculating those functions and some of them can be found
    in this collection.
    As this library was written with trapping phenomena in mind, the basic
    nomenclature used is based on the following idea:

    Picture a point-defect somewhere in your material. This defect offers
    localized states to electrons (and these electrons are thought of as
    "real" electrons, as opposed to the quasiparticles that are usually
    the central actors of semiconductor device simulation). Each of these
    states can either be occupied, or unoccupied by an electron. The
    lineshapes, together with the electronic matrix elemnts, then determine
    the transition between the occupied and the unoccupied state, where
    the process of going from an occupied to an unoccupied state is referred
    to as "oxidation" and the reverse process is referred to as "reduction".
    
"""
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


from tools     import calc_modal_vector, \
                      LineShapeBuffer, \
                      LineShapeEnergyCorrector

from overlaps import DiscreteLineShape, \
                     SmearedLineShape

from classical import ClassicalLineShape
from schmidt   import SchmidtOverlaps
