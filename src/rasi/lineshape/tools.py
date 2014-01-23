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

class LineShapeBuffer(BasicCalculator):
    """
        Samples a lineshape on a given grid. May be used to significantly speed
        up the calculation of numerical integrals.
    """
    def __init__(self, **kwargs):
        self.init_input_variables(
                                  lineshape   = None,
                                  energy_grid = None
                                  )
        self.init_output_variables(
                                  oxidation_values = None,
                                  reduction_values = None,
                                  oxidation        = None,
                                  reduction        = None
                                  )
        self.set_variables(kwargs)

    def do_update(self):
        from scipy.interpolate import interp1d
        ls_changed = self.lineshape.update()

        if self.changed or ls_changed:
            self.internal_oxidation_values = self.lineshape.oxidation(self.energy_grid)
            self.internal_reduction_values = self.lineshape.reduction(self.energy_grid)

            self.internal_oxidation = interp1d(self.energy_grid,self.oxidation_values,bounds_error=False,fill_value=0.)  
            self.internal_reduction = interp1d(self.energy_grid,self.reduction_values,bounds_error=False,fill_value=0.)  
            return True
        return False


class LineShapeEnergyCorrector(BasicCalculator):
    """
        Shifts the energy scale of the line shape by a given amount.
    """
    def __init__(self, **kwargs):
        self.init_input_variables(
                                  lineshape         = None,
                                  correction_energy = 0.
                                 )
        self.init_output_variables()
        self.set_variables(kwargs)

    def do_update(self):
        return self.changed or self.lineshape.update()

    def oxidation(self,E):
        E_c = self.correction_energy
        return self.lineshape.oxidation(E + E_c)

    def reduction(self,E):
        E_c = self.correction_energy
        return self.lineshape.reduction(E + E_c)

def calc_modal_vector(atoms1,atoms2):
    """
        Calculate the 'modal vector', i.e. the difference vector between the two configurations.
        The minimum image convention is applied!
    """
    from scipy.linalg import inv
    from scipy        import array,dot
    from scipy        import sign,floor
    cell1 = atoms1.get_cell()
    cell2 = atoms2.get_cell()

    # The cells need to be the same (otherwise the whole process won't make sense)
    if (cell1 != cell2).any():
        raise ValueError("Encountered different cells in atoms1 and atoms2. Those need to be the same.")
    cell = cell1

    icell = inv(cell)
                                            
    frac1 = atoms1.get_scaled_positions()
    frac2 = atoms2.get_scaled_positions()
    modal_vector_frac = frac1 - frac2
    for i in range(modal_vector_frac.shape[0]):
        for j in range(modal_vector_frac.shape[1]):
            if abs(modal_vector_frac[i,j]) > .5:
                value = modal_vector_frac[i,j]
                vsign = sign(modal_vector_frac[i,j])
                absvalue = abs(value)
                modal_vector_frac[i,j] = value - vsign*floor(absvalue+.5)
    return dot(modal_vector_frac,cell)



############# Everything below this line may be legacy code. ####################
"""
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

""" 
