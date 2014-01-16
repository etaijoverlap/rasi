from overlaps.units import *
from overlaps.constants import *
from tools import write_data

def schroedinger_solve(dx,V,m,):
    from scipy        import zeros
    from scipy.linalg import eigh
    # This is a bit weak, interpolation on the finest grid
    # dx = min( xi-xi_old for xi,xi_old in zip(x[1:],x[0:-1]))
    # xmin = x[0]
    # xmax = x[-1]
    # nsteps = int(floor((xmax-xmin)/dx))
    # xint = [ xmin+dx*i for i in xrange(nsteps) ]
    # Vint = [ interpolate(xv,x,V) for xv in xint ]
    # V = Vint
    # x = xint

    nx = len(V)
    H = zeros([nx,nx])

    for i in range(nx):
        if i > 0:
            H[i,i-1] = -hred**2/(2*dx**2*m)
        H[i,i] = hred**2/(dx**2*m)+V[i]
        if i < nx-1:
            H[i,i+1] = -hred**2/(2*dx**2*m)

    Ew,v = eigh(H)

    Ew = Ew.tolist()
    v  = v.transpose().tolist()
    return Ew,v

def schroedinger_solve_bound(dx,pot,m):
    from scipy import array
    (Ew,wave) = schroedinger_solve(dx,pot,m)
    Vmax = max(pot)
    Ew   = [ E for E in Ew if E < Vmax ]
    wave = [ array(w) for E,w in zip(Ew,wave) if E < Vmax ]
    print "%d bound Eigenvalues: min = %f eV, max = %f eV"%(
                              len(wave),min(Ew)/eV,max(Ew)/eV)
    return (Ew,wave)

def parabolic_overlaps( 
                        # Physical parameters
                        k_i,        # Parabolic const. of initial state
                        k_f,        # Parabolic const. of final state
                        m,          # Mass of the system
                        x_s,        # Separation of parabolas
                        # Method parameters
                        dx,         # Spatial sample width
                        x_i_range,  # Distance of boundary to initial minimum
                        x_f_range,  # Distance of boundary to final minimum
                        x_i_0 = 0., # Position of initial minimum
                        n_write=0,    # Number of wavefunctions to write
                        prefix=None,  # Prefix for output
                        echo=True     # Report status to stdout
                        ):
    from overlaps.tools import gen_grid
    from scipy import matrix
    if echo:
        print 
        print "************************************************************************"
        print " Full numeric solution"
        print 
        print 

    # Set up the calculation grid
    x_f_0 = x_i_0+x_s

    boundaries = [ x_i_0 - x_i_range, x_i_0 + x_i_range, 
                   x_f_0 - x_f_range, x_f_0 + x_f_range ]
    x_left  = min(boundaries)
    x_right = max(boundaries)
    v_x = gen_grid(x_left,x_right,dx)

    if echo:
        print "Setting up potentials..."

    pot_i = [ k_i*(x-x_i_0)**2 for x in v_x ]
    pot_f = [ k_f*(x-x_f_0)**2 for x in v_x ]

    if prefix != None:
        write_data(prefix+"_potentials.dat",v_x,pot_i,pot_f)

    if echo:
        print "Solving Schroedinger equation for initial potential..."
    (Ew_i,wave_i) = schroedinger_solve_bound(dx,pot_i,m)

    if echo:
        print "Solving Schroedinger equation for final potential..."
    (Ew_f,wave_f) = schroedinger_solve_bound(dx,pot_f,m)

    if prefix != None:
        if n_write > 0:
            write_data(prefix+"_wf_i.dat",
                              v_x,*tuple(wave_i[:n_write]))
            write_data(prefix+"_wf_f.dat",
                              v_x,*tuple(wave_f[:n_write]))
        write_data(prefix+"_eigs.dat", Ew_i, Ew_f)
    
    if echo:
        print "Calculating overlaps"
    S = []
    for wave_init in wave_i:
        S_i = []
        for wave_final in wave_f:
            S_i.append((wave_init*wave_final).sum()**2)
        S.append(S_i)
    return Ew_i,Ew_f,S

