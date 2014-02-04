Theory of RASI
==============

The theories implemented in RASI are used to model interactions between the
free carriers in semiconductor devices with point defects.  Extensive research
efforts on the topic of BTI and related phenomena have shown that some of these
defects undergo internal state changes upon carrier capture.
Detailed modeling of this behavior requires to go beyond the 
Shockley-Read-Hall formalism [SRH]_ which is commonly employed in semiconductor
device modeling. Extensions of RASI's defect model compared to SRH include:

1. Detailed modeling of the thermally activated carrier trapping within the
   framework of the *non-radiative multi-phonon transition theory* [NMP]_. This
   essentially leads to an energy-dependent capture cross section that strongly
   varies with the lattice temperature.
2. Multi-state defects. Careful analysis of single-trap capture and emission data
   has established the existence of different internal configurations within
   the point-defects contributing to the recoverable component of the BTI.

Thermally Activated Carrier Trapping
------------------------------------


Bibliography
------------
.. [SRH] Shockley, Read, PUT PAPER REF HERE.
.. [NMP] Some NMP papers
