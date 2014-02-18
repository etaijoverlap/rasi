.. RASI documentation master file, created by
   sphinx-quickstart on Sun Feb  2 22:46:30 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Reliability Analysis of Semiconductor Interfaces (RASI) documentation
=====================================================================

RASI is a software package for the modeling of defects in semiconductors and
(especially) at semiconductor-insulator interfaces. It is implemented as a
library for the `Python programming language <http://www.python.org>`_ within
the `European Union FP7 project MORDRED <http://webhotel2.tut.fi/fys/mordred/>`_.

As of yet, the library implements a post-processing toolflow that brings together 
microscopic defect data from electronic structure calculations and charge
carrier data from macroscopic device models. 

.. graphviz::

    digraph foo {
        es[label="Defect Model\n(Electronic Structure Level)",shape=box];
        model [label="Model Assumptions",shape=box];
        device [label="Charge Carrier Model\n(Macroscopic Device Simultaion)",shape=box];
        meas [label = "Predicted Measurement Data",shape=box];
         es      -> model;
         device  -> model;
         model   -> meas;
    }


Contents:

.. toctree::
   basics
   theory
   implementation
   :maxdepth: 2



.. Indices and tables
   ==================
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`

