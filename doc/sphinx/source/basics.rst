History of RASI
===============

This project was started in late 2013. It is an effort to integrate previously
existing scripts, programs, and code-snippets into a common software framework
which is built as part of the `European Union FP7 project MORDRED
<http://webhotel2.tut.fi/fys/mordred/>`_. It is mainly developed at the 
`Institute for Microelectronics <http://www.iue.tuwien.ac.at>`_ of the 
`Technische Universität Wien <http://www.tuwien.ac.at>`_.

.. Why the Name?
.. -------------

.. RASI is the acronym for *Reliability Analysis of Semiconductor Interfaces*, which
.. describes pretty well what it does. Of course, this acronym was not chosen accidentially.

Basic Implementation Concepts
=============================

The computations performed within RASI combine data from various sources. Due to
their post-processing nature it is possible to represent the calculations performed
by RASI as a hierachy without loops. Every node in this hierachy provides data
to lower-lying elements and receives data from the user or higher-lying elements.
In the implementation, this structure is realized by *Calculator objects*, which
correspond to the calculation steps and which are put together to perform a
specific calculation.

.. graphviz::
    
    digraph {
        c  [label = "Calculator",shape=box];
        p1 [label = "Parameter 1", shape=plaintext];
        p2 [label = "Parameter 2", shape=plaintext];

        c2  [label = "Another\n Calculator",shape=box];
        p21 [label = "Parameter 1", shape=plaintext];
        p22 [label = "Parameter 2", shape=plaintext];

        R [label = "Result", shape=ellipse ] ;
        p1 -> c;
        p2 -> c;
        p21 -> c2;
        p22 -> c2;
        c2 -> c;

        c -> R;
        
    }

Together, these Calculator objects build the mentioned hierachy where the value
in every node of the tree depends on the parameters of the node and
the outputs of the higher-lying nodes. The Calculator objects
implement the necessary logic to ensure that during parameter calibrations
or parameter sweeps only those parts of the calculation tree that have
changed are recalculated.
