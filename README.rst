.. image:: docs/_static/logo_crops_clipart.svg
   :width: 100%
   :align: center

************************************************************************
Cropping and Renumbering Operations for PDB structure and Sequence files
************************************************************************

.. image:: https://readthedocs.org/projects/crops/badge/?version=latest
   :target: http://crops.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://travis-ci.com/rigdenlab/CROPS.svg?branch=master
   :target: https://travis-ci.com/rigdenlab/CROPS
   :alt: CI Status

.. image:: https://codecov.io/gh/rigdenlab/CROPS/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/rigdenlab/CROPS

About CROPS
+++++++++++

CROPS is a useful toolkit developed by the `Rigden <https://github.com/rigdenlab>`_ group at the University of Liverpool that allows to perform a number of customisations on structure and sequence files, aimed at improving the outcomes of a variety of processes. CROPS allows to renumber the index of residues in a PDB structure file in order to match the positional index in the complete sequence. It can also make use of a number of databases to find those residues in the sequence that are either artificially inserted or belong in different biological sequences (e.g. Uniprot sequences), in order to remove artificial insertions that can hamper the search for homologues or decrease the score in a molecular replacement operation. 

This toolkit has been developed to perform three specific tasks:

* renumber residue indices in structure files to make them match residue position in the sequence (provided by a sequence file from the same source);
* crop insertions off sequences in sequence files;
* crop insertions off structures in structure files.

CROPS makes use of a number of databases in order to collect information about segments, insertions, sources, etc. In particular, it is built around the following databases:

* PDB-to-Uniprot `SIFTS database <https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html>`_,
* `Uniclust <https://uniclust.mmseqs.com/>`_ databases,
* custom user-defined interval files.

Supported Formats
+++++++++++++++++

* **Sequences:**  *FASTA* --several header formats recognised. 
* **Structures:** *PDB*. *mmCIF* to be implemented.
* **Databases:**  *CSV*.

Read the Docs !!
++++++++++++++++

You can find all the documentation including practical examples `here <https://crops.readthedocs.io/en/latest/>`_.
