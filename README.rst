.. image:: docs/_static/logo_crops_clipart.svg
   :width: 100%
   :align: center

************************************************************************
Cropping and Renumbering Operations for PDB structure and Sequence files
************************************************************************

.. image:: https://readthedocs.org/projects/crops/badge/?version=latest
   :target: http://crops.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

About CROPS
+++++++++++

CROPS is a useful toolkit developed by the `Rigden group <https://github.com/rigdenlab>`_ at the University of Liverpool, England, that allows to perform a number of customisations on structure and sequence files, aimed at improving the outcomes of a variety of processes. CROPS allows to renumber the index of residues in a PDB structure file in order to match the positional index in the complete sequence and extract specific sequences from a large sequence file. It can also make use of a number of databases to find those residues in the sequence that are either artificially inserted or belong in different biological sequences (e.g. Uniprot sequences), in order to remove artificial insertions that can hamper the search for homologues or decrease the score in a molecular replacement operation.

This toolkit has been developed to perform three specific tasks:

* renumber residue indices in structure files to make them match residue position in the sequence (provided by a sequence file from the same source);
* crop insertions off sequences in sequence files;
* crop insertions off structures in structure files;
* split large sequence files to extract specific sequences from them.

CROPS makes use of a number of databases in order to collect information about segments, insertions, sources, etc. In particular, it is built around the following databases:

* PDB-to-Uniprot `SIFTS database <https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html>`_,
* `Uniclust <https://uniclust.mmseqs.com/>`_ databases,
* custom user-defined interval files.

If your software has any dependencies to versions older than CROPS v0.5.0, please update your software. These versions are now deprecated.

Supported Formats
+++++++++++++++++

* **Sequences:**  *FASTA* --several header formats recognised.
* **Structures:** *PDB*. *mmCIF* to be implemented.
* **Databases:**  *CSV*.

Requirements
++++++++++++

* **Gemmi:** `gemmi.readthedocs.io <https://gemmi.readthedocs.io>`_,
* **BioPython:** `biopython.org <https://biopython.org>`_.

Read the Docs !!
++++++++++++++++

You can find all the documentation including installation instructions and practical examples at `crops.readthedocs.io <https://crops.readthedocs.io/en/latest/>`_.
