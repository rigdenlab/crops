.. _cl_crops_splitseqs:

Sequence file splitting
----------------------

For the separation of fasta sequences in several files grouped per protein ID just type:

.. code-block:: shell-session

   crops-splitseqs PDBall.fasta --output mydir/

This command produces new sequence files of the format ``mydir/PDBID.fasta`` containing the sequences grouped per protein ID.

The output directory argument ``--output`` or ``-o`` is optional. If not provided, the results will be saved in the sequence file's directory by default.

--------------------------------------------------------------

Additionally, the option to separate the sequences by sequence is also available by typing ``--individual`` or ``-i``:

.. code-block:: shell-session

   crops-splitseqs PDBall.fasta --output mydir/ --individual

This command produces new sequence files of the format ``mydir/PDBID_X.fasta`` containing just a single sequence of Protein ID PDBID and (numerical) sequence id X.
