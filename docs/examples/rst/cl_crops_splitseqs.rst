.. _cl_crops_splitseqs:

Sequence file splitting
----------------------

For the separation of fasta sequences in several files grouped per protein ID just type:

.. code-block:: shell-session

   crops-splitseqs PDBall.fasta --output mydir/

This command produces new sequence files of the format ``mydir/PDBID.fasta`` containing the sequences grouped per protein ID.

The output directory argument ``--output`` or ``-o`` is optional. If not provided, the results will be saved in the sequence file's directory by default.

--------------------------------------------------------------

From a large fasta file where only a few sequences are required, the option ``--preselect`` or ``-p`` allows to preselect as many molecule ids as needed:

.. code-block:: shell-session

   crops-splitseqs PDBall.fasta --output mydir/ --preselect 7m6c 4n5b 1o98

--------------------------------------------------------------

Additionally, the option to separate the sequences by sequence is also available by typing ``--individual`` or ``-i``:

.. code-block:: shell-session

   crops-splitseqs PDBall.fasta --output mydir/ --individual

This command produces new sequence files of the format ``mydir/PDBID_X.fasta`` containing just a single sequence of Protein ID PDBID and (numerical) sequence id X.

Options ``--preselect`` and ``--individual`` can be combined to produce individual sequence files only from the selected molecules.
