.. _cl_crops_renumber:

Structure file renumbering
--------------------------

With CROPS you can use a simple command to renumber the residues in a structure file so they follow the position in the original sequence.

.. code-block:: shell-session

   crops-renumber 3org.fasta 3org.pdb -o mydir/

The result of the above call is a new content in the output file ``3org/3org.seqs.pdb`` containing a minimal output of the original ``3org.pdb`` with the residues of all models and chains renumbered according to residue position in the ``3org.fasta`` sequence. Ligands are renumbered with consecutive indices right after the chain ends.

The output directory argument is optional. If not provided, the results will be saved in the sequence file's directory by default.

.. note::

   The residue content and positioning in the sequence and structure files **must** be compatible with each other (i.e. both files should come from the same source). If that is not the case, an ERROR message will appear.
