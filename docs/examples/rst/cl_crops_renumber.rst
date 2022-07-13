.. _cl_crops_renumber:

Structure file renumbering
--------------------------

With CROPS you can use a simple command to renumber the residues in a structure file so they follow the position in the original sequence.

.. code-block:: shell-session

   crops-renumber 3org.fasta 3org.pdb --output mydir/

The result of the above call is a new content in the output file ``mydir/3org/3org.seqs.pdb`` containing a minimal output of the original ``3org.pdb`` with the residues of all models and chains renumbered according to residue position in the ``3org.fasta`` sequence. Ligands are renumbered with consecutive indices right after the chain ends.

The output directory argument is optional. If not provided, the results will be saved to the sequence file's directory ``3org/3org.seqs.pdb`` by default.

--------------------------------------------------------------

If the residue content and positioning in the sequence and the structure files are not compatible (for instance, gaps in the structure are of different length than the same residue segment in the sequence), the use of the option ``--force_alignment`` or ``-f`` will apply the Needleman-Wunsch algorithm to try to bypass these small disagreements between *.fasta* and *.pdb* sequences.

.. code-block:: shell-session

   crops-renumber 5gup.fasta 5gup.pdb --force_alignment --output mydir/

--------------------------------------------------------------

From a large fasta file and a directory containing several *.pdb* structure file, from which only a few are required, the option ``--preselect`` or ``-p`` allows to preselect as many molecule ids as needed:

.. code-block:: shell-session

   crops-splitseqs PDBall.fasta AllPDBs/ --output mydir/ --preselect 7m6c 4n5b 1o98

This command will create new files only for the three pdb ids inserted, regardless of the number of sequences contained within the input *.fasta* file or the number of structures within the *AllPDBs/* directory.
