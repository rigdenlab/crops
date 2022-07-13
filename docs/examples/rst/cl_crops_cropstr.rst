.. _cl_crops_cropstr:

Structure file cropping
------------------------

For the removal of residues out of a structure file, this CROPS command will do it for you:

.. code-block:: shell-session

   crops-cropstr 3org.fasta 3org.pdb dbs/pdb_chain_uniprot.csv --output mydir/

The result of the above call is the output file ``3org/3org.crops.to_uniprot.pdb`` containing a minimal output of the original ``3org.pdb`` with the residues of all models and chains cropped and renumbered according to residue position in the new sequences produced in ``3org.crops.to_uniprot.fasta``. Ligands are renumbered with consecutive indices right after the chain ends. Addional outputs are a renumbered minimal version of the original pdb ``3org/3org.crops.seqs.pdb`` (same as returned by ``crops-renumber``), and a cropped version with the residue numbers being the position in the **old** sequence ``3org/3org.crops.to_uniprot.oldids.pdb``. The interval database ``dbs/pdb_chain_uniprot.csv`` in this case is the SIFTS database mapping each residue to a Uniprot reference or none at all (and hence the *to_uniprot* filetag). When a custom interval database is provided (the custom ``.csv`` database format must be *pdb_ID*, *monomer_ID*, *integer*, *integer*), the filetag name will be *custom* instead.

The output directory argument is optional. If not provided, the results will be saved in the sequence file's directory by default.

.. note::

   The residue content and positioning in the sequence and structure files **must** be compatible with each other (i.e. both files should come from the same source). If that is not the case, an ERROR message will appear.

--------------------------------------------------------------

From a large fasta file and a directory containing several *.pdb* structure file, from which only a few are required, the option ``--preselect`` or ``-p`` allows to preselect as many molecule ids as needed:

.. code-block:: shell-session

   crops-splitseqs PDBall.fasta AllPDBs/ --output mydir/ --preselect 7m6c 4n5b 1o98

This command will create new files only for the three pdb ids inserted, regardless of the number of sequences contained within the input *.fasta* file or the number of structures within the *AllPDBs/* directory.

--------------------------------------------------------------

Additionally, the option to separate the sequence files by unique sequence is also available by typing ``--individual`` or ``-i``:

.. code-block:: shell-session

   crops-splitseqs 3org.fasta 3org.pdb dbs/pdb_chain_uniprot.csv --output mydir/ --individual

This command produces new sequence files of the format ``mydir/PDBID_X.fasta`` containing just a single sequence of Protein ID PDBID and (numerical) sequence id X.

Options ``--preselect`` and ``--individual`` can be combined to produce individual sequence files only from the selected molecules.

--------------------------------------------------------------

Additionally, one of these mutually exclusive conditions can also be imposed:

1. To produce sequences that only discard the *non-Uniprot* (or custom criteria) segments at each of the chains' ends, the option ``--terminals`` or ``-t`` can be added to the command line instruction so only the unwanted parts at the ends are removed:

   .. code-block:: shell-session

      crops-cropstr 3org.fasta 3org.pdb dbs/pdb_chain_uniprot.csv --terminals --output mydir/

   For instance, in a case in which the intervals imported for one particular chain are ``[5,20]`` and ``[90,125]``, this option will tell CROPS to act as if one single interval ``[5,125]`` is provided, therefore preserving the middle part of the sequence and structure that otherwise would be removed.

2. Sometimes, small contributions from Uniprot sequences other than the main one may not be desired in the cropped version. The option ``--uniprot`` or ``-u`` allows to keep Uniprot residues **only** from those Uniprot references that contribute with a percentage of residues above the given threshold:

   .. code-block:: shell-session

      crops-cropstr 3org.fasta 3org.pdb dbs/pdb_chain_uniprot.csv --uniprot 70 uniclust##_yyyy_mm_consensus --terminals --output mydir/

   In the above case, only those Uniprot references that contribute with more than 70% of their original residues are considered.
