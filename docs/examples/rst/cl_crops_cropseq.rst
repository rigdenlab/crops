.. _cl_crops_cropseq:

Sequence file croppping
----------------------

For the removal of residues out of a sequence file containing one or more sequences, this CROPS command will do it for you:

.. code-block:: shell-session

   crops-cropseq 5hea.fasta dbs/pdb_chain_uniprot.csv --output mydir/

This command produces a new sequence file ``mydir/5hea/5hea.crops.to_uniprot.fasta`` containing sequences without the residues cropped as per the interval database ``dbs/pdb_chain_uniprot.csv``. The interval database ``dbs/pdb_chain_uniprot.csv`` in this case is the SIFTS database mapping each residue to a Uniprot reference or none at all (and hence the *to_uniprot* filetag). When a custom interval database is provided instead (the custom ``.csv`` database format must be *pdb_ID*, *monomer_ID*, *integer*, *integer*), the filetag name will be *custom*.

The output directory argument ``--output`` or ``-o`` is optional. If not provided, the results will be saved in the sequence file's directory by default.

--------------------------------------------------------------

From a large fasta file where only a few sequences are required, the option ``--preselect`` or ``-p`` allows to preselect as many molecule ids as needed:

.. code-block:: shell-session

   crops-splitseqs PDBall.fasta --output mydir/ --preselect 7m6c 4n5b 1o98

This command will create new files only for the three pdb ids inserted, regardless of the number of sequences contained within the input *.fasta* file.

--------------------------------------------------------------

As an alternative to the previous one, if the input sequence file contains sequences of different PDB IDs, the user can choose to produce the output with the sequences sorted in different manners by using the ``--sort`` or ``-s`` option followed by one of ``ncrops``, ``percent``, ``ncropsin``, or ``percentin``:

.. code-block:: shell-session

   crops-cropseq 5hea.fasta dbs/pdb_chain_uniprot.csv --sort ncrops --output mydir/

In the above case, the returned list contains the sequences sorted in descendent order by the number of cropped residues, as indicated by argument ``ncrops``. The argument ``percent`` will return the same list sorted by the percentage of residues cropped. Both criteria have a modified version (``ncropsin`` and ``percentin``, respectively). These options will ignore any cropped residues at any end of the chain (terminal segments) when sorting the output sequences.

.. note::

   If the file contains sequences with only one PDB ID, the ``--sort`` option will be ignored.

.. note::

   This option is mutually exclusive with ``--preselect``, and both cannot be used simultaneously.

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

      crops-cropseq 5hea.fasta dbs/pdb_chain_uniprot.csv --terminals --output mydir/

   For instance, in a case in which the intervals imported for one particular chain are ``[5,20]`` and ``[90,125]``, this option will tell CROPS to act as if one single interval ``[5,125]`` is provided, therefore preserving the middle part of the sequence that otherwise would be removed.

2. Sometimes, small contributions from Uniprot sequences other than the main one may not be desired in the cropped version. The option ``--uniprot`` or ``-u`` allows to keep Uniprot residues **only** from those Uniprot references that contribute with a percentage of residues above the given threshold:

   .. code-block:: shell-session

      crops-cropseq 5hea.fasta dbs/pdb_chain_uniprot.csv --uniprot 70 uniclust##_yyyy_mm_consensus --output mydir/

   In the above case, only those Uniprot references that contribute with more than 70% of their original residues are considered.
