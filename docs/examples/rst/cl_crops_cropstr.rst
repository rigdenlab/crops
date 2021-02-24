.. _cl_crops_cropstr:

Structure file cropping
------------------------

For the removal of residues out of a structure file, this CROPS command will do it for you:

.. code-block:: bash

   $> crops-cropstr 3org.fasta 3org.pdb dbs/pdb_chain_uniprot.csv -o mydir/

The result of the above call is a new content in the output file ``3org/3org.crops.to_uniprot.pdb`` containing a minimal output of the original ``3org.pdb`` with the residues of all models and chains renumbered according to residue position in the ``3org.fasta`` sequence. Ligands are renumbered with consecutive indices right after the chain ends. Addional outputs are a renumbered minimal version of the original pdb ``3org/3org.crops.seqs.pdb``, a cropped version with the residue numbers being the position in the old sequence ``3org/3org.crops.to_uniprot.oldids.pdb``, and the cropped sequence ``3org/3org.crops.to_uniprot.fasta``. The interval database ``dbs/pdb_chain_uniprot.csv`` in this case is the SIFTS database mapping each residue to a Uniprot reference or none at all (and hence the *to_uniprot* filetag).

When a custom interval database is provided instead (the custom csv database format must be *pdb_ID*, *monomer_ID*, *integer*, *integer*), the output file name will be ``mydir/5hea/5hea.crops.custom.pdb``.

The output directory argument is optional. If not provided, the results will be saved in the sequence file's directory by default.

.. note::

   The residue content and positioning in the sequence and structure files **must** be compatible with each other (i.e. both files should come from the same source). If that is not the case, an ERROR message will appear.

--------------------------------------------------------------

One of these mutually exclusive conditions can also be imposed:

1. To produce sequences that only discard the *non-Uniprot* (or custom criteria) segments at each of the chains' ends, the option ``--terminal`` or ``-t`` can be added to the command line instruction so only th$

.. code-block:: bash

   $> crops-cropseq 5hea.fasta dbs/pdb_chain_uniprot.csv --terminals --output mydir/

For instance, in a case in which the intervals imported for one particular chain are ``[5,20]`` and ``[90,125]``, this option will tell CROPS to act as if one single interval ``[5,125]`` is provided, therefore $

2. Sometimes, small contributions from Uniprot sequences other than the main one may not be desired in the cropped version. The option ``--uniprot`` or ``-u`` allows to keep Uniprot residues **only** from those$

.. code-block:: bash

   $> crops-cropseq 5hea.fasta dbs/pdb_chain_uniprot.csv --uniprot 70 uniclust##_yyyy_mm_consensus --output mydir/

In the above case, only those Uniprot references that contribute with more than 70% of their original residues are considered.
