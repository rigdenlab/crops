.. _docs_install:

Installing CROPS
----------------

.. note::

   We are working together with `CCP4 <https://www.ccp4.ac.uk/>`_ to include CROPS in their software toolkit. For now, this is an independent software.

You can install a `Pypi <https://pypi.org/>`_ release of this package using ``pip``:

.. code-block:: shell

   pip install crops

``Conda`` packages are not currently supported, but you can install the Pypi package in a ``conda`` environment following the instructions `here <https://docs.conda.io/projects/conda-build/en/latest/user-guide/tutorials/build-pkgs-skeleton.html>`_.

Installing dependencies
-----------------------

There are two software dependencies:

* **Gemmi:** Learn how to install Gemmi at `gemmi.readthedocs.io <https://gemmi.readthedocs.io/en/latest/install.html>`_,
* **BioPython:** Learn how to install Biopython at `biopython.org <https://biopython.org>`_,

We also recommend that you download and keep updated the two databases that CROPS makes use of:

* from `SIFTS <https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html>`_, the **pdb_chain_uniprot.csv** database;
* from `UniClust <https://uniclust.mmseqs.com/>`_ the most recent realease of the Uniclust (consensus) database in *.fasta* format.
