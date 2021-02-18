.. _docs_install:

Installing CROPS
----------------

We are working together with CCP4 to include CROPS in their software toolkit. For now, this is an independent software. 

1. Clone the repository
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

    git clone https://github.com/rigdenlab/crops

2. Install dependencies
^^^^^^^^^^^^^^^^^^^^^^^

The only dependency required is GEMMI. Follow `this link <https://gemmi.readthedocs.io/en/latest/install.html>`_ for installation instructions.

3. Build CROPS
^^^^^^^^^^^^^^

After all the requirements are met, build CROPS and it will be ready for use.

.. code-block:: shell

    mkdir crops/build
    cd crops/build
    cmake ..
    make install


