Introduction
=============

.. _GetStarted:

What is bla?
------------


Installation Guide
-------------------
There are different ways to install the RDKit, depending on your needs and environment.
The easiest way is to use the conda package manager.

.. code-block:: bash

   conda install rdkit

If you are a developer or need the latest features, you can build RDKit from source.
There is a blog entry on `Building RDKit from Source <https://greglandrum.github.io/rdkit-blog/posts/2023-03-17-setting-up-a-cxx-dev-env2.html>`_ that provides detailed instructions.

The Main Objects
-----------------
Most workflows start with a mol object, which represents a molecular structure.
There are different ways to construct a molecules. 

.. code-block:: text

   Mol
     └── Conformer
