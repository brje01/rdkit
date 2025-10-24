Introduction
=============

.. _GetStarted:

What is the RDKit?
-------------------
RDKit stands for Rational Discovery Toolkit.
It is an open-source toolkit for cheminformatics, providing a wide range of functionalities for working with chemical data.

Properties
^^^^^^^^^^
- Business-friendly BSD license
- Core data structures and algorithms in C++
- Python 3.x wrappers generated using Boost.Python
- Java and C# wrappers generated with SWIG
- JavaScript wrappers of most-important functionality

Core Features and Capabilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- 2D and 3D molecular operations
- Descriptor generation for machine learning
- Molecular database cartridge for PostgreSQL
- Cheminformatics nodes for KNIME (distributed from the `KNIME community site <https://www.knime.com/rdkit>`_)


Installation Guide
-------------------
There are different ways to install the RDKit, depending on your needs and environment.
The easiest way is to use the conda package manager.

.. code-block:: bash

   conda install -c conda-forge rdkit 

If you are a developer or need the latest features, you can build RDKit from source.
There is a blog entry on `Building RDKit from Source <https://greglandrum.github.io/rdkit-blog/posts/2023-03-17-setting-up-a-cxx-dev-env2.html>`_ that provides detailed instructions.

