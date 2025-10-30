.. _Main Objects:

The Main Objects
================
Most workflows start with a mol object, which represents a molecular structure.
There are different ways to construct a molecule. 

.. code-block:: text

   Mol
   ├── Atom
   ├── Bond
   └── Conformer (atomic coordinates)

Every Way To Construct a Molecule Object (``Mol``)
--------------------------------------------------
A mol object can be constructured in various ways:

- From a SMILES string
- From a SMARTS pattern
- From a supported file format (Mol2, SDF, PDB, TPL, Marvon, XYZ)
- From a MolBlock (a string representation of a mol file)
- From specialized formats like amino acid/nucleotide sequences, FASTA, RDKit-generated SVG with embedded metadata, PNG images with embedded molecule metadata, SCSR

.. code-block:: python

   from rdkit import Chem
   # mol from SMILES
   mol = Chem.MolFromSmiles('CCO')
   # mol from 3D coordinates
   mol = Chem.MolFromMolFile('molecule.sdf')


Common parameters on the mol object construction are ``sanitize`` and ``removeHs``.

- ``sanitize``: If set to True (default), the molecule will be sanitized upon creation, checking for valence issues, aromaticity, etc.
- ``removeHs``: If set to True (default), all explicit hydrogen atoms will be removed from the molecule upon creation.

Molecules are composed of ``Atom`` and ``Bond`` objects representing the molecular graph structure. 
They can only be accessed through the mol object methods ``GetAtoms()`` and ``GetBonds()``.
See :ref:`Mol object <Mol object>` for more details on the mol object and its methods.

See :ref:`Atom object <Atom object>` for more details on the atom object and its methods.

See :ref:`Bond object <Bond object>` for more details on the bond object and its methods.

Reading In Multiple Molecules
------------------------------


Working with Conformers 
------------------------
The ``Conformer`` object represents the spatial coordinates (positions in 3D space) of the atoms in a molecule.
A mol object can have multiple conformers, which can be accessed through the mol object methods ``GetConformer()`` and ``GetConformers()``.
The ``Conformer`` object only stores postionsal information indexed by an atom ID.
It does not have information about atoms or bonds themselves.
:ref:`Conformer object <Conformer object>` for more details on the conformer object and its methods.

