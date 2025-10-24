.. _Fingerprints:

Fingerprints
============
The fingerprint of a chemical structures identifies a molecule through a special characterisitic, e.g. by the structure or structural keys. 
They come with several advantages such as the code being highly compact and the representation being in bits. 
Disadvantages are that they can be ambiguous (see hash collisions), are not convertible to other representations, and are dependent on the fragment library used to construct them.


Types of fingerprints
---------------------

For each fingerprint type, RDKit exposes a generator (the unified interface for producing fingerprints). The table below summarizes common fingerprint types and typical use cases.

.. list-table:: Fingerprint types and generators
	 :widths: 20 30 50
	 :header-rows: 1

	 * - fp type
		 - generator function
		 - use case
	 * - Morgan fp
		 - ``GetMorganGenerator()``
		 - similarity searching, especially when you need to understand which atoms contribute to similarity
	 * - RDKit topological fp
		 - ``GetRDKitFPGenerator()``
		 - general-purpose similarity searching and substructure screening
	 * - Atom pair fp
		 - ``GetAtomPairGenerator()``
		 - count-based similarity and when you need fully reversible/explainable fingerprints
	 * - Topological torsion fp
		 - ``GetTopologicalTorsionGenerator()``
		 - count-based similarity with focus on 4-atom sequences

Feature generators and methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each generator class implements the same core methods; however, the atom-pair and topological-torsion generators are count-based by default. "Folding" refers to applying a hashing function to guarantee a fixed-length fingerprint (see https://www.daylight.com/dayhtml/doc/theory/theory.finger.html).

The common generator methods and their characteristics are shown below.

.. csv-table:: Fingerprint generator methods
	 :header: "method","length of returned fp","folding","pros","cons"
	 :quote: "

	 "GetFingerprint()","fixed size (size of fp)","yes","fast similarity calculations; memory efficient; works well with similarity search algorithms","hash collisions"
	 "GetSparseFingerprint()","not fixed length","no","memory efficient when most bits are zero; better for very large feature spaces","variable-size outputs complicate downstream tasks"
	 "GetCountFingerprint()","fixed size (size of fp)","yes","preserves count information — better for count-based similarity metrics","larger memory requirements; hash collisions"
	 "GetSparseCountFingerprint()","not fixed length","no","full count information without hashing/folding","slower similarity calculations; variable-size outputs can cause problems for downstream tasks"
