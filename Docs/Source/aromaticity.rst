.. _Aromaticity:

Aromaticity
===========
Aromaticity is one of those unpleasant topics that is simultaneously simple and impossibly complicated. 
Since neither experimental nor theoretical chemists can agree with each other about a definition, it's necessary to pick something arbitrary and stick to it. 
This is the approach taken in the RDKit.

Instead of using patterns to match known aromatic systems, the aromaticity perception code in the RDKit uses a set of rules. 
The rules are relatively straightforward.

Aromaticity is a property of atoms and bonds in rings. 
An aromatic bond must be between aromatic atoms, but a bond between aromatic atoms does not need to be aromatic.

For example the fusing bonds here are not considered to be aromatic by the RDKit: