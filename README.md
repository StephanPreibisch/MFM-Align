# MFM-Align
Software for the alignement of multi-focal microscopy data (Descriptor-based in XY, Information Entropy based in Z)

This software allows the registration of multi-focal microscopy data, which captures several z-planes simultaneously. To achieve aligment we split registration in XY and Z. In XY, we solve a RigidModel2d for all planes over all channels aligning all tiles and interpolating them so they match. In Z, we only determine the correct location in Z, no image interpolation is done. To achieve this, we compute an autofocus-like function over the entire Z-Range of each microscopy plane. Matching these functions using a 1d-Translation Model yields the correct positions in Z. For further details please read the publication on JCB: --Link we available upon final publication --

This software is written in JAVA. The project is Mavenized to allow simple import and running. We suggest cloning this repository into your Eclipse workspace, and then select in Eclipse File>Import>Existing Maven Project ... It will solve all dependency issues and download all requied libraries automatically.

A plugin version of this software is available as Fiji plugin and is hosted on github: https://github.com/bigdataviewer/Descriptor_based_registration - The Fiji-Wiki entry for this plugin can be found here: http://fiji.sc/Descriptor-based_registration_%282d/3d%29
