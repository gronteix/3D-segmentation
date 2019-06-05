# 3D-segmentation

Short repository to temporarily store my spheroid segmentation algorithm.

## Context

Image analysis algorithm to find cell nuclei in confocal image stacks. These are then used to defins a connection grph of the spheroid.

## Use and file prep

The details of the image analysis are customized for our images. The same is true of the file organization. Ours has the following structure:

 -> main folder
 ---> well folder
 -----> time folder
 
 Thus, there are potentially several wells and several time points for each well to analyze. The data is stored accordingly.
