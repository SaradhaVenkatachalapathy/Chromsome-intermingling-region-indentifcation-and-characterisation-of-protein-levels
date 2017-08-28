# Chromsome-intermingling-region-indentifcation-and-characterisation-of-protein-levels

WRITTEN BY: SARADHA VENKATACHALAPATHY

DATE LAST EDITED: April, 2017

Executed in FIJI (ImageJ v 1.51n). Requires Bio-Formats

DESCRIPTION: This script takes a 3D image stack(x,y,z) with four channels: Nucleus, chromosomes (A and B) and Protein, and identifies the overlapping regions between chromosome A and B in 3D. It also measures the amount of protein levels in these overlapping regions. 

ASSUMPTIONS: The input image is a confocal z-stack (slice height=0.5microns) with 4 channels, the channel One contains the nucleus and then there are two channels for the chromosomes and one channel for RNA polymerase2. The channel numbers for the chromosomes and Pol2 are specified. 

ALGORITHM: Two dialog box opens where the user inputs the source (where the images are stored) the results (where results need to be saved) directory. The program opens nucleus channel, uses the Gaussian filter (sigma=0.15, scaled) to the stack and the thresholds the image and asks the user to check the thresholded image (optional). The user can then adjust the threshold for precision. The program then converts the thresholded image to a binary image and these images are stored in the results folder. For the chromosome channels, it identifies nuclear region and then it performs the aforementioned processes. The threshold values are stored in the datafiles directory. The binary chromosome images are opened in combination (2-3) and passed through the AND filter to identify the intermingling regions (IMR) and the images saved in the results directory. Next, the volumes of the binary nucleus, chromosomes and IMRs are measured and saved in the datafiles directory. Then, the amount of pol2 in the total image, nucleus, and the IMR are calculated and stored in the datafiles directory.


