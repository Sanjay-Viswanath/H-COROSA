v2.0: Sanjay Viswanath, Muthuvel Arigovindan, Imaging Systems Lab, EE, IISc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Codes, MAT Files and Scripts for results presented in
Sanjay Viswanath, Manu Ghulyani, Muthuvel Arigovindan, "Structurally Adaptive 
Multi-Derivative Regularization for Image Recovery from Sparse Fourier Samples"
https://arxiv.org/abs/2105.12775
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DAGAN:
        Trained DAGAN networks for reconstruction with 10% and 20% Random trajectories 
used in Table3.
        Because of GitHub restrictions on file-size, the networks have been compressed using 
7-zip software and need to be extracted. The details are present in the inner folders.

PNP:
	Pretrained PNP networks used for reconstructions in Table3 
(RealSN_DnCNN_noise5) and Table4 (RealSN_DnCNN_noise5).

Note: For using these networks in generating results with the datasets used in paper,
file and folder names have to be set in the author provided codes for

DAGAN: https://github.com/tensorlayer/DAGAN

PNP: https://github.com/uclaopt/Provable_Plug_and_Play

These are public codes and redistribution is avoided.

For other datasets, the trained networks provided here may be used with discretion 
regarding the suitability of their training conditions.