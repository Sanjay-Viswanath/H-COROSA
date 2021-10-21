v2.0: Sanjay, Imaging Systems Lab, Dept. of EE, IISc, Bangalore

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Codes, MAT Files and Scripts for results presented in
Sanjay Viswanath, Manu Ghulyani, Muthuvel Arigovindan, "Structurally Adaptive 
Multi-Derivative Regularization for Image Recovery from Sparse Fourier Samples"
https://arxiv.org/abs/2105.12775
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Stored Reconstruction Results corresponding to each Table

Reconstructions_in_Table1: Reconstructions of MRI 1 to MRI 10 images from 10% Spiral
(filename-tag: Spiral_SR_1) and 20% Spiral samples (filename-tag: Spiral_SR_2).

Reconstructions_in_Table2: Reconstructions of MRI 1 to MRI 10 images from 12% TSP
(filename-tag: TSP_SR_1) and 10% Radial samples (filename-tag: Radial_SR_1).

Reconstructions_in_Table3: Reconstructions of T1 to T6 images from 10% Random
(filename-tag: Random_SR_1) and 20% Random samples (filename-tag: Random_SR_2).

Reconstructions_in_Table4: Reconstructions of MRI 1,2,3,8,9,10 images from 10% 
Spiral (filename-tag: Spiral_SR_1), 20% Spiral (filename-tag: Spiral_SR_2), 
12% TSP (filename-tag: TSP_SR_1) and 10% Radial samples (filename-tag: Radial_SR_1).

For each Table, the file naming convention is

Tables 1-3
<Image>_<Trajectory>_SR_<Sampling_Ratio>_Sig_<PSNR>_<Reconstruction_Method>_best_<Tuned
_for_SSIM/SNR>

Table 4
<Image>_<Trajectory>_SR_<Sampling_Ratio>_var_<Noise_Variance>_<Reconstruction_Method>_b
est_<Tuned_for_SSIM/SNR>

Possible values for each field:

<Image>: MRI 1-10 or T1-T6
<Trajectory>: Spiral, TSP, Radial, Random
<Sampling_Ratio>: 1 for 10%/12% and 2 for 20%
<PSNR>: 20 for 20dB
<Noise_Variance>: 15 
<Reconstruction_Method>: hcorosa(proposed method), corosa (COROSA), tgv2 (Second 
Order TGV), tv2 (Second Order TV), hs (HS), DAGAN, PNP
<Tuned_for_SSIM/SNR>: SSIM for SSIM scores, SNR for SNR scores given in paper.