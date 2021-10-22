v2.0: Sanjay Viswanath, Muthuvel Arigovindan, Imaging Systems Lab, EE, IISc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Codes, MAT Files and Scripts for results presented in
Sanjay Viswanath, Manu Ghulyani, Muthuvel Arigovindan, "Structurally Adaptive 
Multi-Derivative Regularization for Image Recovery from Sparse Fourier Samples"
https://arxiv.org/abs/2105.12775
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inside Scripts_Table_Entries:
Select the Image corresponding to the required result
	Table2: MRI_1 to MRI_10
	
Table2_Full: Generate all results for Table2.
The SNR and SSIM scores of the reconstructions by each method for each image are 
displayed.

Inside <IMAGE> folder:
Select the Sampling Trajectory and Sampling Ratio
	<Trajectory_Name>_SR_<Sampling_Ratio>_Sig_<PSNR>

eq: Radial_SR_1_Sig_20 for 10% Radial trajectory samples with 20dB PSNR.

Select the reconstruction method for the required result
	hcorosa,corosa,tgv2,tv2,hs,Full

eg: Run MRI_1_Radial_SR_1_Sig_20_hcorosa_result to generate HCOROSA reconstruction of 
MRI_1 image from 10% Radial Samples using hcorosa code from folder 'Proposed_Method' 
and stored parameter values from 'Lambda_Table2.mat' in 
'Models_and_Measurements\Measurements_for_Table2' folder. 
The SNR and SSIM scores of the reconstructions are displayed.

    Run MRI_1_Radial_SR_1_Sig_20_Full to generate all (HCOROSA,COROSA,TGV2,TV2,HS) 
reconstructions of MRI_1 image from 10% Radial Samples using codes from folders 'Proposed_Method' and 'Compared_Methods/Regularization_Methods' folders.
The codes use stored parameter values from 'Lambda_Table2.mat' in 
'Models_and_Measurements/Measurements_for_Table2' folder. 
The SNR and SSIM scores of the reconstructions by each method are displayed.