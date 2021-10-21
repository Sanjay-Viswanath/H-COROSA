v2.0: Sanjay, Imaging Systems Lab, Dept. of EE, IISc, Bangalore

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Codes, MAT Files and Scripts for results presented in
Sanjay Viswanath, Manu Ghulyani, Muthuvel Arigovindan, "Structurally Adaptive 
Multi-Derivative Regularization for Image Recovery from Sparse Fourier Samples"
https://arxiv.org/abs/2105.12775
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inside Scripts_Table_Entries:
Select the Image corresponding to the required result
	
	Table4: MRI_1,MRI_2,MRI_3,MRI_8,MRI_9,MRI_10

Table4_Full: Generate all results for Table4.
The SNR and SSIM scores of the reconstructions for each image are displayed.

Inside <IMAGE> folder:
Select the Sampling Trajectory and Noise Variance
	<Trajectory_Name>_SR_<Sampling_Ratio>_var_<Noise_Variance>

eg: Spiral_SR_1_var_15 for 10% Spiral trajectory samples with Noise variance of 15.

Select the reconstruction method for the required result
	hcorosa

eg: Run MRI_1_Spiral_SR_1_var_15_hcorosa_result to generate HCOROSA reconstruction of 
MRI_1 image from 10% Spiral Samples using hcorosa code from folder 'Proposed_Method' 
and stored parameter values from 'Lambda_Table4.mat' in 
'Models_and_Measurements\Measurements_for_Table4' folder. 
The SNR and SSIM scores of the reconstructions are displayed.