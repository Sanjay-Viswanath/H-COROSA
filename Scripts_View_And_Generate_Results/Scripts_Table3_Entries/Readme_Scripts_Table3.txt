v2.0: Sanjay, Imaging Systems Lab, Dept. of EE, IISc, Bangalore

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Codes, MAT Files and Scripts for results presented in
Sanjay Viswanath, Manu Ghulyani, Muthuvel Arigovindan, "Structurally Adaptive 
Multi-Derivative Regularization for Image Recovery from Sparse Fourier Samples"
https://arxiv.org/abs/2105.12775
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inside Scripts_Table_Entries:
Select the Image corresponding to the required result
	
	Table3: T1 to T6

Table3_Full: Generate all results for Table3.
The SNR and SSIM scores of the reconstructions by each method for each image are 
displayed.

Inside <IMAGE> folder:
Select the Sampling Ratio (Only Random Trajectory in Table3)
	<Image_Name>_<Sampling_Ratio>per

eg: T1_10per for 10% Random trajectory samples of T1 image.

Select the reconstruction method for the required result
	hcorosa,corosa,Full

eg: Run T1_10per_hcorosa_result to generate HCOROSA reconstruction of T1 image from 10% 
Random Samples using hcorosa code from folder 'Proposed_Method' and stored parameter 
values from 'Lambda_Table3.mat' in 'Models_and_Measurements\Measurements_for_Table3' 
folder. 
The SNR and SSIM scores of the reconstructions are displayed.

    Run T1_10per_Full to generate all (HCOROSA,COROSA) reconstructions of T1 image from 
10% Random Samples using codes from folders 'Proposed_Method' and 'Compared_Methods/Regularization_Methods' folders.
The codes use stored parameter values from 'Lambda_Table3.mat' in 
'Models_and_Measurements/Measurements_for_Table3' folder. 
The SNR and SSIM scores of the reconstructions by each method are displayed.