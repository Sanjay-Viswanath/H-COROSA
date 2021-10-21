v2.0: Sanjay, Imaging Systems Lab, Dept. of EE, IISc, Bangalore

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Codes, MAT Files and Scripts for results presented in
Sanjay Viswanath, Manu Ghulyani, Muthuvel Arigovindan, "Structurally Adaptive 
Multi-Derivative Regularization for Image Recovery from Sparse Fourier Samples"
https://arxiv.org/abs/2105.12775
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How to see final results directly
Tables: Run the scripts from folder 'Scripts_View_Tables'
	eg: View_Table1 to see Table1 as an array
Figures: Run the scripts from folder 'Scripts_View_Figures'
	eg: Figure_10_MRI_10_Spiral_SR_2_Sig_20_Selected_Region to see Figure 10 in 
paper ( Results from reconstruction of MRI 10 from 20% Spiral samples ).

Note that these scripts read stored reconstruction from the MAT_Files in 
'Generated_Reconstruction_Results' folder.

To perform reconstructions again while running these scripts, press y at the initial prompt. 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How to see individual results: 
Inside Scripts_Generate_Results folder:

Select the Table corresponding to the required result
	Scripts_Table1_Entries: Comparison between reconstructions with 10% and 20% 
Spiral trajectories
	Scripts_Table2_Entries: Comparison between reconstructions with 12% TSP and 
10% Radial trajectories
	Scripts_Table3_Entries: Comparison between reconstructions with 10% and 20% 
Random trajectories
	Scripts_Table4_Entries: Comparison between HCOROSA and PNP method with four 
trajectories: Spiral(10% and 20%),TSP(12%),Radial(10%).

For further details including file naming conventions, please see the detailed Readme 
files in these folders.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%