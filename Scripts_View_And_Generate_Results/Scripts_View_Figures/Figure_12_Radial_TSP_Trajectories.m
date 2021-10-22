%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sanjay Viswanath, Manu Ghulyani, Muthuvel Arigovindan,"Structurally Adaptive
% Multi-Derivative Regularization for Image Recovery from Sparse Fourier Samples"
% https://arxiv.org/abs/2105.12775
% v2.0: Sanjay Viswanath, Muthuvel Arigovindan, Imaging Systems Lab, EE, IISc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

cd ..
cd ..
fi = pwd;
cd Scripts_View_And_Generate_Results
cd Scripts_View_Figures

fi1 = fullfile(fi,'Models_and_Measurements','Reference_Models'); 
fi2 = fullfile(fi,'Models_and_Measurements','Measurements_for_Table2'); 

load(fullfile(fi1,'MRI_8.mat'),'Im')
It1 = Im;

load(fullfile(fi2,'MRI_8_Radial_SR_1_Sig_20.mat'),'TF','Imn')
It2 = 255*fftshift(TF);
It3 = abs(Imn);

load(fullfile(fi2,'MRI_8_TSP_SR_1_Sig_20.mat'),'TF','Imn')
It4 = 255*fftshift(TF);
It5 = abs(Imn);


xsp = ones(size(It1,1),10);
n = 4;
ysp = ones(2,n*size(It2,2)+(n+1)*size(xsp,2));
I = [ysp;xsp,It2,xsp,It3,xsp,It4,xsp,It5,xsp; ysp];
Fig1 = figure('Name',...
['Fig.12. Radial sampling trajectory R1 and TSP trajectory T1 along with '...
'Zero filled Fourier inversions of undersampled I8 image with samples from corresponding '...
'trajectories.']...
,'NumberTitle','off');
imshow(I,[],'border','tight')

