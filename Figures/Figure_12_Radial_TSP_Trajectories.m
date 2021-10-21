%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sanjay Viswanath, Manu Ghulyani, Muthuvel Arigovindan,"Structurally Adaptive
% Multi-Derivative Regularization for Image Recovery from Sparse Fourier Samples"
% https://arxiv.org/abs/2105.12775
% v1.0: Sanjay Viswanath, ISL, Dept. of EE, IISc, Bangalore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

cd ..
fi = pwd;
cd Figures;
fi1 = fullfile(fi,'MRI_Data','Reference_Images');
fi2 = fullfile(fi,'MRI_Data','Experiment_2');

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
Fig1 = figure;
imshow(I,[],'border','tight')

