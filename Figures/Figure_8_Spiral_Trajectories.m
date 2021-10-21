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
cd Figures

fi1 = fullfile(fi,'MRI_Data','Experiment_1'); 

N = 256;
c = 240;

load(fullfile(fi1,'MRI_8_Spiral_SR_1_Sig_20.mat'),'TF','Imn')
It11 = 255*imresize(fftshift(TF),[N,N]);
It13 = imresize(abs(Imn(1:240,1:240)),[N,N]);

load(fullfile(fi1,'MRI_8_Spiral_SR_2_Sig_20.mat'),'TF','Imn')
It12 = 255*imresize(fftshift(TF),[N,N]);
It14 = imresize(abs(Imn(1:240,1:240)),[N,N]);

xsp = ones(size(It11,1),2);
n1 = 4;
n2 = n1+1;
ysp = ones(2,n1*size(It11,2)+n2*size(xsp,2));

I2 = [ysp;xsp,It11,xsp,It12,xsp,...
    It13,xsp,It14,xsp;ysp];
Fig = figure;
imshow(I2,[],'border','tight')
