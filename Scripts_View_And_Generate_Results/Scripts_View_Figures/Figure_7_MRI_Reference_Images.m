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

N = 256;
c = 240;

load(fullfile(fi1,'MRI_1.mat'),'Im')
It1 = imresize(Im,[N,N]);

load(fullfile(fi1,'MRI_2.mat'),'Im')
It2 = imresize(Im,[N,N]);

load(fullfile(fi1,'MRI_3.mat'),'Im')
It3 = imresize(Im,[N,N]);

load(fullfile(fi1,'MRI_4.mat'),'Im')
It4 = imresize(Im,[N,N]);

load(fullfile(fi1,'MRI_5.mat'),'Im')
It5 = imresize(Im,[N,N]);

load(fullfile(fi1,'MRI_6.mat'),'Im')
It6 = imresize(Im(1:c,1:c),[N,N]);

load(fullfile(fi1,'MRI_7.mat'),'Im')
It7 = imresize(Im,[N,N]);

load(fullfile(fi1,'MRI_8.mat'),'Im')
It8 = imresize(Im,[N,N]);

load(fullfile(fi1,'MRI_9.mat'),'Im')
It9 = imresize(Im(1:c,1:c),[N,N]);

load(fullfile(fi1,'MRI_10.mat'),'Im')
It10 = imresize(Im(1:c,1:c),[N,N]);

xsp = ones(size(It1,1),2);
n1 = 5;
n2 = n1+1;
ysp = ones(2,n1*size(It1,2)+n2*size(xsp,2));

I = [ysp;xsp,It1,xsp,It2,xsp,It3,xsp,It4,xsp,It5,xsp;ysp;xsp,It6,xsp,It7,...
    xsp,It8,xsp,It9,xsp,It10,xsp;ysp];

Fig = figure('Name','Fig.7.MRI Reference Images (I1-I10)','NumberTitle','off');
imshow(I,[],'border','tight')


