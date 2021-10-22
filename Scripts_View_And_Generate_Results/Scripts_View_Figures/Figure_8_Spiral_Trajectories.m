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

fi1 = fullfile(fi,'Models_and_Measurements','Measurements_for_Table1'); 

N = 256;
c = 240;

fr = fullfile(fi1,'MRI_8_Spiral_SR_1_Sig_20.mat');
if exist(fr,'file')
    load(fr,'TF','Imn')
else
    error(' MRI_Data folder corrupted ');
end
It11 = 255*imresize(fftshift(TF),[N,N]);
It13 = imresize(abs(Imn(1:240,1:240)),[N,N]);

fr = fullfile(fi1,'MRI_8_Spiral_SR_2_Sig_20.mat');
if exist(fr,'file')
    load(fr,'TF','Imn')
else
    error(' MRI_Data folder corrupted ');
end
It12 = 255*imresize(fftshift(TF),[N,N]);
It14 = imresize(abs(Imn(1:240,1:240)),[N,N]);

xsp = ones(size(It11,1),2);
n1 = 4;
n2 = n1+1;
ysp = ones(2,n1*size(It11,2)+n2*size(xsp,2));

I2 = [ysp;xsp,It11,xsp,It12,xsp,...
    It13,xsp,It14,xsp;ysp];
Fig = figure('Name',...
['Fig. 8. Spiral sampling trajectories M1(10%) and M2(20%),'...
'Zero filled Fourier inversions of undersampled I8 with samples from M1 and M2.']...
,'NumberTitle','off');
imshow(I2,[],'border','tight')
