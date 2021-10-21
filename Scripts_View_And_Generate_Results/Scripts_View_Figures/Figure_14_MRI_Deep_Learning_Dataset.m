%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sanjay Viswanath, Manu Ghulyani, Muthuvel Arigovindan,"Structurally Adaptive
% Multi-Derivative Regularization for Image Recovery from Sparse Fourier Samples"
% https://arxiv.org/abs/2105.12775
% v2.0: Sanjay Viswanath, ISL, Dept. of EE, IISc, Bangalore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

cd ..
cd ..
fi = pwd;
cd Scripts_View_And_Generate_Results
cd Scripts_View_Figures
fi2 = fullfile(fi,'Models_and_Measurements','Measurements_for_Table3');

N = 256;

load(fullfile(fi2,'T1_10per.mat'),'Im')
It1 = Im;

load(fullfile(fi2,'T2_10per.mat'),'Im')
It2 = Im;

load(fullfile(fi2,'T3_10per.mat'),'Im')
It3 = Im;

load(fullfile(fi2,'T4_10per.mat'),'Im')
It4 = Im;

load(fullfile(fi2,'T5_10per.mat'),'Im')
It5 = Im;

load(fullfile(fi2,'T6_10per.mat'),'Im')
It6 = Im;

x1 = 35;
x2 = 220;
y1 = 35;
y2 = 220;
Nx = x2-x1+1;
Ny = y2-y1+1;

load(fullfile(fi2,'T1_10per.mat'),'TF')
It7 = 255*imresize(fftshift(TF),[Nx,Ny]);

load(fullfile(fi2,'T2_20per.mat'),'TF')
It8 = 255*imresize(fftshift(TF),[Nx,Ny]);

xsp = ones(size(It1(x1:x2,y1:y2),1),2);
n = 4;
ysp = ones(2,n*size(It1(x1:x2,y1:y2),2)+(n+1)*size(xsp,2));
I = [ysp;xsp,It1(x1:x2,y1:y2),xsp,It2(x1:x2,y1:y2),xsp,...
    It3(x1:x2,y1:y2),xsp,It4(x1:x2,y1:y2),xsp;ysp;xsp,...
    It5(x1:x2,y1:y2),xsp,It6(x1:x2,y1:y2),xsp,...
    It7,xsp,It8,xsp;ysp];
Fig1 = figure('Name',['Fig. 14. MRI Deep Learning Comparison Dataset: T1-T6,',...
    ' Random masks R1 (10%) and R2 (20%)'],'NumberTitle','off');
imshow(I,[],'border','tight')

