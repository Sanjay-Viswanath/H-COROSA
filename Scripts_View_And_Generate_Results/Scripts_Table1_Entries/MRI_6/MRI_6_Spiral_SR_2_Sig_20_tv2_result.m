
cd .. 
cd .. 
cd .. 
close all; 
clear all; 

% Setting Folder Paths 
addpath(genpath('Compared_Methods')); 
fold_in = fullfile('Models_and_Measurements','Measurements_for_Table1'); 
fold_ref = fullfile('Models_and_Measurements','Reference_Models'); 
fold_out = fullfile('Generated_Reconstruction_Results','Reconstructions_in_Table1'); 

% Load Reference Image 
load(fullfile(fold_ref,'MRI_6.mat'),'Im'); 
% Load Measurement 
load(fullfile(fold_in,'MRI_6_Spiral_SR_2_Sig_20.mat'),'Imn','TF');
disp('MRI_6_Spiral_SR_2_Sig_20')

tic;

% TV2 Reconstruction 
disp('TV2 Reconstruction')
% Load Lambda 
load(fullfile(fold_in,'Lambda_Table1.mat'),'Lambda_Table1');
ltv2_snr = Lambda_Table1(29,4); 
ltv2_ssim = Lambda_Table1(29,3); 
NIter = 250; % Number of TV2 iterations 
betaadmm = 1; % Parameter beta in ADMM 
UB = 300; % Upperbound for image values, Effective bounding constraint is [0,UB] 
noisetype = 1; % 1 for AWGN 
lam = ltv2_snr; % Lambda as tuning paramter for regularization, Tuned for SNR 
[Irec, ~] = tv2(TF, Imn, lam, betaadmm,  ...
    NIter, UB, 2, noisetype);
save(fullfile(fold_out,'MRI_6_Spiral_SR_2_Sig_20_tv2_best_snr.mat'),'Irec','lam'); 
% Compute SNR 
temp1 = Im(1:240,1:240); 
temp2 = Irec(1:240,1:240); 
sc1 = sum(temp1(:).*temp2(:))./sum(temp2(:).^2); 
SNR_Score = round(10*log10(sum(temp1(:).^2)/sum((sc1*temp2(:) - temp1(:)).^2)),2); 
Irec_SNR = Irec; 
if ltv2_snr ~= ltv2_ssim 
lam = ltv2_ssim; % Lambda as tuning paramter for regularization, Tuned for SSIM 
[Irec, ~] = tv2(TF, Imn, lam, betaadmm,  ...
    NIter, UB, 2, noisetype);
end 
save(fullfile(fold_out,'MRI_6_Spiral_SR_2_Sig_20_tv2_best_ssim.mat'),'Irec','lam'); 
% Compute SSIM 
temp1 = Im(1:240,1:240); 
temp2 = Irec(1:240,1:240); 
sc2 = sum(temp1(:).*temp2(:))./sum(temp2(:).^2); 
SSIM_Score = round(real(ssim(sc2*temp2,temp1,'Exponents',[1,1,1],'DynamicRange',255)),3); 
Irec_SSIM = Irec; 

clearvars -except Irec_SNR Irec_SSIM SNR_Score SSIM_Score 
toc; 

% Display Reconstruction Scores 
disp(sprintf('SNR=%.2f',SNR_Score)); 
disp(sprintf('SSIM=%.3f',SSIM_Score)); 
rmpath(genpath('Compared_Methods')) 
cd 'Scripts_View_And_Generate_Results' 
cd 'Scripts_Table1_Entries' 
cd 'MRI_6' 
