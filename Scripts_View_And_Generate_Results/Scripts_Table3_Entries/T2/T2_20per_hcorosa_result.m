
cd .. 
cd .. 
cd .. 
close all; 
clear all; 

% Setting Folder Paths 
addpath(genpath('Proposed_Method')); 
fold_in = fullfile('Models_and_Measurements','Measurements_for_Table3'); 
fold_out = fullfile('Generated_Reconstruction_Results','Reconstructions_in_Table3'); 

% Load Reference Image 
load(fullfile(fold_in,'T2_20per.mat'),'Im'); 
% Load Measurement 
load(fullfile(fold_in,'T2_20per.mat'),'Imn','TF');
disp('T2_20per')

tic;

% HCOROSA Reconstruction 
disp('HCOROSA Reconstruction')
% Load Lambda 
load(fullfile(fold_in,'Lambda_Table3.mat'),'Lambda_Table3');
lh_snr = Lambda_Table3(5,4); 
lh_ssim = Lambda_Table3(5,3); 
NIter = 200; % Number of ADMM iterations 
betaadmm = 1; % Parameter c in ADMM 
beta = 1; % Paramter tau in spatial weight penalty 
UB = 300; % Upperbound for image values, Effective bounding constraint is [0,UB] 
Ni = 30; % Number of CG iterations for evaluating restored image 
L = 4; % Number of Multi-Resolution levels 
Nf = 10; % Number of BCD iterations 
noisetype = 1; % 1 for AWGN 
lam = lh_snr; % Lambda as tuning paramter for regularization, Tuned for SNR 
[~, Irec, ~] = hcorosa(TF, Imn, lam, beta, betaadmm, ...
    NIter, Ni, Nf, L, UB, noisetype);
save(fullfile(fold_out,'T2_20per_hcorosa_best_snr.mat'),'Irec','lam'); 
% Compute SNR 
temp1 = Im(1:256,1:256); 
temp2 = Irec(1:256,1:256); 
sc1 = sum(temp1(:).*temp2(:))./sum(temp2(:).^2); 
SNR_Score = round(10*log10(sum(temp1(:).^2)/sum((sc1*temp2(:) - temp1(:)).^2)),2); 
Irec_SNR = Irec; 
if lh_snr ~= lh_ssim 
lam = lh_ssim; % Lambda as tuning paramter for regularization, Tuned for SSIM 
[~, Irec, ~] = hcorosa(TF, Imn, lam, beta, betaadmm, ...
    NIter, Ni, Nf, L, UB, noisetype);
end 
save(fullfile(fold_out,'T2_20per_hcorosa_best_ssim.mat'),'Irec','lam'); 
% Compute SSIM 
temp1 = Im(1:256,1:256); 
temp2 = Irec(1:256,1:256); 
sc2 = sum(temp1(:).*temp2(:))./sum(temp2(:).^2); 
SSIM_Score = round(real(ssim(sc2*temp2,temp1,'Exponents',[1,1,1],'DynamicRange',255)),3); 
Irec_SSIM = Irec; 

clearvars -except Irec_SNR Irec_SSIM SNR_Score SSIM_Score 
toc; 

% Display Reconstruction Scores 
disp(sprintf('SNR=%.2f',SNR_Score)); 
disp(sprintf('SSIM=%.3f',SSIM_Score)); 
rmpath(genpath('Proposed_Method')) 
cd 'Scripts_View_And_Generate_Results' 
cd 'Scripts_Table3_Entries' 
cd 'T2' 
