
cd .. 
cd .. 
cd .. 
close all; 
clear all; 

% Setting Folder Paths 
addpath(genpath('Compared_Methods')); 
fold_in = fullfile('Models_and_Measurements','Measurements_for_Table2'); 
fold_ref = fullfile('Models_and_Measurements','Reference_Models'); 
fold_out = fullfile('Generated_Reconstruction_Results','Reconstructions_in_Table2'); 

% Load Reference Image 
load(fullfile(fold_ref,'MRI_3.mat'),'Im'); 
% Load Measurement 
load(fullfile(fold_in,'MRI_3_Radial_SR_1_Sig_20.mat'),'Imn','TF');
disp('MRI_3_Radial_SR_1_Sig_20')

tic;

% TGV2 Reconstruction 
disp('TGV2 Reconstruction')
% Load Lambda 
load(fullfile(fold_in,'Lambda_Table2.mat'),'Lambda_Table2');
ltgv2_snr = Lambda_Table2(13,2); 
ltgv2_ssim = Lambda_Table2(13,1); 
eps1 = 1e-5; %  Smooth approximation paramter for TV functional 
No = 50; % Number of TGV2 iterations 
Ni = 100; % Number of inner CG iterations 
alpha1 = ones(size(Imn)); % Spatially varying weight for first order term 
alpha2 = ones(size(Imn)); % Spatially varying weight for second order term 
lam = ltgv2_snr; % Lambda2 as tuning paramter for second order term, Tuned for SNR 
lam1 = 10*lam; % Lambda as tuning paramter for first order term, Tuned for SNR 
[X,~,~] = tgv2(Imn, TF, alpha1, alpha2, lam1, lam, eps1, Ni, No);
Irec = X(:,:,1); 
save(fullfile(fold_out,'MRI_3_Radial_SR_1_Sig_20_tgv2_best_snr.mat'),'Irec','lam'); 
% Compute SNR 
temp1 = Im(1:256,1:256); 
temp2 = Irec(1:256,1:256); 
sc1 = sum(temp1(:).*temp2(:))./sum(temp2(:).^2); 
SNR_Score = round(10*log10(sum(temp1(:).^2)/sum((sc1*temp2(:) - temp1(:)).^2)),2); 
Irec_SNR = Irec; 
if ltgv2_snr ~= ltgv2_ssim 
lam = ltgv2_ssim; % Lambda as tuning paramter for second order term, Tuned for SSIM 
lam1 = 10*lam; % Lambda as tuning paramter for first order term, Tuned for SNR 
[X,~,~] = tgv2(Imn, TF, alpha1, alpha2, lam1, lam, eps1, Ni, No);
Irec = X(:,:,1); 
end 
save(fullfile(fold_out,'MRI_3_Radial_SR_1_Sig_20_tgv2_best_ssim.mat'),'Irec','lam'); 
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
rmpath(genpath('Compared_Methods')) 
cd 'Scripts_View_And_Generate_Results' 
cd 'Scripts_Table2_Entries' 
cd 'MRI_3' 
