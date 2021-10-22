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
fo = pwd;
fo1 = fullfile(fo,'Models_and_Measurements','Measurements_for_Table1');
fo2 = fullfile(fo,'Generated_Reconstruction_Results','Reconstructions_in_Table1');
fo3 = fullfile(fo,'Scripts_View_And_Generate_Results','Scripts_Table1_Entries');
cd Scripts_View_And_Generate_Results
cd Scripts_View_Tables
fo4 = pwd;

Ni = 10;
Nt = 5;
Table_1 = zeros(Ni*(Nt+1),4);

Chk_In = input(['Enter y to generate the results instead of using stored reconstructions \n',...
    '(Note: Generating all results need GPU and can take hours depending on the computing hardware) \n Entry (y/n):'],'s');

if Chk_In=='y'||Chk_In=='Y'
    Chk_Run = 1;
else
    Chk_Run = 0;
end

for ind1 = 1:Ni
    
    for ind2 = 1:2
        
        fn1 = ['MRI_',num2str(ind1)];
        fn3 = [fn1,'_Spiral_SR_',num2str(ind2),'_Sig_20'];
        load(fullfile(fo1,[fn3,'.mat']),'Im');
        
        switch ind1
            
            case 7
                cut_ind = 480;
            case {6,8,9,10}
                cut_ind = 240;
            otherwise
                cut_ind = size(Im,1);
                
        end
        Im = Im(1:cut_ind,1:cut_ind);
        
        for ind3 = 1:Nt
            
            switch ind3
                
                case 1
                    fr = fullfile(fo2,[fn3,'_hcorosa_best_ssim.mat']);
                    if exist(fr,'file') && Chk_Run==0
                        load(fr,'Irec');
                    else
                        cd(fullfile(fo3,fn1));
                        save('temp_store.mat');
                        eval([fn3,'_hcorosa_result']);
                        load('temp_store.mat');
                        delete temp_store.mat
                        cd(fo4);
                        load(fr,'Irec');
                    end                        
                    I1 = Irec(1:cut_ind,1:cut_ind);
                    load(fullfile(fo2,[fn3,'_hcorosa_best_snr.mat']),'Irec');
                    I2 = Irec(1:cut_ind,1:cut_ind);
                case 2
                    fr = fullfile(fo2,[fn3,'_corosa_best_ssim.mat']);
                    if exist(fr,'file') && Chk_Run==0
                        load(fr,'Irec');
                    else
                        cd(fullfile(fo3,fn1));
                        save('temp_store.mat');
                        eval([fn3,'_corosa_result']);
                        load('temp_store.mat');
                        delete temp_store.mat
                        cd(fo4);
                        load(fr,'Irec');
                    end
                    I1 = Irec(1:cut_ind,1:cut_ind);
                    load(fullfile(fo2,[fn3,'_corosa_best_snr.mat']),'Irec');
                    I2 = Irec(1:cut_ind,1:cut_ind);
                case 3
                    fr = fullfile(fo2,[fn3,'_tgv2_best_ssim.mat']);
                    if exist(fr,'file') && Chk_Run==0
                        load(fr,'Irec');
                    else
                        cd(fullfile(fo3,fn1));
                        save('temp_store.mat');
                        eval([fn3,'_tgv2_result']);
                        load('temp_store.mat');
                        delete temp_store.mat
                        cd(fo4);
                        load(fr,'Irec');
                    end
                    I1 = Irec(1:cut_ind,1:cut_ind);
                    load(fullfile(fo2,[fn3,'_tgv2_best_snr.mat']),'Irec');
                    I2 = Irec(1:cut_ind,1:cut_ind);
                case 4
                    fr = fullfile(fo2,[fn3,'_tv2_best_ssim.mat']);
                    if exist(fr,'file') && Chk_Run==0
                        load(fr,'Irec');
                    else
                        cd(fullfile(fo3,fn1));
                        save('temp_store.mat');
                        eval([fn3,'_tv2_result']);
                        load('temp_store.mat');
                        delete temp_store.mat
                        cd(fo4);
                        load(fr,'Irec');
                    end
                    I1 = Irec(1:cut_ind,1:cut_ind);
                    load(fullfile(fo2,[fn3,'_tv2_best_snr.mat']),'Irec');
                    I2 = Irec(1:cut_ind,1:cut_ind);
                case 5
                    fr = fullfile(fo2,[fn3,'_hs_best_ssim.mat']);
                    if exist(fr,'file') && Chk_Run==0
                        load(fr,'Irec');
                    else
                        cd(fullfile(fo3,fn1));
                        save('temp_store.mat');
                        eval([fn3,'_hs_result']);
                        load('temp_store.mat');
                        delete temp_store.mat
                        cd(fo4);
                        load(fr,'Irec');
                    end
                    I1 = Irec(1:cut_ind,1:cut_ind);
                    load(fullfile(fo2,[fn3,'_hs_best_snr.mat']),'Irec');
                    I2 = Irec(1:cut_ind,1:cut_ind);
            end
            
            indr1 = (ind1-1)*(Nt+1) + ind3;
            s1 = sum(I1(:).*Im(:))./sum(I1(:).^2);
            s2 = sum(I2(:).*Im(:))./sum(I2(:).^2);
            
            Table_1(indr1,ind2) = round(real(ssim(s1*I1,Im,'Exponents',[1,1,1],'DynamicRange',255)),3);
            Table_1(indr1,ind2+2) = round(10*log10(sum(Im(:).^2)/sum((s2*I2(:) - Im(:)).^2)),2);
            
        end
        
    end
    
end

clearvars -except Table_1
