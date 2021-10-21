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
fo = pwd;
fo1 = fullfile(fo,'Models_and_Measurements','Measurements_for_Table4');
fo2 = fullfile(fo,'Generated_Reconstruction_Results','Reconstructions_in_Table4');
fo3 = fullfile(fo,'Scripts_View_And_Generate_Results','Scripts_Table4_Entries');
cd Scripts_View_And_Generate_Results
cd Scripts_View_Tables
fo4 = pwd;

Ni = 6;
Nt = 4;
Table_4 = zeros(Ni*(Nt+1),4);

Chk_In = input(['Enter y to generate the results instead of using stored reconstructions \n',...
    '(Note: Generating all results need GPU and can take hours depending on the computing hardware) \n Entry (y/n):'],'s');

if Chk_In=='y'||Chk_In=='Y'
    Chk_Run = 1;
else
    Chk_Run = 0;
end

for ind1 = 1:6
    
    for ind2 = 1:Nt
        
        switch ind2
            case 1
                fc2 = '_Spiral_SR_1_var_15';
            case 2
                fc2 = '_Spiral_SR_2_var_15';
            case 3
                fc2 = '_TSP_SR_1_var_15';
            otherwise
                fc2 = '_Radial_SR_1_var_15';
        end
        
        if ind1 < 4
            fn1 = ['MRI_',num2str(ind1)];
            fn3 = [fn1,fc2];
        else
            fn1 = ['MRI_',num2str(ind1+4)];
            fn3 = [fn1,fc2];
        end
        
        load(fullfile(fo1,[fn3,'.mat']),'Im');
        
        switch ind1
            
%             case 7
%                 cut_ind = 480;
%             case {6,8,9,10}
            case {4,5,6}
                cut_ind = 240;
            otherwise
                cut_ind = size(Im,1);
                
        end
        
        Im = Im(1:cut_ind,1:cut_ind);
        
        for ind3 = 1:2
            
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
                    dr = 255;
                otherwise
                    load(fullfile(fo2,[fn3,'_pnp_recon.mat']),'Irec');
                    I1 = Irec(1:cut_ind,1:cut_ind);
                    I2 = Irec(1:cut_ind,1:cut_ind);
                    dr = 255;
            end
            
            indr1 = (ind1-1)*(Nt+1) + ind2;
            s1 = sum(I1(:).*Im(:))./sum(I1(:).^2);
            s2 = sum(I2(:).*Im(:))./sum(I2(:).^2);
            
            Table_4(indr1,ind3) = round(real(ssim(s1*I1,Im,'Exponents',[1,1,1],'DynamicRange',dr)),3);
            Table_4(indr1,ind3+2) = round(10*log10(sum(Im(:).^2)/sum((s2*I2(:) - Im(:)).^2)),2);
            
        end
        
    end
    
end

clearvars -except Table_4
