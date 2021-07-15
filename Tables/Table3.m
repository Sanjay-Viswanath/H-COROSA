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
fo = pwd;
fo1 = fullfile(pwd,'MRI_Data','Experiment_3');
fo2 = fullfile(pwd,'MAT_Files','Experiment_3');
img_ind = [25,28,35,39,45,50];
cd Tables

Ni = 6;
Nt = 4;
Table_3 = zeros(Ni*(Nt+1),4);

for ind1 = 1:Ni
    
    for ind2 = 1:2
        
        if ind2 == 1
            fn3 = ['T',num2str(ind1),'_10per'];
            fn4 = 'DAGAN_10per_noiseless.mat';
        else
            fn3 = ['T',num2str(ind1),'_20per'];
            fn4 = 'DAGAN_20per_noiseless.mat';
        end
        
                
        for ind3 = 1:Nt
            
            switch ind3
                
                case 1
                    load(fullfile(fo1,[fn3,'.mat']),'Im');
                    load(fullfile(fo2,[fn3,'_hcorosa_best_ssim.mat']),'Irec');
                    I1 = Irec;
                    load(fullfile(fo2,[fn3,'_hcorosa_best_snr.mat']),'Irec');
                    I2 = Irec;
                    dr = 255;
                case 2
                    load(fullfile(fo1,[fn3,'.mat']),'Im');
                    load(fullfile(fo2,[fn3,'_corosa_best_ssim.mat']),'Irec');
                    I1 = Irec;
                    load(fullfile(fo2,[fn3,'_corosa_best_snr.mat']),'Irec');
                    I2 = Irec;
                    dr = 255;
                case 3
                    load(fullfile(fo2,fn4),'x_im','x_rec');
                    Im = reshape(x_im(img_ind(ind1),:,:),[256,256]);
                    I1 = double(reshape(x_rec(img_ind(ind1),:,:),[256,256]));
                    I1 = rescale(I1,min(Im(:)),max(Im(:)));
                    I2 = I1;
                    dr = 2;
                otherwise
                    load(fullfile(fo1,[fn3,'.mat']),'Im');
                    load(fullfile(fo2,[fn3,'_pnp_recon.mat']),'Irec');
                    I1 = Irec;
                    I2 = Irec;
                    dr = 255;
            end
            
            indr1 = (ind1-1)*(Nt+1) + ind3;
            s1 = sum(I1(:).*Im(:))./sum(I1(:).^2);
            s2 = sum(I2(:).*Im(:))./sum(I2(:).^2);
            
            Table_3(indr1,ind2) = round(real(ssim(s1*I1,Im,'Exponents',[1,1,1],'DynamicRange',dr)),3);
            Table_3(indr1,ind2+2) = round(10*log10(sum(Im(:).^2)/sum((s2*I2(:) - Im(:)).^2)),2);
            
        end
        
    end
    
end
