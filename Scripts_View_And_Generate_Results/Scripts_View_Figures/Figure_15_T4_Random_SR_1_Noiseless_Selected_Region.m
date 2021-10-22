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
fi1 = fullfile(fi,'Models_and_Measurements','Measurements_for_Table3');
fi2 = fullfile(fi,'Generated_Reconstruction_Results','Reconstructions_in_Table3');

Chk_In = input(['Enter y to generate the results instead of using stored reconstructions \n',...
    '(Note: Generating all results need GPU and can take an hour or more depending on the computing hardware) \n Entry (y/n):'],'s');

if Chk_In=='y'||Chk_In=='Y'
    Chk_Run = 1;
else
    Chk_Run = 0;
end

fn = fullfile(fi1,'T4_10per.mat');
load(fn,'Im','TF','Imn');
It1 = Im;
It2 = 255*fftshift(TF);
It3 = abs(Imn);

fr = fullfile(fi2,'T4_10per_hcorosa_best_ssim.mat');
if exist(fr,'file') && Chk_Run==0
    load(fr,'Irec')
else
    cd ..
    cd Scripts_Table3_Entries
    cd T4
    save('temp_store.mat');
    T4_10per_hcorosa_result
    load('temp_store.mat');
    delete temp_store.mat
    cd ..
    cd ..
    cd Scripts_View_Figures
    load(fr,'Irec')
end
s = sum(Irec(:).*Im(:))/sum(Irec(:).^2);
It4 = s*Irec;

load(fullfile(fi2,'DAGAN_10per_noiseless.mat'),'x_rec');
img_ind = 39;
Irec = double(reshape(x_rec(img_ind,:,:),[256,256]));
Irec = rescale(Irec,min(Im(:)),max(Im(:)));
s = sum(Irec(:).*Im(:))/sum(Irec(:).^2);
It5 = s*Irec;

fn = fullfile(fi2,'T4_10per_pnp_recon.mat');
load(fn,'Irec');
s = sum(Irec(:).*Im(:))/sum(Irec(:).^2);
It6 = s*Irec;

% Itt = 255*ones(size(It12));

x1 = 120; x2 = 170; y1 = 70; y2 = 150;
xsp = 255*ones(y2-y1+1,2);
n = 5;
ysp = 255*ones(2,n*(x2-x1+1)+(n+1)*size(xsp,2));
I = [ysp;xsp,It1(y1:y2,x1:x2),xsp,It3(y1:y2,x1:x2),xsp...
    ,It4(y1:y2,x1:x2),xsp...
    ,It5(y1:y2,x1:x2),xsp,It6(y1:y2,x1:x2),xsp;ysp];
Fig1 = figure('Name',['Fig.15. Reconstruction of T4 from 10% random samples:'...
    '(b) Comparison of selected region from reconstructions along with the scanline location,',...
    ' Left to Right: Original, Measured, HCOROSA, DAGAN, PNP']...
    ,'NumberTitle','off');
imshow(I,[]);
xi = [161,167];
yi = [114,125];
line(2+xi-x1+1,2+yi-y1+1,'color','red');

Imp = improfile(Im,xi,yi);
It1p = improfile(It4,xi,yi);
It2p = improfile(It5,xi,yi);
It3p = improfile(It6,xi,yi);
Fig2 = figure('Name',['Fig.15. Reconstruction of T4 from 10% random samples:'...
    '(c) Intensity profile along a scanline in the selected region.']...
    ,'NumberTitle','off');
p = plot(Imp(4:end),'r','LineWidth',2);
hold on;
plot(It1p(4:end),'b','LineWidth',2);
plot(It2p(4:end),'k','LineWidth',2);
plot(It3p(4:end),'g','LineWidth',2);
xlabel('Pixels along line from left to right','FontSize',18,'FontWeight','Normal','FontName','TimesNewRoman');
ylabel('Pixel Intensity','FontSize',18,'FontWeight','Normal','FontName','TimesNewRoman');
legend('Original Image','H-COROSA','DAGAN','PNP','FontSize',18,'FontWeight','Normal','FontName','TimesNewRoman');

Fig3 = figure('Name','Fig.15. Reconstruction of T4 from 10% random samples: (a) Original T4 with selected region',...
    'NumberTitle','off');
imshow(Im(1:256,1:256),[]);
x1 = 121; x2 = 170; y1 = 70; y2 = 150;
rectangle('Position',[x1,y1,x2-x1+1,y2-y1+1],...
  'EdgeColor', 'r',...
  'LineWidth', 1,...
  'LineStyle','-');
