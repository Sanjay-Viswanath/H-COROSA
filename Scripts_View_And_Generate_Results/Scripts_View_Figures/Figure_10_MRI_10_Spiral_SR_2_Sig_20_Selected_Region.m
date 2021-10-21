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
fi1 = fullfile(fi,'Models_and_Measurements','Measurements_for_Table1');
fi2 = fullfile(fi,'Generated_Reconstruction_Results','Reconstructions_in_Table1');

Chk_In = input(['Enter y to generate the results instead of using stored reconstructions \n',...
    '(Note: Generating all results need GPU and can take an hour or more depending on the computing hardware) \n Entry (y/n):'],'s');

if Chk_In=='y'||Chk_In=='Y'
    Chk_Run = 1;
else
    Chk_Run = 0;
end

fn = 'MRI_10_Spiral_SR_2_Sig_20';
load(fullfile(fi1,[fn,'.mat']),'Im','Imn');
Imn = abs(Imn);

fr = fullfile(fi2,[fn,'_hcorosa_best_ssim.mat']);
if exist(fr,'file') && Chk_Run==0
    load(fr,'Irec')
else
    cd ..
    cd Scripts_Table1_Entries
    cd MRI_10
    save('temp_store.mat');
    MRI_10_Spiral_SR_2_Sig_20_hcorosa_result
    load('temp_store.mat');
    delete temp_store.mat
    cd ..
    cd ..
    cd Scripts_View_Figures
    load(fr,'Irec')
end
s = sum(Im(:).*Irec(:))/sum(Irec(:).^2);
It1 = s*Irec;

fr = fullfile(fi2,[fn,'_corosa_best_ssim.mat']);
if exist(fr,'file') && Chk_Run==0
    load(fr,'Irec')
else
    cd ..
    cd Scripts_Table1_Entries
    cd MRI_10
    save('temp_store.mat');
    MRI_10_Spiral_SR_2_Sig_20_corosa_result
    load('temp_store.mat');
    delete temp_store.mat
    cd ..
    cd ..
    cd Scripts_View_Figures
    load(fr,'Irec')
end
s = sum(Im(:).*Irec(:))/sum(Irec(:).^2);
It2 = s*Irec;


fr = fullfile(fi2,[fn,'_hs_best_ssim.mat']);
if exist(fr,'file') && Chk_Run==0
    load(fr,'Irec')
else
    cd ..
    cd Scripts_Table1_Entries
    cd MRI_10
    save('temp_store.mat');
    MRI_10_Spiral_SR_2_Sig_20_hs_result
    load('temp_store.mat');
    delete temp_store.mat
    cd ..
    cd ..
    cd Scripts_View_Figures
    load(fr,'Irec')
end
s = sum(Im(:).*Irec(:))/sum(Irec(:).^2);
It3 = s*Irec;

x1 = 94; x2 = 174; y1 = 160; y2 = 220;
xsp = 255*ones(y2-y1+1,2);
n = 5;
ysp = 255*ones(2,n*(x2-x1+1)+(n+1)*size(xsp,2));
% I2 = [ysp1;xsp,Im(y1:y2,x1:x2),xsp,Imn(y1:y2,x1:x2),xsp,It1(y1:y2,x1:x2),xsp;ysp1;...
%     xsp,It2(y1:y2,x1:x2),xsp,It3(y1:y2,x1:x2),xsp,It4(y1:y2,x1:x2),xsp;ysp2];
I2 = [ysp;xsp,Im(y1:y2,x1:x2),xsp,Imn(y1:y2,x1:x2),xsp,...
    It1(y1:y2,x1:x2),xsp,It2(y1:y2,x1:x2),xsp,It3(y1:y2,x1:x2),xsp,;ysp];
Fig1 = figure('Name',['Fig.10.Reconstruction of I10 from 20% spiral samples:'...
    '(b) Comparison of selected region from reconstructions along with the scanline location,',...
    ' Left to Right: Original, Measured, HCOROSA, COROSA, HS']...
    ,'NumberTitle','off');
imshow(I2,[]);
xi = [125,134];
yi = [185,191];
line(2+xi-x1+1,2+yi-y1+1,'color','red','LineWidth', 1);

Imp = improfile(Im,xi,yi);
Imnp = improfile(Imn,xi,yi);
It1p = improfile(It1,xi,yi);
It2p = improfile(It2,xi,yi);
It3p = improfile(It3,xi,yi);

Fig2 = figure('Name',['Fig.10.Reconstruction of I10 from 20% spiral samples:'...
    '(c) Intensity profile along a scanline in the selected region.']...
    ,'NumberTitle','off');
temp = 1:length(Imp);
p = plot(temp,Imp,'r','LineWidth',2);
hold on;
plot(temp,It1p,'b','LineWidth',2);
plot(temp,It2p,'g','LineWidth',2);
xlabel('Pixels along line from left to right','FontSize',18,'FontWeight','Bold','FontName','Times');
ylabel('Pixel Intensity','FontSize',18,'FontWeight','Bold','FontName','Times');
legend('Original Image','H-COROSA','COROSA','FontSize',18,'FontWeight','Bold','FontName','Times');
h = gca;
h.FontSize = 18;
h.FontName = 'Times';
h.Box = 'off';

Fig3 = figure('Name','Fig.10.Reconstruction of I10 from 20% spiral samples: (a) Original I2 with selected region',...
    'NumberTitle','off');
imshow(Im(1:240,1:240),[]);
rectangle('Position',[x1,y1,x2-x1+1,y2-y1+1],...
  'EdgeColor', 'r',...
  'LineWidth', 2,...
  'LineStyle','-');
