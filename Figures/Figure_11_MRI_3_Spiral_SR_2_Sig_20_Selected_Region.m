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
fi = pwd;
cd Figures;
fi1 = fullfile(fi,'MRI_Data','Experiment_1');
fi2 = fullfile(fi,'MAT_Files','Experiment_1');

fn = 'MRI_3_Spiral_SR_2_Sig_20';
load(fullfile(fi1,[fn,'.mat']),'Im','Imn');
Imn = abs(Imn);

load(fullfile(fi2,[fn,'_hcorosa_best_ssim.mat']),'Irec')
s = sum(Im(:).*Irec(:))/sum(Irec(:).^2);
It1 = s*Irec;

load(fullfile(fi2,[fn,'_corosa_best_ssim.mat']),'Irec')
s = sum(Im(:).*Irec(:))/sum(Irec(:).^2);
It2 = s*Irec;


load(fullfile(fi2,[fn,'_hs_best_ssim.mat']),'Irec')
s = sum(Im(:).*Irec(:))/sum(Irec(:).^2);
It3 = s*Irec;

x1 = 120; x2 = 160; y1 = 90; y2 = 120;
xsp = 255*ones(y2-y1+1,2);
n = 5;
ysp = 255*ones(2,n*(x2-x1+1)+(n+1)*size(xsp,2));
% I2 = [ysp1;xsp,Im(y1:y2,x1:x2),xsp,Imn(y1:y2,x1:x2),xsp,It1(y1:y2,x1:x2),xsp;ysp1;...
%     xsp,It2(y1:y2,x1:x2),xsp,It3(y1:y2,x1:x2),xsp,It4(y1:y2,x1:x2),xsp;ysp2];
I2 = [ysp;xsp,Im(y1:y2,x1:x2),xsp,Imn(y1:y2,x1:x2),xsp,...
    It1(y1:y2,x1:x2),xsp,It2(y1:y2,x1:x2),xsp,It3(y1:y2,x1:x2),xsp,;ysp];
Fig1 = figure;
imshow(I2,[]);
xi = [140,150];
yi = [98,98];

line(2+xi-x1+1,2+yi-y1+1,'color','red','LineWidth', 1);

Imp = improfile(Im,xi,yi);
Imnp = improfile(Imn,xi,yi);
It1p = improfile(It1,xi,yi);
It2p = improfile(It2,xi,yi);
It3p = improfile(It3,xi,yi);

Fig2 = figure;
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

Fig3 = figure;
imshow(Im,[]);
rectangle('Position',[x1,y1,x2-x1+1,y2-y1+1],...
  'EdgeColor', 'r',...
  'LineWidth', 2,...
  'LineStyle','-');
