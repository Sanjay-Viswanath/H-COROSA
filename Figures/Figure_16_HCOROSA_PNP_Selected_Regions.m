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
fi1 = fullfile(fi,'MRI_Data','Experiment_4');
fi2 = fullfile(fi,'MAT_Files','Experiment_4');

fn = 'MRI_8_TSP_SR_1_var_15';
load(fullfile(fi1,[fn,'.mat']),'Im','Imn');
load(fullfile(fi2,[fn,'_hcorosa_best_ssim.mat']),'Irec');
s = sum(Irec(:).*Im(:))/sum(Irec(:).^2);
It2 = s*Irec;
load(fullfile(fi2,[fn,'_pnp_recon.mat']),'Irec');
s = sum(Irec(:).*Im(:))/sum(Irec(:).^2);
It3 = s*Irec;

x1 = 30; x2 = 90; y1 = 140; y2 = 175;
xsp = 255*ones(y2-y1+1,2);
n = 4;
ysp = 255*ones(2,n*(x2-x1+1)+(n+1)*size(xsp,2));
I1 = [ysp;xsp,Im(y1:y2,x1:x2),xsp,abs(Imn(y1:y2,x1:x2)),xsp,It2(y1:y2,x1:x2),xsp,It3(y1:y2,x1:x2),xsp;ysp];
Fig1 = figure;
imshow(I1,[],'border','tight');
xi = [47,46];
yi = [159,147];
line(xi-x1+1,yi-y1+1,'color','red','LineWidth', 1);

Imp = improfile(Im,xi,yi);
It2p = improfile(It2,xi,yi);
It3p = improfile(It3,xi,yi);
Fig4 = figure;
temp = 1:length(Imp);
p = plot(temp,Imp,'r','LineWidth',2);
hold on;
plot(temp,It2p,'b','LineWidth',2);
plot(temp,It3p,'g','LineWidth',2);
xlabel('Pixels along line from left to right','FontSize',18,'FontWeight','Bold','FontName','Times');
ylabel('Pixel Intensity','FontSize',18,'FontWeight','Bold','FontName','Times');
legend('Original Image','H-COROSA','PNP','FontSize',18,'FontWeight','Bold','FontName','Times');
h = gca;
h.FontSize = 18;
h.FontName = 'Times';
h.Box = 'off';

fn = 'MRI_2_Spiral_SR_1_var_15';
load(fullfile(fi1,[fn,'.mat']),'Im','Imn');
load(fullfile(fi2,[fn,'_hcorosa_best_ssim.mat']),'Irec');
s = sum(Irec(:).*Im(:))/sum(Irec(:).^2);
It2 = s*Irec;
load(fullfile(fi2,[fn,'_pnp_recon.mat']),'Irec');
s = sum(Irec(:).*Im(:))/sum(Irec(:).^2);
It3 = s*Irec;

x1 = 130; x2 = 152; y1 = 145; y2 = 180;
xsp = 255*ones(y2-y1+1,2);
n = 1;
ysp = 255*ones(2,n*(x2-x1+1)+(n+1)*size(xsp,2));
I2 = [ysp;xsp,Im(y1:y2,x1:x2),xsp;ysp;xsp,abs(Imn(y1:y2,x1:x2)),xsp;ysp;...
    xsp,It2(y1:y2,x1:x2),xsp;ysp;xsp,It3(y1:y2,x1:x2),xsp;ysp];
Fig2 = figure;
imshow(I2,[],'border','tight');
xi = [136,145];
yi = [160,160];
line(xi-x1+1,yi-y1+1,'color','red','LineWidth', 1);

Imp = improfile(Im,xi,yi);
It2p = improfile(It2,xi,yi);
It3p = improfile(It3,xi,yi);
Fig5 = figure;
temp = 1:length(Imp);
p = plot(temp,Imp,'r','LineWidth',2);
hold on;
plot(temp,It2p,'b','LineWidth',2);
plot(temp,It3p,'g','LineWidth',2);
xlabel('Pixels along line from left to right','FontSize',18,'FontWeight','Bold','FontName','TimesNewRoman');
ylabel('Pixel Intensity','FontSize',18,'FontWeight','Bold','FontName','Times');
legend('Original Image','H-COROSA','PNP','FontSize',18,'FontWeight','Bold','FontName','TimesNewRoman');
h = gca;
h.FontSize = 18;
h.FontName = 'Times';
h.Box = 'off';

fn = 'MRI_3_Spiral_SR_2_var_15';
load(fullfile(fi1,[fn,'.mat']),'Im','Imn');
load(fullfile(fi2,[fn,'_hcorosa_best_ssim.mat']),'Irec');
s = sum(Irec(:).*Im(:))/sum(Irec(:).^2);
It2 = s*Irec;
load(fullfile(fi2,[fn,'_pnp_recon.mat']),'Irec');
s = sum(Irec(:).*Im(:))/sum(Irec(:).^2);
It3 = s*Irec;

x1 = 120; x2 = 160; y1 = 90; y2 = 120;
xsp = 255*ones(y2-y1+1,2);
n = 2;
ysp1 = 255*ones(15,n*(x2-x1+1)+(n+1)*size(xsp,2));
ysp = 255*ones(2,n*(x2-x1+1)+(n+1)*size(xsp,2));
I3 = [ysp;xsp,Im(y1:y2,x1:x2),xsp,abs(Imn(y1:y2,x1:x2)),xsp;ysp;xsp,It2(y1:y2,x1:x2),xsp,It3(y1:y2,x1:x2),xsp;ysp];
Fig3 = figure;
imshow(I3,[],'border','tight');
xi = [140,150];
yi = [98,98];
line(xi-x1+1,yi-y1+1,'color','red','LineWidth', 1);

Imp = improfile(Im,xi,yi);
It2p = improfile(It2,xi,yi);
It3p = improfile(It3,xi,yi);
Fig6 = figure;
temp = 1:length(Imp);
p = plot(temp,Imp,'r','LineWidth',2);
hold on;
plot(temp,It2p,'b','LineWidth',2);
plot(temp,It3p,'g','LineWidth',2);
xlabel('Pixels along line from left to right','FontSize',18,'FontWeight','Bold','FontName','Times');
ylabel('Pixel Intensity','FontSize',18,'FontWeight','Bold','FontName','Times');
legend('Original Image','H-COROSA','PNP','FontSize',18,'FontWeight','Bold','FontName','Times');
h = gca;
h.FontSize = 18;
h.FontName = 'Times';
h.Box = 'off';
