% Hessian-Schatten Norm Regularization, HS
% Sanjay Viswanath, Muthuvel Arigovindan, Imaging Systems Lab, EE, IISc
%
% Function Output
% X: Restored Image by HS after full set of ADMM iterations
% cost: HS total cost from each iteration to check convergence
%
% Function Values
% TF: Transfer Function in Fourier domain, Sampling Trajectory for MRI
% Imn: Noisy image in restoration or complex Fourier measurement samples in MRI
% lamf: Lambda as tuning paramter for regularization term
% beta: Parameter beta in ADMM 
% NIter: Number of HS iterations 
% UB: Upperbound for image values, Effective bounding constraint is [0,UB]
% soregtype: Second order term, 1 for Hessian-Schattern, 2 for TV2
% noisetype: 1 for AWGN

function [X, cost] = hs(TF, Imn, lamf, beta,  ...
    NIter,   UB, soregtype, noisetype)

if gpuDeviceCount > 0 
    gpu_comp = 1;
else
    gpu_comp = 0;
end

if gpu_comp == 0
    error(' NVIDIA GPU required to run the code ');
end

TF = gpuArray(TF);
Imn = gpuArray(Imn);
lamf = gpuArray(lamf);
beta = gpuArray(beta);
NIter = gpuArray(NIter);
UB = gpuArray(UB);
soregtype = gpuArray(soregtype);
noisetype = gpuArray(noisetype);

Sz = gpuArray(size(Imn));
tmp = zeros(Sz,'gpuArray'); tmp(1) = 1;
Q = real(ifft2(conj(TF).*TF)) +  tmp +  GetQQr2(Sz);
Q = fft2(Q); clear tmp;

X = zeros((Sz),'gpuArray');

lamcb  = zeros(size(Imn),'gpuArray');
lamcm  = zeros(size(Imn),'gpuArray');
lamcr2 = zeros([size(Imn) 3],'gpuArray');

Rb = X;
Rm = (fft2(TF.*fft2(X)));
R2 = GetD2(X);

for ind = 1:NIter
    
    Zb  = proxb(Rb, lamcb,  beta, UB);
    Zm  = proxm(Rm, lamcm,  Imn,  beta, noisetype);
    Zr2 = prox2(R2, lamcr2, lamf, beta,  soregtype);
    
    
    B =   Zb  - lamcb/beta + ...
        real(ifft2(conj(TF).*fft2(Zm  - lamcm/beta))) + ...
        Div2(Zr2-lamcr2/beta);
    
    
    X = real(ifft2(fft2(B)./Q));
    
    
    Rb = X;
    Rm = (ifft2(TF.*fft2(X)));
    R2 = GetD2(X);
    
    lamcb  = lamcb  + beta*(Rb - Zb);
    lamcm  = lamcm  + beta*(Rm - Zm);
    lamcr2 = lamcr2 + beta*(R2 - Zr2);
    
    cost(ind) = CostNLDecon ( X, TF, Imn,  lamf,  ...
        soregtype );
    
    
    
end

X = gather(X);
cost = gather(cost);

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function cost = CostNLDecon ( X, TF, Imn, lamf,  soregtype  )



cost = abs( (ifft2(TF.*fft2(X))) - Imn).^2;


if soregtype == 2
    R2 = sqrt(ApplyDerSqiso2(X));
else
    R2 = GetHS1( X );
end

cost = sum(cost(:))  + lamf*sum(R2(:));


return;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function R = ApplyDerSqiso2(X)



R  = (-2*X + X(:,[end 1:end-1]) + X(:,[2:end 1])).^2;

R  = R + (-2*X + X([end 1:end-1],:) + X([2:end 1],:)).^2;

R  = R +  2*(   X + X([end 1:end-1],[end 1:end-1])  ...
    - X([end 1:end-1],:) - X(:,[end 1:end-1]) ).^2;



return;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function hs1 = GetHS1( X )


a = (-2*X + X(:,[end 1:end-1]) + X(:,[2:end 1]));

d = (-2*X + X([end 1:end-1],:) + X([2:end 1],:));

b =  (   X + X([end 1:end-1],[end 1:end-1])  ...
    - X([end 1:end-1],:) - X(:,[end 1:end-1]) );


lam1 = 0.5*( (a + d) + sqrt((a-d).^2 + 4*b.^2) );
lam2 = 0.5*( (a + d) - sqrt((a-d).^2 + 4*b.^2) );

hs1 = abs(lam1) + abs(lam2);


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function R = GetQQr2(Sz)

X = zeros(Sz);
X(1) = 1;

tmp = (-2*X + X(:,[end 1:end-1]) + X(:,[2:end 1]));
tmp =  -2*tmp + tmp(:,[end 1:end-1]) + tmp(:,[2:end 1]);
R =  tmp;

tmp = (-2*X + X([end 1:end-1],:) + X([2:end 1],:));
tmp =  -2*tmp + tmp([end 1:end-1],:) + tmp([2:end 1],:);
R = R + tmp;

tmp =  (   X + X([end 1:end-1],[end 1:end-1])  ...
    - X([end 1:end-1],:) - X(:,[end 1:end-1]) );
tmp =    tmp + tmp([2:end 1],[2:end 1])  ...
    - tmp([2:end 1],:) - tmp(:,[2:end 1]) ;

R = (R + 2*tmp);


return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D2 = GetD2 ( X )


D2(:,:,1) = -2*X + X(:,[end 1:end-1]) + X(:,[2:end 1]);


D2(:,:,2) = -2*X + X([end 1:end-1],:) + X([2:end 1],:);


D2(:,:,3) = sqrt(2)*(X + X([end 1:end-1],[end 1:end-1])  ...
    - X([end 1:end-1],:) - X(:,[end 1:end-1]));


return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Zm  = proxm(Rm, lamcm,  Imn,  beta, noisetype)

if noisetype == 1
    
    Rm = Rm + lamcm/beta;
    Zm = (Imn + beta*Rm)/(1+beta);
    
end


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Zb  = proxb(Rb, lamcb,  beta, UB)

Rb = Rb + lamcb/beta;

IND1 = Rb < 0;
IND2 = Rb > UB;

Zb = (~IND1).*Rb.*(~IND2) + UB*IND2;

return;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Zr2 = prox2(D2, lamcr2, lamf, beta,  soregtype)

D2 = D2 + lamcr2/beta;

if soregtype == 2
    
    r = sqrt(D2(:,:,1).^2 + D2(:,:,2).^2 + D2(:,:,3).^2);
    rn = r - lamf/beta;
    rn = rn.*(rn >= 0);
    r  = r + (r < 1e-12);
    Zr2 = bsxfun(@times, D2, rn./r);
    
else
    
    Zr2 = proxhs1(D2, lamf/beta);
    
end


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Zr2 = proxhs1( D2, rfac )


a = D2(:,:,1);
d = D2(:,:,2);
b = D2(:,:,3)/sqrt(2);


lam1 = 0.5*( (a + d) + sqrt((a-d).^2 + 4*b.^2) );
lam2 = 0.5*( (a + d) - sqrt((a-d).^2 + 4*b.^2) );


tmp1 = (a - lam1) + b; tmp2 = b + (d-lam1);
E1(:,:,1) = -tmp2;  E1(:,:,2) = tmp1;
n1 = sqrt(E1(:,:,1).^2 + E1(:,:,2).^2);
n1 = n1 + (n1 == 0);
E1 = bsxfun(@times, E1, 1./n1);

tmp1 = (a - lam2) + b; tmp2 = b + (d-lam2);
E2(:,:,1) = -tmp2;  E2(:,:,2) = tmp1;
n2 = sqrt(E2(:,:,1).^2 + E2(:,:,2).^2);
n2 = n2 + (n2 == 0);
E2 = bsxfun(@times, E2, 1./n2);


SGN1 = sign(lam1);
SGN2 = sign(lam2);

lam1 = abs(lam1) - rfac;
lam1 = lam1.*(lam1 > 0);
lam1 = lam1.*SGN1;

lam2 = abs(lam2) - rfac;
lam2 = lam2.*(lam2 > 0);
lam2 = lam2.*SGN2;


a = E1(:,:,1).*E1(:,:,1).*lam1 + E2(:,:,1).*E2(:,:,1).*lam2;
d = E1(:,:,2).*E1(:,:,2).*lam1 + E2(:,:,2).*E2(:,:,2).*lam2;
b = E1(:,:,1).*E1(:,:,2).*lam1 + E2(:,:,1).*E2(:,:,2).*lam2;

Zr2(:,:,1) = a;
Zr2(:,:,2) = d;
Zr2(:,:,3) = sqrt(2)*b;


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function X = Div2(D2)


X =  -2*D2(:,:,1) + D2(:,[end 1:end-1],1) + D2(:,[2:end 1],1);

X  =  X + (-2*D2(:,:,2) + D2([end 1:end-1],:,2) + D2([2:end 1],:,2));

X =   X + sqrt(2)*( D2(:,:,3) + D2([2:end 1],[2:end 1],3)  ...
    - D2([2:end 1],:,3) - D2(:,[2:end 1],3));


return;



