%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COROSA regularization, Implementation of "Image restoration by combined 
% order regularization with optimal spatial adaptation", IEEE TIP 2020
% Sanjay Viswanath, Imaging Systems Lab, EE, IISc
%
% Function Output
% X: Dummy output from Multi-Resolution loop for initializing BCD iterations
% X2: Restored Image by COROSA after full set of BCD iterations
% cost: COROSA total cost from each iteration to check convergence
%
% Function Values
% TF: Transfer Function in Fourier domain, Sampling Trajectory for MRI
% Imn: Noisy image in restoration or complex Fourier measurement samples in MRI
% lamf: Lambda as tuning paramter for regularization
% beta: Paramter tau in spatial weight penalty
% betaadmm: Parameter gamma in ADMM 
% NIter: Number of ADMM iterations in each stage
% N: Number of CG iterations for evaluating restored image
% Nf: Number of BCD iterations
% L: Number of Multi-Resolution levels
% UB: Upperbound for image values, Effective bounding constraint is [0,UB]
% soregtype: Second order term, 1 for Hessian-Schattern, 2 for TV2
% noisetype: 1 for AWGN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, X2, cost] = corosa(TF, Imn, lamf, beta, betaadmm, ...
    NIter, N,  Nf, L,  UB, soregtype, noisetype)

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
betaadmm = gpuArray(betaadmm);
NIter = gpuArray(NIter);
N = gpuArray(N);
Nf = gpuArray(Nf);
L = gpuArray(L);
UB = gpuArray(UB);


Sz = size(Imn);

% Tikhonov Initialization
B = ifft2(fft2(Imn).*conj(TF));
Q = real(ifft2(conj(TF).*TF)) + lamf*(GetQQr2(Sz));

% Downsampling to MR coarsest level 'L'
B = fft2(Reduce(B,L));
Q = fft2(Reduce(Q,L));
X = real(ifft2(B./Q));

% HS regularization at level 'L'
[X, cost{1}]  = admm_sr(X, TF, Imn, lamf, ...
    soregtype, betaadmm, L, UB, NIter, noisetype);

% Iterations from MR levels L to 1
for ind = 1:L
    
    %     ind
    Xe = Expand(X,L-ind+1);
    R1 = sqrt(ApplyDerSqiso1(Xe));
    if soregtype == 2
        R2 = sqrt(ApplyDerSqiso2(Xe)); % If TV2 is selected as 2nd order term
    elseif soregtype == 1
        R2 = GetHS1(Xe);    % If HS is selected as 2nd order term
    end
    
    % Spatial Tau 
    tmp = (Xe.*(Xe > 0));
    tmp = tmp/max(tmp(:));
    betai = exp(-100*tmp.^2);
    betai = betai-min(betai(:));
    betai = betai+(0.01*max(betai(:))/100);
    betai = beta*betai*100/max(betai(:));
    
    % Relative Weight Beta evaluation at level 'ind'
    A = GetRelWt(R1, R2, betai);
    Wr1 =  A;
    Wr2 =  (1-A);
    
    X = Expand(X,1);
    
    % Evaluating restored image with ADMM to initialize level 'ind-1' 
    [X, cost{ind+1}]  = cov_admm_sr(X, TF, Imn,  Wr1, Wr2, lamf, ...
        soregtype, betaadmm, L-ind, UB, NIter, N, noisetype);
    
end

% MR initialization for BCD iterations
X2 = X;

% BCD iterations on final grid
for ind = 1:Nf
    
    
    R1 = sqrt(ApplyDerSqiso1(X2));
    if soregtype == 2
        R2 = sqrt(ApplyDerSqiso2(X2)); % If TV2 is selected as 2nd order term
    elseif soregtype == 1
        R2 = GetHS1(X2);    % If HS is selected as 2nd order term
    end
    
    % Spatial Tau
    tmp = (X2.*(X2 > 0));
    tmp = tmp/max(tmp(:));
    betai = exp(-100*tmp.^2);
    betai = betai-min(betai(:));
    betai = betai+(0.01*max(betai(:))/100);
    betai = betai*100/max(betai(:));
    
    % Relative Weight Beta evaluation
    A = GetRelWt(R1, R2, betai);
    Wr1 =  A;
    Wr2 =  (1-A);
    
    % Restored Image X2 evaluation with ADMM
    [X2, cost{L+1+ind}]  = cov_admm_sr(X2, TF, Imn, Wr1, Wr2, lamf, ...
        soregtype, betaadmm, 0, UB, NIter, N, noisetype);
    
    
    
end

X = gather(X);
X2 = gather(X2);
cost = gather(cost);


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [X, cost]  = admm_sr(X, TF, Imn, lamf, ...
    soregtype, beta, L, UB, NIter, noisetype)


Sz = size(Imn);
tmp = zeros(Sz,'gpuArray'); tmp(1) = 1;
Q = real(ifft2(conj(TF).*TF)) +  tmp +  GetQQr2(Sz);
Q = fft2(ReduceQ(Q,L)); clear tmp;

lamcm  = zeros(size(Imn),'gpuArray');
lamcb  = zeros(size(Imn),'gpuArray');
lamcr2 = zeros([size(Imn) 3],'gpuArray');

Rb = Expand(X,L);
Rm = (ifft2(TF.*fft2(Rb)));
R2 = GetD2(Rb);

for ind = 1:NIter
    
    
    Zb  = proxb(Rb, lamcb,  beta, UB);
    Zm  = proxm(Rm, lamcm,  Imn,  beta, noisetype);
    Zr2 = prox2(R2, lamcr2, lamf, beta,  soregtype);
    
    
    B =   Zb  - lamcb/beta + ...
        real(ifft2(conj(TF).*fft2(Zm  - lamcm/beta))) + ...
        Div2(Zr2-lamcr2/beta);
    B =   Reduce(B,L);
    
    
    X = real(ifft2(fft2(B)./Q));
    
    
    Rb = Expand(X,L);
    Rm = (ifft2(TF.*fft2(Rb)));
    R2 = GetD2(Rb);
    
    
    lamcb  = lamcb  + beta*(Rb - Zb);
    lamcm  = lamcm  + beta*(Rm - Zm);
    lamcr2 = lamcr2 + beta*(R2 - Zr2);
    
    cost(ind) = CostNLDecon ( X,  TF, Imn, lamf, ...
        soregtype, L );
    
    
    
end


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [X, cost]  = cov_admm_sr(X, TF, Imn, Wr1, Wr2, lamf, ...
    soregtype, beta, L, UB, NIter, N, noisetype)


cst = sum(abs(Imn(:)).^2);

lamcb  = zeros(size(Imn),'gpuArray');
lamcm  = zeros(size(Imn),'gpuArray');
lamcr1 = zeros([size(Imn) 2],'gpuArray');
lamcr2 = zeros([size(Imn) 3],'gpuArray');

Rb = Expand(X,L);
Rm = (ifft2(TF.*fft2(Rb)));
R1 = bsxfun(@times, GetD1(Rb), Wr1);
R2 = bsxfun(@times, GetD2(Rb), Wr2);



for ind = 1:NIter
    
    
    Zb  = proxb(Rb, lamcb,  beta, UB);
    Zm  = proxm(Rm, lamcm,  Imn,  beta, noisetype);
    Zr1 = prox1(R1, lamcr1, lamf, beta);
    Zr2 = prox2(R2, lamcr2, lamf, beta, soregtype);
    
    
    B =   Zb  - lamcb/beta + ...
        real(ifft2(conj(TF).*fft2(Zm  - lamcm/beta))) + ...
        Div1(bsxfun(@times, Zr1-lamcr1/beta, Wr1)) + ...
        Div2(bsxfun(@times, Zr2-lamcr2/beta, Wr2));
    B =   Reduce(B,L);
    
    
    X = MinQ(X, TF, Wr1, Wr2, B, L,   N, cst);
    
    
    
    Rb = Expand(X,L);
    Rm = (ifft2(TF.*fft2(Rb)));
    R1 = bsxfun(@times, GetD1(Rb), Wr1);
    R2 = bsxfun(@times, GetD2(Rb), Wr2);
    
    lamcb  = lamcb  + beta*(Rb - Zb);
    lamcm  = lamcm  + beta*(Rm - Zm);
    lamcr1 = lamcr1 + beta*(R1 - Zr1);
    lamcr2 = lamcr2 + beta*(R2 - Zr2);
    
    
    
    cost(ind) = CostNLDeconf ( X, TF, Imn,  lamf,  Wr1, Wr2, ...
        soregtype,  L );
    
    
    
end


return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function X = MinQ(X, TF,  Wr1, Wr2, B, L,   N, cst)


grad = ApplyM(X, TF,  Wr1, Wr2,  L) - B;
gradp = grad;
gamma = sum(grad(:).*gradp(:));

qcost = X.*(grad - B);
qcost = sum(qcost(:)) + cst;

dir = - gradp;


for ind = 1:N
    
    dirm = ApplyM( dir, TF,  Wr1, Wr2,  L);
    kappa =  sum(dirm(:).*dir(:));
    
    if abs(kappa) < 1e-8
        break;
    end
    
    alpha = gamma/kappa;
    
    X = X + alpha*dir;
    grad = grad + alpha*dirm;
    
    if norm(alpha*dir(:))/norm(X(:)) < 1e-6
        break;
    end
    
    qcost = X.*(grad - B);
    qcost = sum(qcost(:)) + cst;
    
    
    gradp = grad;
    
    gamma2 = sum(grad(:).*gradp(:));
    
    if abs(gamma) < 1e-8
        break;
    end
    
    beta = gamma2/gamma;
    
    
    
    dir = -gradp +  beta*dir;
    
    gamma = gamma2;
    
end


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function R = ApplyM(X, TF,  Wr1, Wr2,  L)

Xe = Expand(X,L);

R  = Reduce(real(ifft2(conj(TF).*TF.*fft2(Xe))), L);

R  = R + Reduce(Xe, L)     + ...
    Reduce(Div1(bsxfun(@times, GetD1(Xe), Wr1.^2)), L) + ...
    Reduce(Div2(bsxfun(@times, GetD2(Xe), Wr2.^2)), L);



return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Xc = Expand(X, L)

Xc = X;
Szc = size(Xc);

for ind = 1:L
    
    p1 = 6*Xc + Xc([end 1:end-1],:,:) + Xc([2:end 1],:,:);
    p2 = 4*(Xc + Xc([2:end 1],:,:));
    Szc = Szc.*[2 1];
    Xc = zeros(Szc,'gpuArray');
    Xc(1:2:end-1,:,:) = p1/8;
    Xc(2:2:end,:,:)  = p2/8;
    
end


for ind = 1:L
    
    p1 = 6*Xc + Xc(:,[end 1:end-1],:) + Xc(:,[2:end 1],:);
    p2 = 4*(Xc + Xc(:,[2:end 1],:));
    Szc = Szc.*[1 2];
    Xc = zeros(Szc,'gpuArray');
    Xc(:,1:2:end-1,:) = p1/8;
    Xc(:,2:2:end,:)  = p2/8;
    
end


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Xc = Reduce(X, L)

Xc = X;


for ind = 1:L
    
    p1 = Xc(1:2:end-1,:,:);
    p2 = Xc(2:2:end,  :,:);
    
    p1 = 6*p1 + p1([end 1:end-1],:,:) + p1([2:end 1],:,:);
    p2 = 4*(p2 + p2([end 1:end-1],:,:));
    
    Xc = (p1 + p2)/8;
    
end


for ind = 1:L
    
    p1 = Xc(:,1:2:end-1,:);
    p2 = Xc(:,2:2:end,  :);
    
    p1 = 6*p1 + p1(:,[end 1:end-1],:) + p1(:,[2:end 1],:);
    p2 = 4*(p2 + p2(:,[end 1:end-1],:));
    
    Xc = (p1 + p2)/8;
    
end


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function Qc = ReduceQ(Q, L)

tmp = zeros(size(Q)./[2^L 2^L],'gpuArray');
tmp(1) = 1;
Edel = fftn(Expand(tmp, L));
Qc = ifftn(conj(Edel).*Q.*Edel);
Qc = fftn(Qc(1:2^L:end,1:2^L:end,:));

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function cost = CostNLDecon ( X, TF, Imn, lamf,  soregtype,  L )


Xf = Expand(X, L);

cost = abs((ifft2(TF.*fft2(Xf))) - Imn).^2;


if soregtype == 2
    R = (ApplyDerSqiso2(Xf)).^0.5;
else
    R = GetHS1( Xf );
end

cost = sum(cost(:)) + lamf*sum(R(:));


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function cost = CostNLDeconf ( X,  TF, Imn,  lamf,  Wr1, Wr2, ...
    soregtype,  L )


Xf = Expand(X, L);

cost = abs((ifft2(TF.*fft2(Xf))) - Imn).^2;

if soregtype == 1
    R2 = sqrt(ApplyDerSqiso2(Xf));
else
    R2 = GetHS1( Xf );
end
R1 = sqrt(ApplyDerSqiso1(Xf));

R =  Wr1.*sqrt(R1) + Wr2.*R2;

cost = sum(cost(:)) + lamf*sum(R(:));


return;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function R = ApplyDerSqiso2(X)



R  = (-2*X + X(:,[end 1:end-1]) + X(:,[2:end 1])).^2;

R  = R + (-2*X + X([end 1:end-1],:) + X([2:end 1],:)).^2;

R  = R +  2*(   X + X([end 1:end-1],[end 1:end-1])  ...
    - X([end 1:end-1],:) - X(:,[end 1:end-1]) ).^2;



return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function R = ApplyDerSqiso1(X)



R  = (X - X(:,[end 1:end-1])).^2;
R  = R + (X - X([end 1:end-1],:)).^2;




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

X = zeros(Sz,'gpuArray');
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


function D1 = GetD1( X )

D1(:,:,1) = X - X(:,[end 1:end-1],:);
D1(:,:,2) = X - X([end 1:end-1],:,:);


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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Zr1 = prox1(D1, lamcr1, lamf, beta)

D1 = D1 + lamcr1/beta;

r = sqrt(D1(:,:,1).^2 + D1(:,:,2).^2);
rn = r - lamf/beta;
rn = rn.*(rn >= 0);
r  = r + (r < 1e-12);
Zr1 = bsxfun(@times, D1, rn./r);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Zr2 = prox2(D2, lamcr2, lamf, beta, soregtype)

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function X = Div1(D1)


X =  D1(:,:,1) -  D1(:,[2:end 1],1);

X  =  X + D1(:,:,2) -  D1([2:end 1],:,2);



return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [A,ErrC] = GetRelWt(R1, R2, beta)

 
d = R1 - R2;
 
Des = sqrt( (d+2*beta).^2 - 4*(d.*beta) );

C0 = ((d+2*beta) + Des)./(2*d);
C1 = ((d+2*beta) - Des)./(2*d);

IND0 = (C0 <=1) & (C0 >= 0);
IND1 = ( (C1 <=1) & (C1 >= 0) );

CHK = IND0 | IND1;
if min(CHK(:)) == 0
   ErrC = 1;
else
   ErrC = 0;
end


A = IND0.*C0 + IND1.*C1.*(~IND0);


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
