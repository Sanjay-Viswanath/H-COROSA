% Second Order Total Generalized Variation, TGV2 Smooth Version
% Imaging System Lab, EE, IISc
%
% Function Output
% X: Restored Image by TGV2 after full set of ADMM iterations
% cost: TV2 total cost from each iteration to check convergence
%
% Function Values
% Imn: Noisy image in restoration or complex Fourier measurement samples in MRI
% TF: Transfer Function in Fourier domain, Sampling Trajectory for MRI
% alpha1: Spatially varying weight for first order term
% alpha2: Spatially varying weight for second order term
% lam1: Tuning parameter for first order term
% lam2: Tuning parameter for second order term
% eps1: Smooth approximation paramter for TV functional 
% Ni: Number of inner CG iterations
% No: Number of TGV2 iterations 

function [X,cost,costq] = tgv2(Imn, TF, alpha1, alpha2, lam1, lam2, epsl, Ni, No)


[X, B, PC] = InitTGVapp(Imn, TF, lam1, lam2);

cst = sum(Imn(:).^2);

for ind = 1:No
    
    cost(ind) = CostTGVapp(Imn, TF, X, alpha1, alpha2, lam1, lam2, epsl);
    
    [W1, W2] = GetW(X,epsl);
    [Xn,costq] = MMCG(X, TF, B, lam1*alpha1.*W1, lam2*alpha2.*W2, PC, Ni, cst);
    X  = Xn;
    
end


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [X, B, PC] = InitTGVapp(Imn, TF, lam1, lam2)

Sz = size(Imn);

B = conj(TF).*fft2(Imn);

tmp = zeros([Sz 3]);
tmp(1,1,1) = 1;
PC(:,:,:,1) = ApplyM(tmp, TF, lam1*ones(Sz), lam2*ones(Sz));

tmp = zeros([Sz 3]);
tmp(1,1,2) = 1;
PC(:,:,:,2) = ApplyM(tmp, TF, lam1*ones(Sz), lam2*ones(Sz));

tmp = zeros([Sz 3]);
tmp(1,1,3) = 1;
PC(:,:,:,3) = ApplyM(tmp, TF, lam1*ones(Sz), lam2*ones(Sz));


for ind = 1:3
    for ind2 = 1:3
        M(:,:,ind,ind2) =  fft2(PC(:,:,ind,ind2));
    end
end


for ind = 1:Sz(1)
    for ind2 = 1:Sz(2)
        PC(ind,ind2,:,:) = inv(squeeze(M(ind,ind2,:,:)));
    end
end


[Xtv2, ~] = tv2(TF, Imn, lam2, 1,  ...
                  250,   300, 2, 1);
X(:,:,1) = Xtv2;
X(:,:,2) = Xtv2 - Xtv2(:,[end 1:end-1]);
X(:,:,3) = Xtv2 - Xtv2([end 1:end-1],:);

B = real(ifft2(B));
B(:,:,2:3) = zeros([Sz 2]);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [W1, W2] = GetW(X, epsl)

tmp = (X(:,:,1) - X(:,[end 1:end-1],1) - X(:,:,2)).^2 + ...
            (X(:,:,1) - X([end 1:end-1],:,1) - X(:,:,3)).^2;
W1 = 0.5./sqrt(epsl+tmp);

tmp = (X(:,:,2) - X(:,[end 1:end-1],2)).^2 + ...
      (X(:,:,3) - X([end 1:end-1],:,3)).^2 + ...
      0.5*( X(:,:,2) - X([end 1:end-1],:,2) + X(:,:,3) - X(:,[end 1:end-1],3)).^2;

W2 = 0.5./sqrt(epsl+tmp);


return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
function [X,qcosta] = MMCG (X, TF, B, W1, W2, PC, Ni, cst)      

   grad = ApplyM(X, TF, W1, W2) - B;
   gradp = grad; 
   gamma = sum(grad(:).*gradp(:));

   qcost = X.*(grad - B);
   qcosta(1) = sum(qcost(:)) + cst;
   
   dir = - gradp;
   
   
   for ind = 1:Ni
       
       dirm = ApplyM( dir, TF, W1, W2 );
       kappa =  sum(dirm(:).*dir(:));
       
       if abs(kappa) < 1e-12
           break;
       end
       
       alpha = gamma/kappa;
       
       X = X + alpha*dir;
       grad = grad + alpha*dirm;
       
       if norm(alpha*dir(:))/norm(X(:)) < 1e-6
           break;
       end
       
       qcost = X.*(grad - B);
       qcosta(1+ind) = sum(qcost(:)) + cst;

       
       gradp = grad; 
       %gradp  = ApplyPC(PC,  grad);
       gamma2 = sum(grad(:).*gradp(:));
       
       if abs(gamma) < 1e-12
           break;
       end
       
       beta = gamma2/gamma;
       
      
       
       dir = -gradp +  beta*dir;
       
       gamma = gamma2;
   
  end
   
   
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

function  R = ApplyM(X, TF, W1, W2)

R = real(ifft2(conj(TF).*TF.*fft2(X(:,:,1))));
R(:,:,2:3) = zeros(size(X(:,:,2:3)));

tmp = W1.*(X(:,:,1) - X(:,[end 1:end-1],1) - X(:,:,2));
R(:,:,1) = R(:,:,1) +  tmp - tmp(:,[2:end 1]);
R(:,:,2) = R(:,:,2)  - tmp;


tmp = W1.*(X(:,:,1) - X([end 1:end-1],:,1) - X(:,:,3));
R(:,:,1) = R(:,:,1) + tmp - tmp([2:end 1],:);
R(:,:,3) = R(:,:,3)  - tmp;

tmp = W2.*(X(:,:,2) - X(:,[end 1:end-1],2));
R(:,:,2) = R(:,:,2) +  tmp - tmp(:,[2:end 1]);

tmp = W2.*(X(:,:,3) - X([end 1:end-1],:,3));
R(:,:,3) = R(:,:,3) +  tmp - tmp([2:end 1],:);


tmp = 0.5*( X(:,:,2) - X([end 1:end-1],:,2) + X(:,:,3) - X(:,[end 1:end-1],3));
tmp = W2.*tmp;
R(:,:,2) = R(:,:,2) + (tmp - tmp([2:end 1],:));
R(:,:,3) = R(:,:,3) + (tmp - tmp(:,[2:end 1]));


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function R = ApplyPC(PC,X)

X = fft(fft(X,[],1),[],2);

for ind = 1:3
  R(:,:,ind) = PC(:,:,ind,1).*X(:,:,1) + ...
      PC(:,:,ind,2).*X(:,:,2) + PC(:,:,ind,3).*X(:,:,3);
end

R = real(ifft(ifft(R,[],1),[],2));

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function cost = CostTGVapp(Imn, TF, X, alpha1, alpha2, lam1, lam2, epsl)

Wd = abs((ifft2(fft2(X).*TF)) - Imn).^2;

tmp = (X(:,:,1) - X(:,[end 1:end-1],1) - X(:,:,2)).^2 + ...
            (X(:,:,1) - X([end 1:end-1],:,1) - X(:,:,3)).^2;
W1 = alpha1.*sqrt(epsl+tmp);

tmp = (X(:,:,2) - X(:,[end 1:end-1],2)).^2 + ...
      (X(:,:,3) - X([end 1:end-1],:,3)).^2 + ...
      0.5*( X(:,:,2) - X([end 1:end-1],:,2) + X(:,:,3) - X(:,[end 1:end-1],3)).^2;

W2 = alpha2.*sqrt(epsl+tmp);

cost = sum(Wd(:)) + lam1*sum(W1(:)) + lam2*sum(W2(:));

return;
    
