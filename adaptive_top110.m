clear;clc;close all

nelx = 120; nely = 40; volfrac = 0.5; penal = 3.0; rmin = 5; ft = 3.; eta = 0.5; s=0;
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
beta = 1;
if ft == 1 || ft == 2
    xPhys = x;
elseif ft == 3
    xTilde = x;
    xPhys = (tanh(beta*eta)+tanh(beta*(xTilde-eta)))/(tanh(beta*eta)+tanh(beta*(1-eta)));    %% Tangent hyperbolic projection
end
loopbeta = 0;
loop = 0;
change = 1;
%% START ITERATION
while change > 0.01
    loopbeta = loopbeta+1;
    loop = loop+1;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dv = ones(nely,nelx);
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    elseif ft == 3
        dx = (beta*(sech(beta*(xTilde-eta))).^2)./(tanh(beta*eta)+tanh(beta*(1-eta)));    %% Chain Term for Tangent hyperbolic projection
        dc(:) = H*(dc(:).*dx(:)./Hs);
        dv(:) = H*(dv(:).*dx(:)./Hs);
    end
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    l1 = 0; l2 = 1e9; move = 0.2;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
        if ft == 1
            xPhys = xnew;
        elseif ft == 2
            xPhys(:) = (H*xnew(:))./Hs;
        elseif ft == 3
            xTilde(:) = (H*xnew(:))./Hs;
            xPhys = (tanh(beta*eta)+tanh(beta*(xTilde-eta)))/(tanh(beta*eta)+tanh(beta*(1-eta)));
        end
        if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
    end
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
        mean(xPhys(:)),change);
    %% PLOT DENSITIES
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ADAPTIVE BETA UPDATE SCHEME%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    betamax = 1000;  mndobj = 5; mndmin = 0.01;
    if ft == 3 && beta < 1000
        M_ndTilde=graylevelindicator(xTilde(:));
        
        M_ndPhys=graylevelindicator(xPhys(:));
        M_nd_init =graylevelindicator(volfrac);
        alpha = log(100)/log(M_nd_init);
        if M_ndPhys>mndobj
            beta = min(betamax,(M_nd_init/max(M_ndPhys,mndobj))^alpha); %% PHASE 1 UPDATE
            if loop>1 
                decrease_ratio=  (M_ndPhys)/M_ndPhys_old; %% REDUCTION RATE
                decrease_ratios(loop-1)=decrease_ratio; %% STORE REDUCTION RATE
                avg_decrease_ratio = mean(decrease_ratios); %% AVERAGE REDUCTION RATE
                if M_ndTilde_old-M_ndTilde<0.2 %% DETECT CONGETION 
                    s=s+1;
                end
            end
            if s>3 %% SWITCH TO PHASE 2
                beta_i = logspace(0,3,100); %% LOG-SPACED BETA SMAPLES
                xPhys_Predict = (tanh(beta_i'*eta)+tanh(beta_i'*(xTilde(:)'-eta)))./(tanh(beta_i'*eta)+tanh(beta*(1-eta))); %% PREDICTED DENSITY PROJECTION
                M_ndPhys_Predict=graylevelindicator(xPhys_Predict');%% PREDICTED GRAY LEVEL INDICATOR
                M_ndTraget=avg_decrease_ratio*M_ndPhys; %% TARGET GRAY LEVEL INDICATOR
                beta=min(betamax,interp1(M_ndPhys_Predict,beta_i,M_ndTraget)); %% PHASE 2 UPDATE
            end
        end
        M_ndPhys_old=M_ndPhys;
        M_ndTilde_old=M_ndTilde;
    end
    fprintf('Parameter beta increased to %g.\n',beta);
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mnd] =graylevelindicator(x)
mnd = (mean(x.*(1-x)))*400;
end