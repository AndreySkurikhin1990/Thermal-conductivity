function [POUT] = GRAYDIFF(iclsd, A, EPS, HO, F, ID, PIN)
% Function graydiff
%
%******************************************************************************
% Routine to fill view factor matrix and solve for missing surface temperatures/fluxes
%
% INPUT:
% N     = number of surfaces in enclosure
% iclsd = closed or open configuration identifier
%         iclsd=1: configuration is closed; diagonal Fii evaluated from summation rule
%         iclsd/=1: configuration has openings; Fii must be specified
% A(N)  = vector containing surface areas, [m2]
% EPS(N)= vector containing surface emittances
% HO(N) = vector containing external irradiation, in [W/m2]
% F(N,N)= vector containing view factors; on input only Fij with j>i (iclsd=1) or
%         j>=i (iclsd/=1) are required; remainder are calculated
% ID(N) = vector containing surface identifier:
%         ID=0: surface heat flux is specified, in [W/m2]
%         ID=1: surface temperature is specified, in [K]
% PIN(N)= vector containing surface emissive powers (id=1) and fluxes (id=2)
% OUTPUT: 
% POUT(N)= vector containing unknown surface fluxes (for surfaces with id=1) and 
%           emissive powers (for surfaces with id=0)



% Compute missing view factors
% Lower left triangle by reciprocity
N = size(PIN,2);
for i = 2:N
    for j = 1:i-1
        F(i,j) = (A(j)/A(i))*F(j, i);
    end
end
% If closed configuration, need to also calculate diagonal terms by 
% summation rule
if(iclsd==1) 
    for i = 1:N
        F(i, i) = 1.;
        for  j = 1:N
            if i~=j
                F(i, i)=F(i,i) - F(i, j);
            end
        end
    end
end

% Fill q- and e-coefficient matrices
for  i = 1:N
    for  j = 1:N
        ikr = KronD(i,j);
        %ikr = (i/j)*(j/i)         % Kronecker delta_ij
        qm(i, j) = ikr/EPS(j)-(1./EPS(j)-1.)*F(i, j);
        em(i, j) = ikr - F(i, j);
    end
end
% Fill POUT-coefficient matrix and RHS
for i=1:N
    
    B(i) = -HO(i);
    for  j = 1:N
        pm(i, j) = qm(i, j)*ID(j) - em(i, j)*(1-ID(j));
        B(i) = B(i) +(em(i, j)*ID(j) - qm(i,j)*(1-ID(j)))*PIN(j);
    end
end

            
% Invert POUT-coefficient matrix and multiply by RHS to obtain POUT
POUT = GAUSS(pm, B);

%*******************************************************************
function [X] = GAUSS(A, B)
N = size(B,2);
for I=1:N
    L(I) = I;
    SMAX = 0.0;
    for J = 1:N
        if abs(A(I, J)) > SMAX
            SMAX = abs(A(I, J));
        end
    end
    S(I) = SMAX;
end
for  K=1:N-1 
     RMAX=0.0;
        for I = K:N 
        R = abs(A(L(I), K))/S(L(I));
        if (R<=RMAX)
            continue;
        else
            J = I;
            RMAX=R;
        end
        end
        LK = L(J);
        L(J)= L(K);
        L(K)=LK;
        for I = K+1:N
            XMULT = A(L(I), K)/A(LK, K);
            for  J = K+1:N 
                A(L(I), J) = A(L(I), J) - XMULT*A(LK, J);
            end
            A(L(I),K) = XMULT;
        end
end

for K = 1:N-1
    for I = K+1:N 
        B(L(I)) = B(L(I)) - A(L(I),K)*B(L(K));
    end
end
X(N) = B(L(N))/A(L(N), N);
for  I= N-1:-1:1
    SUM = B(L(I));
    for J = I+1:N
        SUM = SUM - A(L(I), J)*X(J);
    end
    X(I) = SUM/A(L(I), I);
end


function [d] = KronD(j,k)

if j == k
    d = 1;
else
    d = 0;
end