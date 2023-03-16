function Semigrxchdf
% Program Semigrxchdf

N = 2;
% Dimensions
L1       = 0.8;
L2       = 0.6;

% Surface 1 (bottom)
A(1)     = L1;                % per unit depth
HO(1)    = 866;               % solar irradiation including specular 
                              % reflection from A2
EPS(1,1) = 0.8;
EPS(2,1) = 0.1;
id(1)    = 1;                 % T specified
T(1)     = 350.;              % T in K
% Surface2 (left)
A(2)     = L2;                % per unit depth
HO(2)    = 500.;              % direct solar irradiation
EPS(1,2) = 0.1;
EPS(2,2) = 0.8;
id(2)=0;                      % q specified
q(2)= 0.;                     % q in W/m2
% View Factors; since configuration is open (iclsd=2), 
% diagonal terms are also needed
iclsd    = 2;
% Range 1&2 are the same in this problem (no specular reflected 
% energy hits another surface)
F(1,1)   = 0.;
F(1,2)   = 0.25;
F(2,2)   = 0.;
% Solve system of equations
[q2, T2]   = SEMIGRAYDF(iclsd, A, EPS, HO, F, id, q, T);


% Output
fprintf('\n surface       T             q [W/m2]         Q [W] \n');
for i = 1:N
    fprintf('%5d %14f %14f %14f \n', i, T2(i), q2(i), q2(i)*A(i));
end




%==========================================================================
%
%==========================================================================
function [q, T]  =  SEMIGRAYDF(iclsd, A, EPS, HO, F, ID, q, T)
% Function Semigraydf
%******************************************************************************
% Routine to solve for missing surface temperatures/fluxes for the special case that
% all surfaces are diffuse reflectors; calling GRAYDIFF for 2 spectral ranges
% INPUT:
% N       =  number of surfaces in enclosure
% iclsd   =  closed or open configuration identifier
%         iclsd = 1: configuration is closed; diagonal Fii evaluated from summation rule
%         iclsd/ = 1: configuration has openings; Fii must be specified
% A(N)    =  vector containing surface areas, (m2)
% EPS(2,N) =  matrix containing surface emittances for 2 spectral ranges
% HO(N)   =  vector containing external irradiation, in (W/m2)
% F(N,N)  =  matrix containing view factors; on input only Fsij with j>i (iclsd = 1) or
%          j> = i (iclsd/ = 1) are required; remainder are calculated
% ID(N)   =  vector containing surface identifier:
%         ID = 0: surface heat flux is specified, in (W/m2)
%         ID = 1: surface temperature is specified, in (K)
% q(N)    =  vector containing known surface fluxes (only for surfaces with id = 2)
% T(N)    =  vector containing known surface temperatures (only for surfaces with id = 1)
% OUTPUT:
% q(N)    =  vector containing surface fluxes (for all surfaces)
% T(N)    =  vector containing surface temperatures (for all surfaces)


sigma  = 5.670E-8;
N = size(EPS,2);
% First calculate fluxes for external irradiation (range 1)
for i  =  1:N
    id1(i) = 1;   %#ok<*AGROW> % all surfaces have T = 0
    PIN(i) = 0.0;
end

% Set local values for eps and Fs
for i  =  1:N
    epsl(i) = EPS(1,i);
end
q1  =  GRAYDIFF(iclsd,A,epsl,HO,F,id1,PIN);

% Calculate unknown fluxes and temperatures by setting q2 = q-q1;
% the q1 terms give rise to a "negative irradiation" term
for i = 1:N
    % Set local values for eps
    epsl(i) = EPS(2,i);
    HO2(i)  = -q1(i)/EPS(2,i);
    for j = 1:N
        HO2(i) = HO2(i)+(1./EPS(2,j)-1.)*F(i,j)*q1(j);
    end
    % Fill PIN array with q and T
    if(ID(i) == 0)
        PIN(i) = q(i);
    else
        PIN(i) = sigma*pow(T(i),4.0);    % Convert temperatures to emissive powers
    end
end
POUT = GRAYDIFF(iclsd,A,epsl,HO2,F,ID,PIN);

% Convert emissive powers to temperatures
for i  =  1:N
    if(ID(i) == 0)
        T(i) = pow((POUT(i)/sigma),.25);
    else
        q(i) = POUT(i);
    end
end

%__________________________________________________________________________
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
        F(i,j) = A(j)/A(i)*F(j, i);
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

%__________________________________________________________________________            
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


%__________________________________________________________________________
function y = pow(x, n)
y = x.^n;

function [d] = KronD(j,k)

if nargin < 2
    error('Too few inputs.  See help KronD')
elseif nargin > 2
    error('Too many inputs.  See help KronD')
end

if j == k
    d = 1;
else
    d = 0;
end
