function Bandmxch
% Program Bandmxch

N = 2;
M = 2;

% Spectral band limit wavelengths
Bandwl(1) = 4.0;                    % wavelength in micron

% Dimensions
L1 = 0.8;
L2 = 0.6;


% Surface 1 (bottom)
A(1)       = L1;                         % per unit depth
HOs(1,1)   = 1203.5*bbfn(4.0*5762);
HOs(2,1)   = 1203.5 - HOs(1,1);                   % solar irradiation including specular reflection from A2
EPS(1,1)   = 0.8;
RHOs(1,1)  = 0.0;
EPS(2,1)   = 0.1;
RHOs(2,1)  = 0.0;
id(1)      = 1;                         % T specified
T(1)       = 350.;                       % T in K
% Surface2 (left)
A(2)      = L2;                         % per unit depth
HOs(1,2)  = 500.0*bbfn(4.*5762);
HOs(2,2)  = 500.0 - HOs(1,2);                     % direct solar irradiation
EPS(1,2)  = 0.1;
RHOs(1,2) = 0.9;
EPS(2,2)  = 0.8;
RHOs(2,2) = 0.2;
id(2)     = 0;                        % q specified
q(2)      = 0.;                        % q in W/m2
% View Factors; since configuration is open (iclsd=2), diagonal terms are also needed
iclsd = 2;
% Range 1&2 are the same in this problem (no specular reflected energy hits another surface)
Fs(1,1,1) = 0.;
Fs(1,1,2) = 0.25;
Fs(1,2,2) = 0.;
Fs(2,1,1) = 0.;
Fs(2,1,2) = 0.25;
Fs(2,2,2) = 0.;
% Solve system of equations
[q2, T2] = BANDAPP(iclsd, A, Bandwl, EPS, RHOs, HOs, Fs, id, q, T);
% Output
fprintf('\n surface       T             q [W/m2]         Q [W] \n');
for i = 1:N
    fprintf('%5d %14f %14f %14f \n', i, T2(i), q2(i), q2(i)*A(i));
end



%==========================================================================
%
%==========================================================================
function [q, T] = BANDAPP(iclsd, A, Bandwl, EPS, RHOs, HOs, Fs, id, q, T)
% Function Bandapp
%******************************************************************************
% Routine to solve for missing surface temperatures/fluxes calling GRAYDIFSPEC
% for M spectral ranges
% INPUT:
% M      = number of spectral bands
% N      = number of surfaces in enclosure
% iclsd  = closed or open configuration identifier
%          iclsd=1: configuration is closed; diagonal Fsii evaluated from 
%          summation rule
%          iclsd/=1: configuration has openings; Fsii must be specified
% A(N)   = vector containing surface areas, (m2)
% Bandwl(M-1) = vector containing wavelengths separating the M spectral 
%               bands, (micro-m)
% EPS(M,N) = matrix containing surface emittances for 2 spectral ranges
% RHOs(M,N)= vector containing surface specular reflectance components
% HOs(M,N) = vector containing external irradiation, in (W/m2)
% Fs(M,N,N)= matrix containing view factors; on input only Fsij with 
%           j>i (iclsd=1) or
%           j>=i (iclsd/=1) are required; remainder are calculated
% ID(N)  = vector containing surface identifier:
%          ID=0: surface heat flux is specified, in (W/m2)
%          ID=1: surface temperature is specified, in (K)
% q(N)   = vector containing known surface fluxes 
%           (only for surfaces with id=2)
% T(N)   = vector containing known surface temperatures 
%          (only for surfaces with id=1)
% OUTPUT:
% q(N)   = vector containing surface fluxes (for all surfaces)
% T(N)   = vector containing surface temperatures (for all surfaces)

sigma = 5.670E-8;
Tresset = 0.1;     %temperature residual setpoint in K

% First guess temperature for surfaces for which flux is specified
Tmin=1.e8;
Tmax=0.;
iqsp=0;
N = size(EPS,2);
M = size(EPS,1);
for i = 1:N
    if(id(i)~=0)
        if(T(i)>Tmax)
            Tmax=T(i);
        end
        if(T(i)<Tmin)
            Tmin=T(i);
        end
    end
end
for i = 1:N
    if(id(i)~=1)
        iqsp=iqsp+1;             % number of surfaces with q specified
        T(i)=.5*(Tmin+Tmax); %#ok<*SAGROW>
    end
end
Bandwl(M)=1.e5;
% Scan over spectral bands
counter0 = 0;
Tres = 10;
while(Tres>Tresset)
    for i = 1:N
        fracl(i)=0.;
    end
    % Initialize (calculated) total fluxes and emissive powers
    for i = 1:N
        qc(i)=0.;
    end
    for i = 1:N
        Eb(i)=sigma*pow(T(i),4.0);
    end
    for j = 1:M
        % Set local values for eps, rhos, hos, pin and Fs
        for i = 1:N
            epsl(i)=EPS(j,i);
            rhosl(i)=RHOs(j,i);
            hosl(i)=HOs(j,i);
            frach(i) = bbfn(Bandwl(j)*T(i));  % bbfn at max wavelength for band
            PIN(i)=Eb(i)*(frach(i)-fracl(i)); % fraction of emissive power in band
            for k = i:N
                FSl(i,k)=Fs(j,i,k);
            end
        end
        for o = 1:N
            IDl(o)=1;        % as used here ALL T's are specified
            fracl(o)=frach(o);
        end                % bbfn at min wavelength for (next) band
        
        POUT = GRAYDIFFSPEC(iclsd,A,epsl,rhosl,hosl,FSl,IDl,PIN);
        % Adding fluxes of all bands
        for o = 1:N
            qc(o)=qc(o)+POUT(o);
        end
    end
    
    % Update temperatures for surfaces for which flux is specified
    % This is done by assuming irradiation H(r) to be correct in eq. 
    % (6.10-top; modified
    % for bandmodel), ie, FZ= (q - eps*Eb)_last - (q_specified - eps*Eb(T_new))
    % (this REALLY should be in double precision!)
    Tres=0.;             % initialize temperature error residual
    for i = 1:N
        if(id(i)~=1)
            % Calculate unchangeable part of zero function 
            % (all except eps*Eb(T_new))
            FZb=qc(i)-q(i);
            ebfrl=0;
            for j = 1:M
                ebfrh=bbfn(Bandwl(j)*T(i));
                FZb=FZb-EPS(j,i)*Eb(i)*(ebfrh-ebfrl);
                ebfrl=ebfrh;
            end
            
            if(qc(i)>q(i))
                % Temperature is too high
                Ttop=T(i);
                Tbot=Ttop;
                counter = 0;
                FZ=10;
                
                while (FZ>0.)
                    Tbot=Tbot-10.;
                    FZ=FZb;
                    
                    ebfrl=0;
                    Ebbot=sigma*pow(Tbot,4.0);
                    for j = 1:M
                        ebfrh=bbfn(Bandwl(j)*Tbot);
                        FZ=FZ+EPS(j,i)*Ebbot*(ebfrh-ebfrl);
                        ebfrl=ebfrh;
                    end
                    if (counter > 100)
                        break;
                    end
                    counter = counter + 1;
                end
                 % not low enough yet
            else
                % Temperature is too low
                Tbot=T(i);
                Ttop=Tbot;
                FZ=FZb;
                counter = 0;
                while(FZ<0.)
                    Ttop=Ttop+10.;
                    % Calculate remainder of zero function
                    FZ=FZb;
                    ebfrl=0;
                    Ebtop=sigma*pow(Ttop,4.0);
                    for j = 1:M
                        ebfrh=bbfn(Bandwl(j)*Ttop);
                        FZ=FZ+EPS(j,i)*Ebtop*(ebfrh-ebfrl);
                        ebfrl=ebfrh;
                    end
                    if (counter > 100)
                        break;
                    end
                    counter = counter + 1;
                end
            end   % not high enough yet
            counter = 0;
            Tnew=0.5*(Ttop+Tbot);
            while (Ttop-Tbot>Tresset)
                Tnew=0.5*(Ttop+Tbot);
                % Calculate remainder of zero function
                FZ=FZb;
                ebfrl=0;
                Ebnew=sigma*pow(Tnew,4.0);
                for j = 1:M
                    ebfrh=bbfn(Bandwl(j)*Tnew);
                    FZ=FZ+EPS(j,i)*Ebnew*(ebfrh-ebfrl);
                    ebfrl=ebfrh;
                end
                % Check for convergence (really need double precision here!)
                if(FZ<0.)
                    Tbot=Tnew;   % too low
                end
                if(FZ>0.)
                    Ttop=Tnew;   % too high
                end
                counter = counter + 1;
                if (counter > 100)
                    break;
                end
            end
            Tres=Tres+abs(T(i)-Tnew);
            T(i)=Tnew;
        end
    end % i
    Tres=Tres/iqsp;
    
    if (counter0 > 100)
        break;
    end
    counter0 = counter0 + 1;
end
for o = 1:N
    q(o)=qc(o);
end


%__________________________________________________________________________
function [BBFN] =  bbfn(X)
%    ********************************************************************
%     *  This subroutine calculates the fractional blackbody             *
%     *  emissive power f(n*lambda*T), where X=n*lambda*T in (micro-m*K) *
%     ********************************************************************

CC = 15/pi^4;           % [-]
C2 = 1.4388e4;            % [cm K]
eps = 1e-16;            %


%
V = C2/X;               % integration lower limit
EX = exp(V);
%
M = 1;
BBFN = 0;
EM = 1;
VM = M*V;
while (VM^3*EM > eps)
    BM = (6.0 + VM*(6.0 + VM*(3.0 + VM)))/M^4;
    EM = EM/EX;
    BBFN = BBFN + BM*EM;
    M = M + 1;
    VM = M*V;
    if M > 1000000
        fprintf('M = %s \n', num2str(M));
        break;
    end
end
BBFN = CC*BBFN;

%__________________________________________________________________________
function [POUT] =  GRAYDIFFSPEC(iclsd,A,EPS,RHOs,HOs,Fs,ID,PIN)
% Function graydifspec
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
% RHOs(N)= vector containing surface specular reflectance components
% HOs(N) = vector containing external irradiation, in [W/m2]
% Fs(N,N)= vector containing view factors; on input only Fsij with j>i (iclsd=1) or
%         j>=i (iclsd/=1) are required; remainder are calculated
% ID(N) = vector containing surface identifier:
%         ID=0: surface heat flux is specified, in [W/m2]
%         ID=1: surface temperature is specified, in [K]
% PIN(N)= vector containing surface emissive powers (id=1) and fluxes (id=2)
% OUTPUT:
% POUT(N)= vector containing unknown surface fluxes (for surfaces with id=1) and
%           temperatures (for surfaces with id=2)

% Compute missing view factors
% Lower left triangle by reciprocity
N = size(A,2);
for  i = 2:N
    for  j=1:(i-1)
        Fs(i,j) = A(j)/A(i)*Fs(j,i);
    end
end
% If closed configuration, need to also calculate diagonal terms by summation rule
if(iclsd == 1)
    for i = 1:N
        Fs(i,i) = 1.;
        for  j = 1:N
            if(j~=i)
                Fs(i,i) = Fs(i,i)-(1.-RHOs(j))*Fs(i,j);
            end
        end
        Fs(i,i) = Fs(i,i)/(1.-RHOs(i));
    end
end
% Fill q- and e-coefficient matrices
for  i = 1:N
    for  j = 1:N
        ikr = KronD(i,j);
        %ikr=(i/j)*(j/i);         % Kronecker delta_ij
        qm(i,j) = ikr/EPS(j) - ((1.-RHOs(j))/EPS(j)-1.)*Fs(i,j);
        em(i,j) = ikr - (1.-RHOs(j))*Fs(i,j);
    end
end

% Fill POUT-coefficient matrix and RHS
for  i = 1:N
    B(i) = -HOs(i);
    for  j = 1:N
        pm(i,j) = qm(i,j)*ID(j)-em(i,j)*(1-ID(j));
        B(i)    = B(i)+(em(i,j)*ID(j)-qm(i,j)*(1-ID(j)))*PIN(j);
    end
end

% Invert POUT-coefficient matrix and multiply by RHS to obtain POUT
POUT = GAUSS(pm,B);         %  can also do: POUT = pm\B', or inv(pm)*B';

%_________________________________________________________________________
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
function [y] = pow(x, n)
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


%__________________________________________________________________________






