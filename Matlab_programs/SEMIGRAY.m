function [q, T] = SEMIGRAY(iclsd, A, EPS, RHOs, HOs, Fs, ID, q, T)
% function Semigray
%******************************************************************************
% Routine to solve for missing surface temperatures/fluxes calling GRAYDIFSPEC
% for 2 spectral ranges
% INPUT:
% N     = number of surfaces in enclosure
% iclsd = closed or open configuration identifier
%         iclsd=1: configuration is closed; diagonal Fii evaluated from 
%         summation rule
%         iclsd/=1: configuration has openings; Fii must be specified
% A(N)  = vector containing surface areas, [m2]
% EPS(2,N)= matrix containing surface emittances for 2 spectral ranges
% RHOs(2,N)= vector containing surface specular reflectance components
% HOs(N) = vector containing external irradiation, in [W/m2]
% Fs(2,N,N)= matrix containing view factors; on input only Fsij with j>i (iclsd=1) or
%          j>=i (iclsd/=1) are required; remainder are calculated
% ID(N) = vector containing surface identifier:
%         ID=0: surface heat flux is specified, in [W/m2]
%         ID=1: surface temperature is specified, in [K]
% q(N)   = vector containing known surface fluxes (only for surfaces with id=2)
% T(N)   = vector containing known surface temperatures (only for surfaces with id=1)
% OUTPUT:
% q(N)   = vector containing surface fluxes (for all surfaces)
% T(N)   = vector containing surface temperatures (for all surfaces)


sigma = 5.670E-8;

% First calculate fluxes for external irradiation (range 1)
N = size(EPS, 2);
for i = 1:N
    id1(i) = 1;   %#ok<*AGROW> % all surfaces have T=0
    PIN(i) = 0.0;
end

% Set local values for eps, rhos and Fs
for i = 1:N    
    epsl(i)   = EPS(1,i);
    rhosl(i)  = RHOs(1,i);
    for j = 1:N
        FSl(i,j) = Fs(1, i, j);
    end
end
q1 = GRAYDIFFSPEC(iclsd,A,epsl,rhosl,HOs,FSl,id1,PIN);
% Calculate unknown fluxes and temperatures by setting q2=q-q1;
% the q1 terms give rise to a "negative irradiation" term
for i = 1:N
    % Set local values for eps, rhos and Fs
    % Calculate missing viewfactors for determination of HOs2
    epsl(i)  = EPS(2,i);
    rhosl(i) = RHOs(2,i);
    for j = i:N
        FSl(i,j)  = Fs(2,i,j);
        FSl(j,i)  = A(i)/A(j)*FSl(i,j);
    end
    % If closed configuration, need to also calculate diagonal terms 
    % by summation rule
    if(iclsd == 1)
        FSl(i,i)  = 1.0;
        for j = 1:N
            if(j ~= i)
                FSl(i,i) = FSl(i, i) - (1.0 - rhosl(j))*FSl(i,j);
            end
        end
    end
    HOs2(i) = -q1(i)/EPS(2,i);
    for j = 1:N
        HOs2(i) = HOs2(i) + ((1. - rhosl(j))/EPS(2,j) - 1.0)*FSl(i,j)*q1(j);
    end
    % Fill PIN array with q and T
    if(ID(i) == 0)
        PIN(i) = q(i);
    else
         % Convert temperatures to emissive powers
        PIN(i) = sigma*pow(T(i),4.0); 
    end
end
POUT = GRAYDIFFSPEC(iclsd,A,epsl,rhosl,HOs2,FSl,ID,PIN);

% Convert emissive powers to temperatures
for i = 1:N
    if(ID(i) == 0)
        T(i) = pow((POUT(i)/sigma), 0.25);
    else
        q(i) = POUT(i);
    end
end



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

