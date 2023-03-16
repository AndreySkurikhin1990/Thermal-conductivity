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

