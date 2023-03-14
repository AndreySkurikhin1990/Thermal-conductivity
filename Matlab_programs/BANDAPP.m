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



