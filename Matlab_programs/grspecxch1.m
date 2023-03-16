function t = grspecxch1()
%Program Grspecxch
format longg;
N = 6;
sigma=5.670E-8;

% Dimensions
ww=0.4;
hh=0.3;
ll=1.0e6;

% Surfaces 1&3 (bottom and top)
A(1)=ww*ll;
A(3)=ww*ll;          % ww x ll
HOs(1)=0.;
HOs(3)=0.;           % no external irradiation
EPS(1)=0.3;
RHOs(1)= 0.7;
EPS(3)=0.3;
RHOs(3)= 0.;
id(1)=1;
id(3)=1;             % T specified
T(1)=1000.;
T(3)=1000.;          % T in K

% Surfaces 2&4 (left and right)
A(2)=hh*ll;
A(4)=hh*ll;         % hh x ll
HOs(2)=0.;
HOs(4)=0.;          % no external irradiation
EPS(2)=0.8;
RHOs(2)= 0.2;
EPS(4)=0.8;
RHOs(4)= 0.;
id(2)=1;
id(4)=1;            % T specified! Fs(1,2)=F1-2
T(2)= 600.;
T(4)= 600.;         % T in K

% Surfaces 5&6 (front and back)
A(5)=hh*ww;
A(6)=hh*ww;         % hh x ww
HOs(5)=0.;
HOs(6)=0.;          % no external irradiation
EPS(5)=0.8;
RHOs(5)= 0.;
EPS(6)=0.8;
RHOs(6)= 0.;
id(5)=1;
id(6)=1;            % T specified
T(5)= 600.;         % Fs(1,2)=F1-2
T(6)= 600.;         % T in K

%fill PIN array with q ant T
for i=1:N
    if (id(i)==0)
        PIN(i)=q(i);
    else
        PIN(i)=sigma*T(i)^4;  %convert temperatures to emmissive powers
    end
end

% view factors; since configuration is closed (iclsd=1), diagonal terms are not needed
iclsd=1;

%F12
arg(1)=hh;
arg(2)=ll;
arg(3)=ww;		%h,l,w
Fs(1,2)=view(39,3,arg,0); %Fs(1,2)=F1-2

% Fs13
arg(1)=ww;
arg(2)=ll;
arg(3)=hh;     % a,b,c

r=view(38,3,arg,0);
t=parlplates(ww,ww,2.*ww,ll,0.,ll,hh);
b=RHOs(2);

Fs(1,3) =r+b*t;  % Fs1-3=F1-3+rhos2*F1(2)-3

b=RHOs(2);
t=perpplates(ww, 2.0*ww, 0.0, hh, ll, 0.0, ll);
Fs(1,4) = Fs(1,2)+b*t;         % Fs1-4=F1-4+rhos2*F1(2)-4
%Fs(1,4)=0;
Fs(1,5) = 0.5*(1.-(1.-RHOs(2))*Fs(1,2)-Fs(1,3)-Fs(1,4));
Fs(1,6) = Fs(1,5);

% Fs23
arg(1)=ww;
arg(2)=ll;
arg(3)=hh;    % h,l,w
Fs(2,3)=view(39,3,arg,0)+RHOs(1)*perpplates(hh,2.*hh,0.,ww,ll,0.,ll);  % Fs2-3=F2-3+rhos1*F1(1)-3
%Fs(2,3)=0;

% Fs24
arg(1)=ll;
arg(2)=hh;
arg(3)=ww;     % a,b,c
Fs(2,4)=view(38,3,arg,0)+RHOs(1)*parlplates(hh,hh,2.*hh,ll,0.,ll,ww);  % Fs2-4=F2-4+rhos1*F2(1)-4
%Fs(2,4)=0;

% Fs25
Fs(2,5)=.5*(1.-(1.-RHOs(1))*A(1)*Fs(1,2)/A(2)-Fs(2,3)-Fs(2,4));
Fs(2,6)=Fs(2,5);

% Fs33=rhos1*F3(1)-3 + rhos1*rhos2*F3(12+21)-3
arg(1)=ww;
arg(2)=ll;
arg(3)=2.0*hh;  % a,b,c
Fs(3,3)=RHOs(1)*(view(38,3,arg,0)+RHOs(2)*parlplates(ww,ww,2.*ww,ll,0.,ll,2.*hh));
%Fs(3,3)=0;

Fs(3,4)=Fs(1,4)+RHOs(1)*(perpplates(0.,ww,hh,2.*hh,ll,0.,ll) ...       % Fs3-4=Fs1-4+rhos1*F3(1)-3
+RHOs(2)*perpplates(ww,2.*ww,hh,2.*hh,ll,0.,ll));     %       +rhos1*rhos2*F3(12+21)-4
%Fs(3,4)=0;    

Fs(3,5)=.5*(1.-(1.-RHOs(1))*A(1)*Fs(1,3)/A(3)-(1.-RHOs(2))*A(2)*Fs(2,3)/A(3)-Fs(3,3)-Fs(3,4));
Fs(3,6)=Fs(3,5);

% Fs44=rhos2*F4(2)-4 + rhos1*rhos2*F4(12+21)-4
arg(1)=hh;
arg(2)=ll;
arg(3)=2.*ww;  % a,b,c
Fs(4,4)=RHOs(2)*(view(38,3,arg,0)+RHOs(1)*parlplates(hh,hh,2.*hh,ll,0.,ll,2.*ww));
%Fs(4,4)=0;

Fs(4,5)=.5*(1.-((1.-RHOs(1))*A(1)*Fs(1,4)+(1.-RHOs(2))*A(2)*Fs(2,4)+A(3)*Fs(3,4))/A(4)-Fs(4,4));
Fs(4,6)=Fs(4,5);

% Fs56
Fs(5,6)=.5*(1.-((1.-RHOs(1))*A(1)*Fs(1,5)+(1.-RHOs(2))*A(2)*Fs(2,5)+A(3)*Fs(3,5)+A(4)*Fs(4,5))/A(5));

A=A;
Fs=Fs;

%Solve system of equations
POUT = GRAYDIFFSPEC(iclsd,A,EPS,RHOs,HOs,Fs,id,PIN);
%Output
%Convert emissive powers to temperatures
sumq = 0;                              %Check total flux=0?
for i=1:N
    if (id(i)==0)
        T(i)=pow((POUT(i)/sigma),.25);
    else
        q(i)=POUT(i);
    end
    QA(i)=q(i)*A(i);
    sumq=sumq+QA(i);
end
%Output
fprintf('\n surface       T             q [W/m2]         Q [W] \n');
for i = 1:N
    fprintf('%5d %14f %14f %14f \n', i, T(i), q(i), QA(i)/ll);
end
fprintf('Sum of all fluxes = %f \n', sumq);
t=0;
end
%==========================================================================
function [POUT] =  GRAYDIFFSPEC(iclsd,A,EPS,RHOs,HOs,Fs,ID,PIN)
% Function graydifspec
%******************************************************************************
% Routine to fill view factor matrix and solve for missing surface temperatures/fluxes
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
Fs=Fs;
for  i = 2:N
    for  j=1:(i-1)
        Fs(i,j) = A(j)/A(i)*Fs(j,i);
    end
end
Fs=Fs
% If closed configuration, need to also calculate diagonal terms by summation rule
if (iclsd == 1)
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
t=0;
for k=1:N
    for j=1:N
        if (k==j)
            t(k)=Fs(k,j);
        end
    end
end
t=t;
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
B=B;
% Invert POUT-coefficient matrix and multiply by RHS to obtain POUT
POUT = GAUSS(pm,B);         %  can also do: POUT = pm\B', or inv(pm)*B';
end
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
end
%_________________________________________________________________________
function viewFF = view(NO, NARG, ARG, fl)
% Function VIEWFACTORS
%******************************************************************************
% View Factor routine
%
%  *************************************************************************
%  *  THIS SUBROUTINE EVALUATES ALL VIEW FACTORS LISTED IN APPENDIX D      *
%  *  NO = CONFIGURATION NUMBER (1 THROUGH 51),                            *
%  *  NARG = NUMBER OF NON-DIMENSIONAL ARGUMENTS FOR VIEW FACTOR,          *
%  *  ARG = ARRAY CONTAINING THE NARG DIFFERENT ARGUMENTS (in alphabetical *
%  *          order as given in App.D, angles given in degrees),           *
%  *  VIEW = VIEW FACTOR RETURNED BY THE SUBROUTINE                        *
%  *      (FOR NUMBERS 1-9, DF/DL IS GIVEN, WHERE DL IS NONDIMENSIONAL     *
%  *       INFINITESIMAL DIMENSION OF RECEIVING SURFACE).                  *
%  *************************************************************************
A=0; B=0; C=0;
if ((NO < 1) || (NO > 51))
    viewF = warning_i(NO);
end
PI=pi;
epsmin=1e-15;
%PI = 3.1415926E0;
% % %
    if (NO==38)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        A  = ARG(1);
        B  = ARG(2);
        C  = ARG(3);
        if (abs(C)>epsmin)
        X  = A/C;
        Y  = B/C;
        else
            X=0.0; Y=0.0;
        end
        if ((fl==4) || (fl==13) || (fl==15))
        %X=X'
        %Y=Y'
        %A=A'
        %B=B'
        %C=C'
        end
        RTX = sqrt(1.0 + X^2);
        RTY = sqrt(1.0 + Y^2);
        RT  = sqrt(1.0 + X^2 + Y^2);
        if ((abs(X)>epsmin) && (abs(Y)>epsmin))
        viewF = (log(RTX*RTY/RT)+X*RTY*atan(X/RTY)+Y*RTX*atan(Y/RTX)-X*atan(X)-Y*atan(Y))*2.0/(pi*X*Y);
        else
            viewF = 0.0;
        end
        end
    elseif (NO==39)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        H   = ARG(1);
        L   = ARG(2);
        W   = ARG(3);
        HH  = H/L;
        WW  = W/L;
        W2  = WW*WW;
        H2  = HH*HH;
        HW2 = H2+W2;
        HW  = sqrt(H2+W2);
        H12 = H2+1.;
        W12 = W2+1.;
        C1  = W12*H12/(H12+W2);
        C2  = W2*(H12+W2)/W12/HW2;
        C3  = H2*(H12+W2)/H12/HW2;
        viewF = (WW*atan(1./WW)+HH*atan(1./HH)-HW*atan(1./HW) ...
            +.25*(log(C1)+W2*log(C2)+H2*log(C3)))/(PI*WW);
        end
    end
    viewFF=viewF;
end

function [viewf] = warning_i(NO)
fprintf('***ILLEGAL VALUE FOR NO (NO ="\d"), MUST BE 1<=NO>=51 *** \n', NO);
viewf = 0.0;
end

function [viewf] = warning_ii(NO)
fprintf('*** WRONG VALUE OF NARG (NARG =" \d ") FOR NO ="<<NO<<"*** \n', NO);
viewf = 0.0;
end

function [y] = pow(x, n)
y = x.^n;
end

function [d] = KronD(j,k)
if (nargin < 2)
    error('Too few inputs.  See help KronD')
elseif (nargin > 2)
    error('Too many inputs.  See help KronD')
end
if (j == k)
    d = 1;
else
    d = 0;
end
end

%Function parlplates
%***********************************************************************************
function [parlpltf] = parlplates(X1, X2, X3, Y1, Y2, Y3, C)
%  *************************************************************************
%  *  THIS SUBROUTINE EVALUATES THE VIEW FACTOR BETWEEN TWO PARALLEL,      *
%  *  RECTANGULAR pltf OF SIZE X1xY1 AND (X3-X2)x(Y3-Y2), DISPLACED        *
%  *  FROM ANOTHER BY C IN THE Z-DIRECTION, BY X2 IN THE X-DIRECTON, AND   *
%  *  BY Y2 IN THE Y-DIRECTION, IN TERMS OF VIEW FACTORS OF DIRECTLY       *
%  *  OPPOSITE, IDENTICAL RECTANGLES (CONFIG. 38), AS SHOWN IN FIG.4-15b   *
%  *************************************************************************
F=0;
fl=0;
epsmin=1e-15;
NARG = 3;
ARG(1) = X3; ARG(2) = Y3; ARG(3) = C;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = A*view(38,NARG,ARG,fl);
end
%disp('1par');
%F=F' %1
%--------
ARG(2) = Y2;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG,fl);
end
%disp('2par');
%F=F' %2
%--------
ARG(2) = Y3-Y1;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG,fl);
end
%disp('3par');
%F=F' %3
%--------
ARG(2) = Y2-Y1;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F+A*view(38,NARG,ARG,fl);
end
%disp('4par');
%F=F' %4
%--------
ARG(1) = X2;
ARG(2) = Y3;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG,fl);
end
%disp('5par');
%F=F' %5
%--------
ARG(2) = Y2;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F+A*view(38,NARG,ARG,fl);
end
%disp('6par');
%F=F' %6
%--------
ARG(2) = Y3-Y1;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F+A*view(38,NARG,ARG,fl);
end
%disp('7par');
%F=F' %7
%--------
ARG(2) = Y2-Y1;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG,fl);
end
%disp('8par');
%F=F' %8
%--------
ARG(1) = X3-X1;
ARG(2) = Y3;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG,fl);
end
%disp('9par');
%F=F' %9
%--------
ARG(2) = Y2;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F+A*view(38,NARG,ARG,fl);
end
%disp('10par');
%F=F' %10
%--------
ARG(2) = Y3-Y1;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F+A*view(38,NARG,ARG,fl);
end
%disp('11par');
%F=F' %11
%--------
ARG(2) = Y2-Y1;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG,fl);
end
%disp('12par');
%F=F' %12
%--------
ARG(1) = X2-X1;
ARG(2) = Y3;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F+A*view(38,NARG,ARG,fl);
end
%disp('13par');
%F=F' %13
%--------
ARG(2) = Y2;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG,fl);
end
%disp('14par');
%F=F' %14
%--------
ARG(2) = Y3-Y1;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG,fl);
end
%disp('15par');
%F=F' %15
%--------
ARG(2) = Y2-Y1;
A = ARG(1)*ARG(2); fl=fl+1;
if (abs(A)>epsmin)
    F = F+A*view(38,NARG,ARG,fl);
end
%disp('16par');
%F=F' %16
%--------
parlpltf = F/(4.*(X1*Y1));
end

%Function perplates
%***********************************************************************************
function [perppltf] = perpplates(X1, X2, Y1, Y2, Z1, Z2, Z3)

%  *************************************************************************
% *  THIS SUBROUTINE EVALUATES THE VIEW FACTOR BETWEEN TWO PERPENDICULAR,   *
%  *  RECTANGULAR pltf OF SIZE (X2-X1)xZ1 AND (Y2-Y1)x(Z3-Z2), DISPLACED    *
%  *  FROM ANOTHER BY Z2 IN THE Z-DIRECTION, BY X1 IN THE X-DIRECTON, AND   *
%  *  BY Y1 IN THE Y-DIRECTION, IN TERMS OF VIEW FACTORS OF DIRECTLY        *
%  *  ADJACENT, IDENTICAL RECTANGLES (CONFIG. 39), AS SHOWN IN FIG.4-15a    *
%  *************************************************************************
NARG = 3;
F = 0.0;
fl=0;
ARG(1) = Y2;
ARG(2) = Z3;
ARG(3) = X2;
A = X2*Z3; fl=fl+1;
epsmin=1e-15;
if (abs(A*Y2)>epsmin)
    F=A*view(39,NARG,ARG,fl);
end
%disp('1per');
%F=F'
ARG(1) = Y1;
A = X2*Z3; fl=fl+1;
if (abs(A*Y1)>epsmin)
    F = F-A*view(39,NARG,ARG,fl);
end
%disp('2per');
%F=F'
ARG(1) = Y2;
ARG(3) = X1;
A = X1*Z3; fl=fl+1;
if (abs(A*Y2)>epsmin)
    F = F-A*view(39,NARG,ARG,fl);
end
%disp('3per');
%F=F'
ARG(1) = Y1;
A = X1*Z3; fl=fl+1;
if (abs(A*Y1)>epsmin)
    F = F+A*view(39,NARG,ARG,fl);
end
%disp('4per');
%F=F'
ARG(1) = Y2;
ARG(2) = Z2;
A = X1*Z2; fl=fl+1;
if (abs(A*Y2)>epsmin)
    F=F+A*view(39,NARG,ARG,fl);
end
%disp('5per');
%F=F'
ARG(1) = Y1;
A = X1*Z2; fl=fl+1;
if (abs(A*Y1)>epsmin)
    F = F-A*view(39,NARG,ARG,fl);
end
%disp('6per');
%F=F'
ARG(1) = Y2;
ARG(3) = X2;
A = X2*Z2; fl=fl+1;
if (abs(A*Y2)>epsmin)
    F = F-A*view(39,NARG,ARG,fl);
end
%disp('7per');
%F=F'
ARG(1) = Y1;
A = X2*Z2; fl=fl+1;
if (abs(A*Y1)>epsmin)
    F = F+A*view(39,NARG,ARG,fl);
end
%disp('8per');
%F=F'
ARG(1) = Y2;
ARG(2) = (Z3-Z1);
A = X2*(Z3-Z1); fl=fl+1;
if (abs(A*Y2)>epsmin)
    F = F-A*view(39,NARG,ARG,fl);
end
%disp('9per');
%F=F'
ARG(1) = Y1;
A = X2*(Z3-Z1); fl=fl+1;
if (abs(A*Y1)>epsmin)
    F = F+A*view(39,NARG,ARG,fl);
end
%disp('10per');
%F=F'
ARG(1) = Y2;
ARG(3) = X1;
A = X1*(Z3-Z1); fl=fl+1;
if (abs(A*Y2)>epsmin)
    F = F+A*view(39,NARG,ARG,fl);
end
%disp('11per');
%F=F'
ARG(1) = Y1;
A = X1*(Z3-Z1); fl=fl+1;
if (abs(A*Y1)>epsmin)
    F = F-A*view(39,NARG,ARG,fl);
end
%disp('12per');
%F=F'
ARG(1) = Y2;
ARG(2) = (Z2-Z1);
ARG(3) = X2;
A = X2*(Z2-Z1); fl=fl+1;
if (abs(A*Y2)>epsmin)
    F = F+A*view(39,NARG,ARG,fl);
end
%disp('13per');
%F=F'
ARG(1) = Y1;
A = X2*(Z2-Z1); fl=fl+1;
if (abs(A*Y1)>epsmin)
    F = F-A*view(39,NARG,ARG,fl);
end
%disp('14per');
%F=F'
ARG(1) = Y2;
ARG(3) = X1;
A = X1*(Z2-Z1); fl=fl+1;
if (abs(A*Y2)>epsmin)
    F = F-A*view(39,NARG,ARG,fl);
end
%disp('15per');
%F=F'
ARG(1) = Y1;
A = X1*(Z2-Z1); fl=fl+1;
if (abs(A*Y1)>epsmin)
    F = F+A*view(39,NARG,ARG,fl);
end
%disp('16per');
%F=F'
perppltf = F/(2.*(X2-X1)*Z1);
end

function viewF = viewfactors(NO, NARG, ARG, fl)
if (NO==1)
        if (NARG~=1)
            [ viewF ] = warning_ii(NO);
        else
        PHI   = ARG(1);
        viewF = 0.5*cos(PHI*PI/180.);
        end
    elseif (NO==2)
        if (NARG~=2)
            [ viewF ] = warning_ii(NO);
        else
        L  = ARG(1);
        R  = ARG(2);
        RR = R/L;
        viewF = 2.*RR/pow((1.+RR*RR),2);
        end
    elseif (NO==3)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        R   = ARG(1);
        X   = ARG(2);
        XX  = X/R;
        viewF = 2.*XX/pow((1.+XX*XX),2);
        end
    elseif (NO==4)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        R1  = ARG(1);
        R2  = ARG(2);
        Z   = ARG(3);
        RR  = R2/R1;
        if (RR > 1.)
            fprintf('R2 MUST BE LESS THAN R1 \n');
            viewF = 0.;
        else
        ZZ = Z/R1;
        X  = 1.+RR*RR+ZZ*ZZ;
        viewF = 2.*ZZ*(X-2.*RR*RR)*RR/pow((X*X-4.*RR*RR),1.5);
        end
        end
    elseif (NO==5)
        if (NARG~=3)
            [ viewF ] = warning_ii(NO);
        else
        X   = ARG(1);
        Y   = ARG(2);
        PHI = ARG(3);
        YY  = Y/X;
        PHI = PHI*PI/180.;
        viewF = 0.5*YY*pow(sin(PHI),2)/pow((1.+YY*YY-2.*YY*cos(PHI)),1.5);
        end
    elseif (NO==6)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        B   = ARG(1);
        R   = ARG(2);
        PHI = ARG(3);
        BB  = B/R;
        PHI = PHI*PI/180.;
        viewF = atan(BB)*cos(PHI)/PI;
        end
    elseif (NO==7)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        L   = ARG(1);
        R1  = ARG(2);
        R2  = ARG(3);
        RR  = R1/R2;
        BL  = L/R1;
        AA  = BL*BL+RR*RR+1.0;
        viewF = 2.0*RR*BL*BL*AA/pow((AA*AA-4.*RR*RR),1.5);
        end
    elseif (NO==8)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        R1  = ARG(1);
        R2  = ARG(2);
        X   = ARG(3);
        RR  = R1/R2;
        XX  = X/R2;
        AA  = XX*XX-RR*RR+1.;
        BB  = XX*XX+RR*RR+1.;
        viewF = 2.*XX*AA/pow((BB*BB-4.*RR*RR),1.5);
        end
    elseif (NO==9)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        R  = ARG(1);
        X  = ARG(2);
        XX = 0.5*X/R;
        if (XX>0.5E6)
            viewF = 0;
        else
            viewF  = 1.-0.5*XX*(2.*XX*XX+3.)/pow((XX*XX+1.),1.5);
        end
        end
    elseif (NO==10)
        if(NARG~=3)
            [viewF] = warning_ii(NO);
        else
        A=ARG(1);
        B=ARG(2);
        C=ARG(3);
        AA=A/C;
        BB=B/C;
        RTA=sqrt(1.+AA*AA);
        RTB=sqrt(1.+BB*BB);
        viewF = (AA/RTA*atan(BB/RTA)+BB/RTB*atan(AA/RTB))/(2.*PI);
        end
    elseif (NO==11)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        A  = ARG(1);
        B  = ARG(2);
        C  = ARG(3);
        X  = A/B;
        Y  = C/B;
        RT = sqrt(X*X+Y*Y);
        viewF  = (atan(1.0/Y)-Y/RT*atan(1.0/RT))/(2.0*PI);
        end
    elseif (NO==12)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        H  = ARG(1);
        R  = ARG(2);
        HR = H/R;
        viewF  = 1./(1.+HR*HR);
        end
    elseif (NO==13)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        A  = ARG(1);
        H  = ARG(2);
        R  = ARG(3);
        HH = H/A;
        RR = R/A;
        ZZ = 1. + HH*HH+RR*RR;
        viewF = .5*(1.-(ZZ-2.*RR*RR)/sqrt(ZZ*ZZ-4.*RR*RR));
        end
    elseif (NO==14)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        H  = ARG(1);
        L  = ARG(2);
        R  = ARG(3);
        HH = H/L;
        RR = R/L;
        ZZ = 1.0+HH*HH+RR*RR;
        viewF = 0.5*HH*(ZZ/sqrt(ZZ*ZZ-4.*RR*RR)-1.);
        end
    elseif (NO==15)
        if (NARG~=3)
            [ viewF ] = warning_ii(NO);
        else
        H  = ARG(1);
        L  = ARG(2);
        R  = ARG(3);
        HH = H/R;
        LL = L/R;
        X  = pow((1.+HH),2)+LL*LL;
        Y  = pow((1.-HH),2)+LL*LL;
        viewF = LL/HH*(1./LL*atan(LL/sqrt(HH*HH-1.)) ...
            +(X-2.*HH)/sqrt(X*Y)*atan(sqrt(X/Y*(HH-1.)/(HH+1.))) ...
            -atan(sqrt((HH-1.)/(HH+1.))))/PI;
        end
    elseif (NO==16)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        H  = ARG(1);
        R  = ARG(2);
        viewF = pow((R/H),2);
        end
    elseif (NO==17)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        H  = ARG(1);
        R  = ARG(2);
        HH = H/R;
        viewF = (atan(sqrt(1./(HH*HH-1.)))-sqrt(HH*HH-1.)/HH/HH)/PI;
        end
    elseif (NO==18)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        H   = ARG(1);
        R   =ARG(2);
        PHI =ARG(3);
        PHI = PHI*PI/180.;
        PHIMAX = acos(R/H);
        if (PHI>PHIMAX)
            viewF = 0.; 
            fprintf(' PHI MUST BE LESS THAN ACOS(R/H),"<<180.*PHIMAX/PI\n ');
        else
        viewF = pow((R/H),2)*cos(PHI);
        end
        end
    elseif (NO==19)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        H = ARG(1);
        L = ARG(2);
        R = ARG(3);
        HH = H/R;
        LL = L/R;
        if (HH<-1.0)
            viewF = 0.0;
        elseif (HH<+1.0)
            X = LL*LL+HH*HH;
            Y = sqrt(X-1.);
            viewF = (HH/pow(X,1.5)*acos(-HH/LL/Y)-Y*sqrt(1.-HH*HH)/X ...
                -asin(Y/LL/LL)+PI/2.)/PI;
        else
            X = LL*LL+HH*HH;
            viewF = pow(H/X,1.5);
        end
        end
    elseif (NO==20)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        H = ARG(1);
        R = ARG(2);
        Z = ARG(3);
        if ((Z<0.0) || (Z>H)) 
            fprintf(' INVALID Z"\n ');
            viewF = 0;
        else
        HH = 0.5*H/R;
        ZZ = 0.5*Z/R;
        viewF = 1.+HH-(ZZ*ZZ+0.5)/sqrt(ZZ*ZZ+1.)-(pow((HH-ZZ),2)+0.5)/sqrt(pow((HH-ZZ),2)+1.);
        end
        end
    elseif (NO==21)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        R  = ARG(1);
        Z  = ARG(2);
        if (Z<0.0)
            fprintf(' INVALID Z"\n ');
            viewF = 0.;
        else
        ZZ = Z/R;
        if (ZZ>1.E4)
            viewF = pow(ZZ,-3.);
        else
            viewF = 0.5*((ZZ*ZZ+2.)/sqrt(ZZ*ZZ+4.)-ZZ);
        end
        end
        end
    elseif (NO==22)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        R1 = ARG(1);
        R2 = ARG(2);
        Z  = ARG(3);
        RR = R2/R1;
        if (RR>1.0)
            fprintf(' R2 MUST BE LESS THAN R1"\n ');
            viewF = 0.0;
        else
        ZZ = Z/R1;
        X  = 1.+RR*RR+ZZ*ZZ;
        viewF = 0.5*ZZ*(X/sqrt(X*X-4.*RR*RR)-1.);
        end
        end
    elseif (NO==23)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        L  = ARG(1);
        X  = ARG(2);
        XX = X/L;
        PHI = ARG(3);
        PHI = PHI*PI/180.;
        viewF = 0.5*(1.+(cos(PHI)-XX)/sqrt(1.+XX*XX-2.*XX*cos(PHI)));
        end
    elseif (NO==24)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        PHI1 = ARG(1);
        PHI2 = ARG(2);
        if (PHI2<PHI1)
            fprintf(' WARNING, PHI2 MUST BE GREATER THAN PHI1"\n ');
            PHI  = PHI1;
            PHI1 = PHI2;
            PHI2 = PHI;
        end
        if (PHI1>90.)
            PHI1 = 90.;
            fprintf(' PHI1 MUST BE LESS THAN 90 OR GREATER THAN -90"\n ');
        elseif (PHI1<-90.)
            PHI1 = -90.;
            fprintf(' PHI1 MUST BE LESS THAN 90 OR GREATER THAN -90"\n ');
        end
        PHI1 = PHI1*PI/180.;
        PHI2 = PHI2*PI/180.;
        viewF = 0.5*(sin(PHI2)-sin(PHI1));
        end
    elseif (NO==25)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        A = ARG(1);
        B = ARG(2);
        R = ARG(3);
        AA = A/R;
        BB = B/R;
        viewF = AA/(AA*AA+BB*BB);
        end
    elseif (NO==26)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        A  = ARG(1);
        B  = ARG(2);
        C  = ARG(3);
        X  = A/C;
        Y  = B/C;
        XX = sqrt(1.+X*X);
        YY = sqrt(1.+Y*Y);
        viewF = (YY*atan(X/YY)-atan(X)+X*Y/XX*atan(Y/XX))/(PI*Y);
        end
    elseif (NO==27)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        A  = ARG(1);
        B  = ARG(2);
        C  = ARG(3);
        X  = A/B;
        Y  = C/B;
        if (Y>1.E6)
            viewF = 0.;
        else
        XY = X*X+Y*Y;
        YY = 1.+Y*Y;
        viewF = (atan(1/Y)+.5*Y*log(Y*Y*(XY+1.)/YY/XY)-Y/sqrt(XY)*atan(1./sqrt(XY)))/PI;
        end
        end
    elseif (NO==28)
        if (NARG~=4)
            [viewF] = warning_ii(NO);
        else        
        H  = ARG(1);
        R  = ARG(2);
        S  = ARG(3);
        X  = ARG(4);
        HH = H/R;
        SS = S/R;
        XX = X/R;
        C  = SS*SS+XX*XX;
        CC = sqrt(C);
        A  = HH*HH+C-1.;
        B  = HH*HH-C+1.;
        viewF = SS/C*(1.-((1./PI)*(acos(B/A)-0.5*(sqrt(A*A+4.*HH*HH)*acos(B/A/CC)-B*asin(1./CC))/HH)) ...
            -0.25*A/HH);
        end
    elseif (NO==29)
        if (NARG~=1)
            [viewF] = warning_ii(NO);
        else
        PHI  = ARG(1);
        viewF = 0.5*(1.+cos(PHI*PI/180.));
        end
    elseif (NO==30)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        A   = ARG(1);
        R1  = ARG(2);
        R2  = ARG(3);
        RR1 = R1/A;
        RR2 = R2/A;
        if (RR2>1.0 )
            fprintf(' R@ MUST BE LESS THAN A"\n ');
            viewF = 0.;
        else
        viewF = RR2*RR2/pow((1.+RR1*RR1),1.5);
        end
        end
    elseif (NO==31)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        % IDENTICAL TO #21
        R  = ARG(1);
        X  = ARG(2);
        XX = 0.5*X/R;
        if (XX > 65000.0)
            viewF = 0.;
        else
        viewF = (XX*XX+0.5)/sqrt(1.+XX*XX)-XX;
        end
        end
    elseif (NO==32)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        H  = ARG(1);
        W  = ARG(2);
        HW = H/W;
        viewF = sqrt(1.+HW*HW)-HW;
        end
    elseif (NO==33)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        H  = ARG(1);
        W  = ARG(2);
        HW = H/W;
        viewF = 0.5*(1.+HW-sqrt(1.+HW*HW));
        end
    elseif (NO==34)
        if (NARG~=1)
            [viewF] = warning_ii(NO);
        else
        if (ARG(1)>180.)
            viewF = 0.;
        else
        ALPHA2 =ARG(1)*PI/360.;
        viewF = 1.-sin(ALPHA2);
        end
        end
    elseif (NO==35)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        R  = ARG(1);
        S  = ARG(2);
        X  = 1+0.5*S/R;
        viewF = (asin(1./X)+sqrt(X*X-1.)-X)/PI;
        end
    elseif (NO==36)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        R1 = ARG(1);
        R2 = ARG(2);
        S  = ARG(3);
        RR = R2/R1;
        SS = S/R1;
        C  = 1.+RR+SS;
        viewF = (PI+sqrt(C*C-pow((RR+1.),2))-sqrt(C*C-pow((RR-1.),2)) ...
            +(RR-1)*acos((RR-1.)/C)-(RR+1)*acos((RR+1.)/C))/(2.*PI);
        end
    elseif (NO==37)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        A   = ARG(1);
        B1  = ARG(2);
        B2  = ARG(3);
        BB1 = B1/A;
        BB2 = B2/A;
        viewF = (atan(BB1)-atan(BB2))/(2.*PI);
        end
    elseif (NO==40)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        A   = ARG(1);
        R1  = ARG(2);
        R2  = ARG(3);
        RR1 = R1/A;
        RR2 = R2/A;
        X   = 1.+(1.+RR2*RR2)/(RR1*RR1);
        viewF = 0.5*(X-sqrt(X*X-4.*pow((R2/R1),2)));
        end
    elseif (NO==41)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        L  = ARG(1);
        R1 = ARG(2);
        R2 = ARG(3);
        RR = R1/R2;
        if (RR>1e0)
            fprintf(' R2 MUST BE LESS THAN R1"\n ');
            viewF = 0.;
        else
        BL = L/R2;
        if (BL>1000.0)
            viewF = 0.;
        else
        AA = BL*BL+RR*RR-1.;
        BB = BL*BL-RR*RR+1.;
        viewF = BB/(8.*RR*BL)+(acos(AA/BB) ...
            -sqrt(pow(((AA+2.)/RR),2)-4.)/(2.*BL)*acos(AA*RR/BB) ...
            -AA/(2.*RR*BL)*asin(RR))/(2.*PI);
        end
        end
        end
    elseif (NO==42)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        H  = ARG(1);
        R  = ARG(2);
        HH = 0.5*H/R;
        viewF = 1.+HH-sqrt(1.+HH*HH);
        end
    elseif (NO==43)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        H  = ARG(1);
        R  = ARG(2);
        HH = 0.5*H/R;
        viewF = 2.*HH*(sqrt(1.+HH*HH)-HH);
        end
    elseif (NO==44)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        H  = ARG(1);
        R1 = ARG(2);
        R2 = ARG(3);
        RR = R2/R1;
        if (RR<1.0)
            fprintf(' R2 MUST BE GREATER THAN R1"\n ');
            viewF = 0.;
        else
        HH = H/R1;
        AA = HH*HH+4.*RR*RR;
        BB = 4.0*(RR*RR-1.0);
        viewF = 1.-1./RR-(sqrt(AA)-HH)/(4.*RR)+(2./RR*atan(sqrt(BB)/HH) ...
            -HH/(2.*RR)*(sqrt(AA)/HH*asin((HH*HH*(1.-2./RR/RR)+BB)/(AA-4.)) ...
            -asin((RR*RR-2.)/RR/RR)))/PI;
        end
        end
    elseif (NO==45)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        H  = ARG(1);
        R1 = ARG(2);
        R2 = ARG(3);
        RR = R2/R1;
        if (RR<1e0)
            fprintf(' R2 MUST BE GREATER THAN R1"\n ');
            viewF = 0.;
        else
        HH = H/R1;
        if (HH>6250.)
            viewF = 1.;
        else
        AA = HH*HH+RR*RR-1.0;
        BB = HH*HH-RR*RR+1.0;
        CC = HH*HH+RR*RR+1.0;
        viewF = 1.0-AA/(4.0*HH)-(acos(BB/AA) ...
            -sqrt(CC*CC-4.*RR*RR)/(2.*HH)*acos(BB/RR/AA) ...
            -BB/(2.0*HH)*asin(1.0/RR))/PI;
        end
        end
        end
    elseif (NO==46)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        H  = ARG(1);
        R1 = ARG(2);
        R2 = ARG(3);
        RR = R1/R2;
        if (RR>1.0)
            fprintf(' R2 MUST BE GREATER THAN R1"\n ');
            viewF = 0.;
        else
        HH = H/R1;
        AA = 1.0-RR*RR-HH*HH;
        BB = 1.0-RR*RR+HH*HH;
        CC = 1.0+HH*HH+RR*RR;
        X  = sqrt(1.0-RR*RR);
        Y  = RR*AA/BB;
        viewF = (RR*(atan(X/HH)-atan(2.*X/HH))+0.25*HH*(asin(2.*RR*RR-1.)-asin(RR)) ...
            +0.25*X*X/HH*(0.5*PI+asin(RR))-0.25*sqrt(CC*CC-4.*RR*RR)/HH*(0.5*PI+asin(Y))  ...
            +0.25*sqrt(4.+HH*HH)*(0.5*PI+asin(1.-2.*RR*RR*HH*HH/(4.*X*X+HH*HH))))/PI;
        end
        end
    elseif (NO==47)
        if (NARG~=3)
            [viewF] = warning_ii(NO);
        else
        D  = ARG(1);
        L1 = ARG(2);
        L2 = ARG(3);
        D1 = D/L1;
        D2 = D/L2;
        viewF = 0.25*atan(1./sqrt(D1*D1+D2*D2+D1*D1*D2*D2))/PI;
        end
    elseif (NO==48)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        A=ARG(1);
        R=ARG(2);
        RR=R/A;
        viewF = 0.5*(1.-1./sqrt(1.+RR*RR));
        end
    elseif (NO==49)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        A=ARG(1);
        R=ARG(2);
        RR=R/A;
        viewF = 1./sqrt(1.+RR*RR);
        end
    elseif (NO==50)
        if (NARG~=4)
            [viewF] = warning_ii(NO);
        else
        R1 = ARG(1);
        R2 = ARG(2);
        S  = ARG(3);
        OM = ARG(4);
        SS = S/R1;
        RR = R2/R1;
        A  = asin(1./(SS+1.))*180./PI;
        if (OM<A)
            fprintf(' Formula only valid for OM > "<<A<<"deg"\n ');
        else
        OM  = OM*PI/180.;
        SR  = 1.+SS+RR/tan(OM);
        viewF = 0.5*(1.0-SR/sqrt(SR*SR+RR*RR));
        end
        end
    elseif (NO==51)
        if (NARG~=2)
            [viewF] = warning_ii(NO);
        else
        D  = ARG(1);
        S  = ARG(2);
        SS = D/S;
        if (SS>1.0)
            fprintf(' D MUST BE LESS THAN S"\n ');
            viewF = 0.;
        else
        viewF = SS*acos(SS)+1.-sqrt(1.-SS*SS);   
        end
        end
end
end