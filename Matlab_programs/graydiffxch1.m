function graydiffxch1
% Program Graydiffxch

N = 6;
sigma = 5.670E-8;

%dimensions
ww = 1E1;
hh = 0.3;
ll = 1E1;

%surfaces 1&3(bottom and top)
A(1) = ww*ll;	% ww x ll
A(3) = ww*ll;
HO(1) = 0.;	%no external irradiation
HO(3) = 0.;
EPS(1) = 0.3;
EPS(3) = 0.3;
id(1) = 1;	%T specified
id(3) = 1;
T(1) = 400.;	%T in K
T(3) = 400.;

%surfaces 2&4(left and right)
A(2) = hh*ll;	% ww x ll
A(4) = hh*ll;
HO(2) = 0.;	%no external irradiation
HO(4) = 0.;
EPS(2) = 0.8;
EPS(4) = 0.8;
id(2) = 1;	%T specified
id(4) = 1;
T(2) = 405.;	%T in K
T(4) = 405.;

%surfaces 5&6(front and back)
A(5) = hh*ww;	% ww x ll
A(6) = hh*ww;
HO(5) = 0.0;	%no external irradiation
HO(6) = 0.0;
EPS(5) = 0.8;
EPS(6) = 0.8;
id(5) = 1;	%T specified
id(6) = 1;
T(5) = 405.;	%T in K
T(6) = 405.;

%fill PIN array with q ant T
for  i = 1:N
    if id(i) == 0
        PIN(i) = q(i);
    else
        PIN(i) = sigma*T(i)^4;
    end
end
% view factors; since configuration is closed (iclsd = 1), diagonal terms are not needed
iclsd = 1;

%F12
arg(1) = hh;
arg(2) = ll;
arg(3) = ww;	%h,l,w
F(1,2) = view(39,3,arg);

%F13
arg(1) = ww;
arg(2) = ll;
arg(3) = hh;	%a,b,c
F(1,3) = view(38,3,arg);
F(1,4) = F(1,2);

%F15
arg(1) = hh;
arg(2) = ww;
arg(3) = ll;	%h,l,w
F(1,5) = view(39,3,arg);
F(1,6) = F(1,5);

%F23
arg(1) = ww;
arg(2) = ll;
arg(3) = hh;	%h,l,w
F(2,3) = view(39,3,arg);

%F24
arg(1) = ll;
arg(2) = hh;
arg(3) = ww;	%a,b,c
F(2,4) = view(38,3,arg);

%F25
arg(1) = ww;
arg(2) = hh;
arg(3) = ll;	%h,l,w
F(2,5) = view(39,3,arg);
F(2,6) = F(2,5);
F(3,4) = F(1,4);
F(3,5) = F(1,5);
F(3,6) = F(1,6);
F(4,5) = F(2,5);
F(4,6) = F(2,6);

%F56
arg(1) = ww;
arg(2) = hh;
arg(3) = ll;	%a,b,c
F(5,6) = view(38,3,arg);

%Solve system of equations
POUT  =  GRAYDIFF(iclsd, A, EPS, HO, F, id, PIN);
%Output
%Convert emissive powers to temperatures
sumq = 0;                              %Check total flux = 0?
for i = 1:N
    if  id(i) == 0
        %T(i) = (POUT(i)/sigma)^(1/4);
    else
        q(i) = POUT(i);
    end
    QA(i) = q(i)*A(i);
    sumq = sumq + QA(i);
end

% Output
fprintf('\n surface       T             q [W/m2]         Q [W] \n');
for i = 1:N
    fprintf('%5d %14f %14f %14f \n', i, T(i), q(i), QA(i)/ll);
end
fprintf('Sum of all fluxes  =  %f \n', sumq);
	
	
    
%==========================================================================
%
%==========================================================================
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

%for k=1:N
    %s=0.0;
    %for j=1:N
        %s=s+F(k,j);
    %end
%end

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

%__________________________________________________________________________
function viewF = view(NO, NARG, ARG)
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


if (NO < 1 || NO > 51)
    viewF = warning_i(NO);
end

PI = 3.1415926E0;
% % %
switch (NO)
    case 1
        if(NARG~=1)
            [viewF] = warning_ii(NO);
            return
        end
        PHI   = ARG(1);
        viewF = 0.5*cos(PHI*PI/180.);
    case 2
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return
        end
        L  = ARG(1);
        R  = ARG(2);
        RR = R/L;
        viewF = 2.*RR/pow((1.+RR*RR),2);
    case 3
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return
        end
        R   = ARG(1);
        X   = ARG(2);
        XX  = X/R;
        viewF = 2.*XX/pow((1.+XX*XX),2);
    case 4
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        R1  = ARG(1);
        R2  = ARG(2);
        Z   = ARG(3);
        RR  = R2/R1;
        if (RR > 1.)
            fprintf('R2 MUST BE LESS THAN R1 \n');
            viewF = 0.;
            return;
        end
        ZZ = Z/R1;
        X  = 1.+RR*RR+ZZ*ZZ;
        viewF = 2.*ZZ*(X-2.*RR*RR)*RR/pow((X*X-4.*RR*RR),1.5);
    case 5
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        X   = ARG(1);
        Y   = ARG(2);
        PHI = ARG(3);
        YY  = Y/X;
        PHI = PHI*PI/180.;
        viewF = 0.5*YY*pow(sin(PHI),2)/pow((1.+YY*YY-2.*YY*cos(PHI)),1.5);
    case 6
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        B   = ARG(1);
        R   = ARG(2);
        PHI = ARG(3);
        BB  = B/R;
        PHI = PHI*PI/180.;
        viewF = atan(BB)*cos(PHI)/PI;
    case 7
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        L   = ARG(1);
        R1  = ARG(2);
        R2  = ARG(3);
        RR  = R1/R2;
        BL  = L/R1;
        AA  = BL*BL+RR*RR+1.;
        viewF = 2.*RR*BL*BL*AA/pow((AA*AA-4.*RR*RR),1.5);
    case 8
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        R1  = ARG(1);
        R2  = ARG(2);
        X   = ARG(3);
        RR  = R1/R2;
        XX  = X/R2;
        AA  = XX*XX-RR*RR+1.;
        BB  = XX*XX+RR*RR+1.;
        viewF = 2.*XX*AA/pow((BB*BB-4.*RR*RR),1.5);
    case 9
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        R  = ARG(1);
        X  = ARG(2);
        XX = 0.5*X/R;
        if(XX>.5E6)
            viewF = 0;
        else
            viewF  = 1.-0.5*XX*(2.*XX*XX+3.)/pow((XX*XX+1.),1.5);
        end
    case 10
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        A=ARG(1);
        B=ARG(2);
        C=ARG(3);
        AA=A/C;
        BB=B/C;
        RTA=sqrt(1.+AA*AA);
        RTB=sqrt(1.+BB*BB);
        viewF = (AA/RTA*atan(BB/RTA)+BB/RTB*atan(AA/RTB))/(2.*PI);
    case 11
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        A  = ARG(1);
        B  = ARG(2);
        C  = ARG(3);
        X  = A/B;
        Y  = C/B;
        RT = sqrt(X*X+Y*Y);
        viewF  = (atan(1./Y)-Y/RT*atan(1./RT))/(2.*PI);
    case 12
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        H  = ARG(1);
        R  = ARG(2);
        HR = H/R;
        viewF  = 1./(1.+HR*HR);
    case 13
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        A  = ARG(1);
        H  = ARG(2);
        R  = ARG(3);
        HH = H/A;
        RR = R/A;
        ZZ = 1. + HH*HH+RR*RR;
        viewF = .5*(1.-(ZZ-2.*RR*RR)/sqrt(ZZ*ZZ-4.*RR*RR));
    case 14
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        H  = ARG(1);
        L  = ARG(2);
        R  = ARG(3);
        HH = H/L;
        RR = R/L;
        ZZ = 1.+HH*HH+RR*RR;
        viewF = .5*HH*(ZZ/sqrt(ZZ*ZZ-4.*RR*RR)-1.);
    case 15
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
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
    case 16
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        
        H  = ARG(1);
        R  = ARG(2);
        viewF = pow((R/H),2);
    case 17
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        H  = ARG(1);
        R  = ARG(2);
        HH = H/R;
        viewF = (atan(sqrt(1./(HH*HH-1.)))-sqrt(HH*HH-1.)/HH/HH)/PI;
    case 18
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        H   = ARG(1);
        R   =ARG(2);
        PHI =ARG(3);
        PHI = PHI*PI/180.;
        PHIMAX = acos(R/H);
        if(PHI>PHIMAX)
            viewF = 0.;
            
            fprintf(' PHI MUST BE LESS THAN ACOS(R/H),"<<180.*PHIMAX/PI\n ');
            return
        end
        
        viewF = pow((R/H),2)*cos(PHI);
    case 19
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        H = ARG(1);
        L = ARG(2);
        R = ARG(3);
        HH = H/R;
        LL = L/R;
        if(HH<-1.)
            viewF = 0.;
        elseif(HH<+1.)
            X = LL*LL+HH*HH;
            Y = sqrt(X-1.);
            viewF = (HH/pow(X,1.5)*acos(-HH/LL/Y)-Y*sqrt(1.-HH*HH)/X ...
                -asin(Y/LL/LL)+PI/2.)/PI;
        else
            X = LL*LL+HH*HH;
            viewF = pow(H/X,1.5);
        end
    case 20
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        H = ARG(1);
        R = ARG(2);
        Z = ARG(3);
        if(Z<0. || Z>H)
            
            fprintf(' INVALID Z"\n ');
            viewF = 0;
            return
        end
        HH = .5*H/R;
        ZZ = .5*Z/R;
        viewF = 1.+HH-(ZZ*ZZ+0.5)/sqrt(ZZ*ZZ+1.)-(pow((HH-ZZ),2)+0.5)/sqrt(pow((HH-ZZ),2)+1.);
        
    case 21
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        R  = ARG(1);
        Z  = ARG(2);
        if(Z<0.)
            fprintf(' INVALID Z"\n ');
            viewF = 0.;
            return
        end
        ZZ = Z/R;
        if(ZZ>1.E4)
            viewF = pow(ZZ,-3.);
        else
            viewF = 0.5*((ZZ*ZZ+2.)/sqrt(ZZ*ZZ+4.)-ZZ);
        end
    case 22
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        R1 = ARG(1);
        R2 = ARG(2);
        Z  = ARG(3);
        RR = R2/R1;
        if(RR>1.)
            fprintf(' R2 MUST BE LESS THAN R1"\n ');
            viewF = 0.;
            return
        end
        ZZ = Z/R1;
        X  = 1.+RR*RR+ZZ*ZZ;
        viewF = 0.5*ZZ*(X/sqrt(X*X-4.*RR*RR)-1.);
    case 23
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        L  = ARG(1);
        X  = ARG(2);
        XX = X/L;
        PHI = ARG(3);
        PHI = PHI*PI/180.;
        viewF = 0.5*(1.+(cos(PHI)-XX)/sqrt(1.+XX*XX-2.*XX*cos(PHI)));
    case 24
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        PHI1 = ARG(1);
        PHI2 = ARG(2);
        if(PHI2<PHI1)
            fprintf(' WARNING, PHI2 MUST BE GREATER THAN PHI1"\n ');
            PHI  = PHI1;
            PHI1 = PHI2;
            PHI2 = PHI;
        end
        if(PHI1>90.)
            PHI1 = 90.;
            fprintf(' PHI1 MUST BE LESS THAN 90 OR GREATER THAN -90"\n ');
        elseif(PHI1<-90.)
            PHI1 = -90.;
            fprintf(' PHI1 MUST BE LESS THAN 90 OR GREATER THAN -90"\n ');
        end
        PHI1 = PHI1*PI/180.;
        PHI2 = PHI2*PI/180.;
        viewF = 0.5*(sin(PHI2)-sin(PHI1));
    case 25
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        A = ARG(1);
        B = ARG(2);
        R = ARG(3);
        AA = A/R;
        BB = B/R;
        viewF = AA/(AA*AA+BB*BB);
        
    case 26
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        A  = ARG(1);
        B  = ARG(2);
        C  = ARG(3);
        X  = A/C;
        Y  = B/C;
        XX = sqrt(1.+X*X);
        YY = sqrt(1.+Y*Y);
        viewF = (YY*atan(X/YY)-atan(X)+X*Y/XX*atan(Y/XX))/(PI*Y);
    case 27
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        A  = ARG(1);
        B  = ARG(2);
        C  = ARG(3);
        X  = A/B;
        Y  = C/B;
        if (Y>1.E6)
            viewF = 0.;
            return;
        end
        XY = X*X+Y*Y;
        YY = 1.+Y*Y;
        viewF = (atan(1/Y)+.5*Y*log(Y*Y*(XY+1.)/YY/XY)-Y/sqrt(XY)*atan(1./sqrt(XY)))/PI;
        
    case 28
        if(NARG~=4)
            [viewF] = warning_ii(NO);
            return;
        end
        
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
    case 29
        if(NARG~=1)
            [viewF] = warning_ii(NO);
            return;
        end
        PHI  = ARG(1);
        viewF = 0.5*(1.+cos(PHI*PI/180.));
    case 30
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        A   = ARG(1);
        R1  = ARG(2);
        R2  = ARG(3);
        RR1 = R1/A;
        RR2 = R2/A;
        if(RR2>1.)
            fprintf(' R@ MUST BE LESS THAN A"\n ');
            viewF = 0.;
            return
        end
        viewF = RR2*RR2/pow((1.+RR1*RR1),1.5);
    case 31
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        % IDENTICAL TO #21
        R  = ARG(1);
        X  = ARG(2);
        XX = 0.5*X/R;
        if(XX > 65000.0)
            viewF = 0.;
            return
        end
        viewF = (XX*XX+0.5)/sqrt(1.+XX*XX)-XX;
        
    case 32
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        H  = ARG(1);
        W  = ARG(2);
        HW = H/W;
        viewF = sqrt(1.+HW*HW)-HW;
    case 33
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        H  = ARG(1);
        W  = ARG(2);
        HW = H/W;
        viewF = 0.5*(1.+HW-sqrt(1.+HW*HW));
        
    case 34
        if(NARG~=1)
            [viewF] = warning_ii(NO);
            return;
        end
        if(ARG(1)>180.)
            viewF = 0.;
            return
        end
        ALPHA2 =ARG(1)*PI/360.;
        viewF = 1.-sin(ALPHA2);
    case 35
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        R  = ARG(1);
        S  = ARG(2);
        X  = 1+0.5*S/R;
        viewF = (asin(1./X)+sqrt(X*X-1.)-X)/PI;
    case 36
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        R1 = ARG(1);
        R2 = ARG(2);
        S  = ARG(3);
        RR = R2/R1;
        SS = S/R1;
        C  = 1.+RR+SS;
        viewF = (PI+sqrt(C*C-pow((RR+1.),2))-sqrt(C*C-pow((RR-1.),2)) ...
            +(RR-1)*acos((RR-1.)/C)-(RR+1)*acos((RR+1.)/C))/(2.*PI);
    case 37
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        A   = ARG(1);
        B1  = ARG(2);
        B2  = ARG(3);
        BB1 = B1/A;
        BB2 = B2/A;
        viewF = (atan(BB1)-atan(BB2))/(2.*PI);
    case 38
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        A  = ARG(1);
        B  = ARG(2);
        C  = ARG(3);
        X  = A/C;
        Y  = B/C;
        if(X<0.002 || Y<0.002)
            viewF = 0.0;
            return
        end
        
        RTX = sqrt(1.+X*X);
        RTY = sqrt(1.+Y*Y);
        RT  = sqrt(1.+X*X+Y*Y);
        viewF = (log(RTX*RTY/RT)+X*RTY*atan(X/RTY)+Y*RTX*atan(Y/RTX) ...
            -X*atan(X)-Y*atan(Y))*2./(PI*X*Y);
        
    case 39
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
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
        
    case 40
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        A   = ARG(1);
        R1  = ARG(2);
        R2  = ARG(3);
        RR1 = R1/A;
        RR2 = R2/A;
        X   = 1.+(1.+RR2*RR2)/(RR1*RR1);
        viewF = 0.5*(X-sqrt(X*X-4.*pow((R2/R1),2)));
        
    case 41
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        L  = ARG(1);
        R1 = ARG(2);
        R2 = ARG(3);
        RR = R1/R2;
        if(RR>1.)
            fprintf(' R2 MUST BE LESS THAN R1"\n ');
            viewF = 0.;
            return
        end
        BL = L/R2;
        if(BL>1000.)
            viewF = 0.;
            return
        end
        AA = BL*BL+RR*RR-1.;
        BB = BL*BL-RR*RR+1.;
        viewF = BB/(8.*RR*BL)+(acos(AA/BB) ...
            -sqrt(pow(((AA+2.)/RR),2)-4.)/(2.*BL)*acos(AA*RR/BB) ...
            -AA/(2.*RR*BL)*asin(RR))/(2.*PI);
    case 42
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        H  = ARG(1);
        R  = ARG(2);
        HH = 0.5*H/R;
        viewF = 1.+HH-sqrt(1.+HH*HH);
    case 43
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        H  = ARG(1);
        R  = ARG(2);
        HH = 0.5*H/R;
        viewF = 2.*HH*(sqrt(1.+HH*HH)-HH);
    case 44
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        H  = ARG(1);
        R1 = ARG(2);
        R2 = ARG(3);
        RR = R2/R1;
        if(RR<1.)
            fprintf(' R2 MUST BE GREATER THAN R1"\n ');
            viewF = 0.;
            return
        end
        HH = H/R1;
        AA = HH*HH+4.*RR*RR;
        BB = 4.*(RR*RR-1.);
        viewF = 1.-1./RR-(sqrt(AA)-HH)/(4.*RR)+(2./RR*atan(sqrt(BB)/HH) ...
            -HH/(2.*RR)*(sqrt(AA)/HH*asin((HH*HH*(1.-2./RR/RR)+BB)/(AA-4.)) ...
            -asin((RR*RR-2.)/RR/RR)))/PI;
    case 45
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        H  = ARG(1);
        R1 = ARG(2);
        R2 = ARG(3);
        RR = R2/R1;
        if(RR<1.)
            fprintf(' R2 MUST BE GREATER THAN R1"\n ');
            viewF = 0.;
            return
        end
        HH = H/R1;
        if(HH>6250.)
            viewF = 1.;
            return
        end
        AA = HH*HH+RR*RR-1.;
        BB = HH*HH-RR*RR+1.;
        CC = HH*HH+RR*RR+1.;
        viewF = 1.-AA/(4.*HH)-(acos(BB/AA) ...
            -sqrt(CC*CC-4.*RR*RR)/(2.*HH)*acos(BB/RR/AA) ...
            -BB/(2.*HH)*asin(1./RR))/PI;
    case 46
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        H  = ARG(1);
        R1 = ARG(2);
        R2 = ARG(3);
        RR = R1/R2;
        if(RR>1.)
            fprintf(' R2 MUST BE GREATER THAN R1"\n ');
            viewF = 0.;
            return
        end
        HH = H/R1;
        AA = 1.-RR*RR-HH*HH;
        BB = 1.-RR*RR+HH*HH;
        CC = 1.+HH*HH+RR*RR;
        X  = sqrt(1.-RR*RR);
        Y  = RR*AA/BB;
        viewF = (RR*(atan(X/HH)-atan(2.*X/HH))+0.25*HH*(asin(2.*RR*RR-1.)-asin(RR)) ...
            +0.25*X*X/HH*(0.5*PI+asin(RR))-0.25*sqrt(CC*CC-4.*RR*RR)/HH*(0.5*PI+asin(Y))  ...
            +0.25*sqrt(4.+HH*HH)*(0.5*PI+asin(1.-2.*RR*RR*HH*HH/(4.*X*X+HH*HH))))/PI;
        
    case 47
        if(NARG~=3)
            [viewF] = warning_ii(NO);
            return;
        end
        D  = ARG(1);
        L1 = ARG(2);
        L2 = ARG(3);
        D1 = D/L1;
        D2 = D/L2;
        viewF = 0.25*atan(1./sqrt(D1*D1+D2*D2+D1*D1*D2*D2))/PI;
        
    case 48
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        A=ARG(1);
        R=ARG(2);
        RR=R/A;
        viewF = 0.5*(1.-1./sqrt(1.+RR*RR));
        
    case 49
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        A=ARG(1);
        R=ARG(2);
        RR=R/A;
        viewF = 1./sqrt(1.+RR*RR);
        
    case 50
        if(NARG~=4)
            [viewF] = warning_ii(NO);
            return;
        end
        R1 = ARG(1);
        R2 = ARG(2);
        S  = ARG(3);
        OM = ARG(4);
        SS = S/R1;
        RR = R2/R1;
        A  = asin(1./(SS+1.))*180./PI;
        if(OM<A)
            fprintf(' Formula only valid for OM > "<<A<<"deg"\n ');
        end
        OM  = OM*PI/180.;
        SR  = 1.+SS+RR/tan(OM);
        viewF = 0.5*(1.-SR/sqrt(SR*SR+RR*RR));
    case 51
        if(NARG~=2)
            [viewF] = warning_ii(NO);
            return;
        end
        D  = ARG(1);
        S  = ARG(2);
        SS = D/S;
        if(SS>1)
            fprintf(' D MUST BE LESS THAN S"\n ');
            viewF = 0.;
            return
        end
        viewF = SS*acos(SS)+1.-sqrt(1.-SS*SS);
        
end

%__________________________________________________________________________

function [viewf] = warning_i(NO)
fprintf('***ILLEGAL VALUE FOR NO (NO = %d), MUST BE 1<=NO>=51 *** \n', NO);
viewf = 0.0;

%__________________________________________________________________________
function [viewf] = warning_ii(NO)
fprintf('*** WRONG VALUE OF NARG (NARG = %d) FOR NO = "<<NO<<"*** \n', NO);
viewf = 0.0;

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
  