function viewfactors
%    PROGRAM VIEWFACTORS
%  ******************************************************************************
%  *   Program to calculate view factors for all configurations in Appendix D   *
%  *   as well as between parallel and perpendicular rectangles (Fig. 4-15)     *
%  *   Input: NO =   viewfactor number (1-51)                                   *
%  *          NARG = number of arguments required for viewfactor                *
%  *          ARG(1..NARG) = arguments (in alphabetical order as given in App.D,*
%  *                         angles given in degrees)                           *
%  ******************************************************************************

x = 1;
while x~=0.
    NO = input('view factor number? (number in App. D or 98 for perppltf and 99 for parlpltf) \n');
    NARG = input('number of arguments? \n');
    for n=1:NARG
        ARG(n) = input('ARG("<<n<<")?= \n');
    end
    if(NO<1)
        fprintf('Illegal view factor number! \n');
    elseif(NO<52)
        F=view(NO,NARG,ARG);
    elseif(NO==98)
        if(NARG~=7)
            fprintf('Illegal number of arguments for PERPPLTF! \n');
        end
        X1=ARG(1);
        X2=ARG(2);
        Y1=ARG(3);
        Y2=ARG(4);
        Z1=ARG(5);
        Z2=ARG(6);
        Z3=ARG(7);
        F=perppltf(X1,X2,Y1,Y2,Z1,Z2,Z3);
        
    elseif(NO==99)
        if(NARG~=7)
            fprintf('Illegal number of arguments for PARLPLTF! \n');
        end
        X1=ARG(1);
        X2=ARG(2);
        X3=ARG(3);
        Y1=ARG(4);
        Y2=ARG(5);
        Y3=ARG(6);
        Z=ARG(7);
        F=parlpltf(X1,X2,X3,Y1,Y2,Y3,Z);
    else
        fprintf('Illegal view factor number! \n');
    end    
    fprintf('F = %f \n', F);
    x = input('hit 0 or enter to exit, anything other number to continue \n');
end


%==========================================================================
%
%==========================================================================
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
        
        RTX = sqrt(1.0 + X*X);
        RTY = sqrt(1.0 + Y*Y);
        RT  = sqrt(1.0 + X*X + Y*Y);
        viewF = (log(RTX*RTY/RT)     ...
                 + X*RTY*atan(X/RTY) ...
                 + Y*RTX*atan(Y/RTX) ...
                 - X*atan(X)- Y*atan(Y) ...
                )*2./(pi*X*Y);
        
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



function [viewf] = warning_i(NO)
fprintf('***ILLEGAL VALUE FOR NO (NO = %d), MUST BE 1<=NO>=51 *** \n', NO);
viewf = 0.0;

function [viewf] = warning_ii(NO)
fprintf('*** WRONG VALUE OF NARG (NARG = %d) FOR NO ="<<NO<<"*** \n', NO);
viewf = 0.0;

function [y] = pow(x, n)
y = x^n;




