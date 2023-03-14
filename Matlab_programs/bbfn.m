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





