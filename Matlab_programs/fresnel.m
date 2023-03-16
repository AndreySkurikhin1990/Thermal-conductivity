function [RHOS, RHOP, RHO] = fresnel(AN, AK, TH)  

C1 = AN^2 - AK^2 - (sin(TH))^2;
P2 = .5E0*(sqrt(C1^2 + (2.E0*AN*AK)^2) + C1);
Q2 = P2 - C1;
CT = cos(TH);
C2 = (AN^2-AK^2)*CT;
C3 = 2.E0*AN*AK*CT;
P  = sqrt(P2);
Q  = sqrt(Q2);
RHOP = ((P-C2)^2 + (Q-C3)^2)/((P+C2)^2 + (Q+C3)^2);
RHOS = ((CT-P)^2 + Q2)/((CT+P)^2 + Q2);
RHO  =  .5E0*(RHOP + RHOS);


