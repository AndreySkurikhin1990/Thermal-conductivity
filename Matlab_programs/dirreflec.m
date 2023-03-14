function dirreflec()

A = input('INPUT INDEX OF REFRACTION: [N, K]:  ');
AN = A(1); 
AK = A(2);
fprintf('N = %f, K = %f \n', AN, AK);
ISTEP = input('SPACING FOR OUTPUT, IN DEGREES:   ');


fprintf('     THETA     RHOS      RHOP      RHO       EPS \n');
PI2=2.E0*atan(1.E0);
j = 1;
for I = 1:ISTEP:90
    TH(j) = (I-1)*PI2/9.E1;
    PH(j) = (I-1);
    [RHOS(j),RHOP(j),RHO(j)] = fresnel(AN,AK,TH(j));
    EPS(j) = 1.E0 - RHO(j);
    fprintf(' %g %g %g %g %g \n' , ...
            PH(j), RHOS(j), RHOP(j),RHO(j), EPS(j));
    j = j + 1;
end
h = polar(repmat(TH',1,4), [RHOS', RHOP', RHO', EPS']);
set(h, 'LineWidth', 2);
legend('RHOS', 'RHOP', 'RHO', 'EPS');

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


