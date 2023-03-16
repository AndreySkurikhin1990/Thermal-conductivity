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
F = 0.;
ARG(1) = Y2;
ARG(2) = Z3;
ARG(3) = X2;
A = X2*Z3;
if(abs(A*Y2)>1E-6)
    F=A*view(39,NARG,ARG);
end
ARG(1) = Y1;
A = X2*Z3;
if(abs(A*Y1)>1E-6)
    F = F-A*view(39,NARG,ARG);
end
ARG(1) = Y2;
ARG(3) = X1;
A = X1*Z3;
if(abs(A*Y2)>1E-6)
    F = F-A*view(39,NARG,ARG);
end
ARG(1) = Y1;
A = X1*Z3;
if(abs(A*Y1)>1E-6)
    F = F+A*view(39,NARG,ARG);
end
ARG(1) = Y2;
ARG(2) = Z2;
A = X1*Z2;
if(abs(A*Y2)>1E-6)
    F=F+A*view(39,NARG,ARG);
end
ARG(1) = Y1;
A = X1*Z2;
if(abs(A*Y1)>1E-6)
    F = F-A*view(39,NARG,ARG);
end
ARG(1) = Y2;
ARG(3) = X2;
A = X2*Z2;
if(abs(A*Y2)>1E-6)
    F = F-A*view(39,NARG,ARG);
end
ARG(1) = Y1;
A = X2*Z2;
if(abs(A*Y1)>1E-6)
    F = F+A*view(39,NARG,ARG);
end
ARG(1) = Y2;
ARG(2) = (Z3-Z1);
A = X2*(Z3-Z1);
if(abs(A*Y2)>1E-6)
    F = F-A*view(39,NARG,ARG);
end
ARG(1) = Y1;
A = X2*(Z3-Z1);
if(abs(A*Y1)>1E-6)
    F = F+A*view(39,NARG,ARG);
end
ARG(1) = Y2;
ARG(3) = X1;
A = X1*(Z3-Z1);
if(abs(A*Y2)>1E-6)
    F = F+A*view(39,NARG,ARG);
end
ARG(1) = Y1;
A = X1*(Z3-Z1);
if(abs(A*Y1)>1E-6)
    F = F-A*view(39,NARG,ARG);
end
ARG(1) = Y2;
ARG(2) = (Z2-Z1);
ARG(3) = X2;
A = X2*(Z2-Z1);
if(abs(A*Y2)>1E-6)
    F = F+A*view(39,NARG,ARG);
end
ARG(1) = Y1;
A = X2*(Z2-Z1);
if(abs(A*Y1)>1E-6)
    F = F-A*view(39,NARG,ARG);
end
ARG(1) = Y2;
ARG(3) = X1;
A = X1*(Z2-Z1);
if(abs(A*Y2)>1E-6)
    F = F-A*view(39,NARG,ARG);
end
ARG(1) = Y1;
A = X1*(Z2-Z1);
if(abs(A*Y1)>1E-6)
    F = F+A*view(39,NARG,ARG);
end
perppltf = F/(2.*(X2-X1)*Z1);
