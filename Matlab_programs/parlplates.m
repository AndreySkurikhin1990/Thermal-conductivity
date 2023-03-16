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
epsmin=1e-15;
NARG = 3;
ARG(1) = X3; ARG(2) = Y3; ARG(3) = C;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = A*view(38,NARG,ARG);
end
disp('1');
F=F' %1
%--------
ARG(2) = Y2;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG);
end
disp('2');
F=F' %2
%--------
ARG(2) = Y3-Y1;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG);
end
disp('3');
F=F' %3
%--------
ARG(2) = Y2-Y1;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    t=view(38,NARG,ARG)
    F = F+A*view(38,NARG,ARG);
end
disp('4');
F=F' %4
%--------
ARG(1) = X2;
ARG(2) = Y3;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG);
end
disp('5');
F=F' %5
%--------
ARG(2) = Y2;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = F+A*view(38,NARG,ARG);
end
disp('6');
F=F' %6
%--------
ARG(2) = Y3-Y1;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = F+A*view(38,NARG,ARG);
end
disp('7');
F=F' %7
%--------
ARG(2) = Y2-Y1;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG);
end
disp('8');
F=F' %8
%--------
ARG(1) = X3-X1;
ARG(2) = Y3;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG);
end
disp('9');
F=F' %9
%--------
ARG(2) = Y2;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = F+A*view(38,NARG,ARG);
end
disp('10');
F=F' %10
%--------
ARG(2) = Y3-Y1;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = F+A*view(38,NARG,ARG);
end
disp('11');
F=F' %11
%--------
ARG(2) = Y2-Y1;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG);
end
disp('12');
F=F' %12
%--------
ARG(1) = X2-X1;
ARG(2) = Y3;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = F+A*view(38,NARG,ARG);
end
disp('13');
F=F' %13
%--------
ARG(2) = Y2;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG);
end
disp('14');
F=F' %14
%--------
ARG(2) = Y3-Y1;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = F-A*view(38,NARG,ARG);
end
disp('15');
F=F' %15
%--------
ARG(2) = Y2-Y1;
A = ARG(1)*ARG(2)
if (abs(A)>epsmin)
    F = F+A*view(38,NARG,ARG);
end
disp('16');
F=F' %16
%--------
parlpltf = F/(4.*(X1*Y1));
end