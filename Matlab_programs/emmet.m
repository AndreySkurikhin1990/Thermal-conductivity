% function emmet

function emmet_ = emmet(n, k)
% Calculates spectral, hemispherical emittance of a metal, eq(3.77)
% (accurate for k>=n only; do not use for dielectrics (k=0)!!


c1=n*n+k*k;
c2=n*n-k*k;
c3=(n+1)*(n+1)+k*k;
%if (k<1.E-5)
   % k = 1.E-5;  % guard against k=0
%end
emparl = 8.E0*n/c1*(1.E0-n/c1*log(c3)+c2/(k*c1)*atan(k/(n+1.E0)));
emperp = 8.E0*n*(1.E0-n*log(c3/c1)+c2/k*atan(k/(n+c1)));

emmet_ = .5E0*(emparl+emperp);


