% Function Emdiel


function [emdiel_] = emdiel(n)

% Calculates spectral, hemispherical emittance of dielectric, eq(3.82)
emdiel_ = (4.E0*n+2.E0)/3.E0/pow((n+1.E0),2)+2.E0*pow(n,3)*(n*n+2.E0*n-1.E0) ...
    /(n*n+1.E0)/(pow(n,4)-1.E0)-pow((n*(n*n-1.E0)),2)/pow((n*n+1.E0),3) ...
    *log((n-1.E0)/(n+1.E0)) -8.E0*pow(n,4)*(pow(n,4)+1.E0)*log(n)/ ...
    (n*n+1.E0)/pow((pow(n,4)-1.E0),2);


function y = pow(x,n)
y = x^n;
