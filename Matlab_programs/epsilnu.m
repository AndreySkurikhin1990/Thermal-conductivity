function epsnu = epsilnu(n)
eps=(4*n+2)/3/((n+1)^2);
eps=eps+2*(n^3)*(n^2+2*n-1)/(n^2+1)/(n^4-1);
eps=eps-8*(n^4)*(n^4+1)*log(n)/(n^2+1)/((n^4-1)^2);
eps=eps-(n^2)*log((n-1)/(n+1))*((n^2-1)^2)/((n^2+1)^3);
epsnu=eps;
end