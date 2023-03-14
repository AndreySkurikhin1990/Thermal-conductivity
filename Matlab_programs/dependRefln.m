function ronu = dependRefln(n)
ref=0;
ref=0.5+((n-1)*(3*n+1))/(6*(n+1)^2);
ref=ref-(2*(n^3)*(n^2+2*n-1))/((n^2+1)*(n^4-1)); 
podln=(n-1)/(n+1);
podln=abs(podln);
ref=ref+(8*(n^4)*((n^4)+1)*log(n))/((n^2+1)*((n^4-1)^2));
ref=ref+(n^2)*((n^2-1)^2)*log(podln)/((n^2+1)^3);
ronu=ref;
end