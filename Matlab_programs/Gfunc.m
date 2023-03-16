function Gzn = Gfunc(y,Tz,teta0,tetatau0,y0,T0,Tk,nv,tol,alsred,Ref,koeftep,Ns)
sigma=5.668e-8;
te0=273.15; 
%temv = [20 50 100 250 500]; temv = temv+273.15; tepv = [0.045 0.047 0.051 0.065 0.096]; tepv = tepv*4184/3600; kotepv = koef(temv,tepv,length(tepv))';koeftep = polyval(kotepv,Tz);
tau0=alsred*y0;
tau=alsred*y;
N=(koeftep*alsred)/(4*sigma*(Tz^3));
%Tr=1;for k=1:Ns Tr=Tr*exp(-alsred*tol); Tr=Tr*(1-Rfl); end;Tr=Tr*exp(-tau0);Ab=1-Tr;Ab=0;lam0=2898/(T0+te0);pi=3.1415926535897932;kap=(lam0*alsred)/(4*pi);Rfl=((nv-1)^2+kap^2)/((nv+1)^2+kap^2);eps=(0.6e-8)*4184/3600; stch=eps/sigma;Ns=y0/tol;
stch=1-Ref;
R0=stch*sigma*((teta0*T0)^4)+Ref*sigma*((tetatau0*Tk)^4);
Rtau0=stch*sigma*((tetatau0*Tk)^4)+Ref*sigma*((teta0*T0)^4);
tm=(sigma*(Tz^4));
beta0=R0/tm;
betatau0=Rtau0/tm;
taus=alsred*y;
a=tau; 
tm=integroexpon(4,a); 
if (a==0)
    tm=1/3; 
end;
a=tau0; 
tr=integroexpon(4,a);
a=tau0-tau; 
tu=integroexpon(4,a); 
if (a==0) 
    tu=1/3; 
end;
G1 = beta0*(-tm+(tau/tau0)*tr+(1-tau/tau0)/3);
G2 = betatau0*((1-tau/tau0)*tr-tu+tau/(3*tau0));
G3 = 2*N*(teta0+(tau/tau0)*(tetatau0-teta0));
Gzn = (G1+G2+G3)/(2*N);
end