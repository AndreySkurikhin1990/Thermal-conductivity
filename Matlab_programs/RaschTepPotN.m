function koefte = RaschTepPotN(nv,tau,y0,tol,tem,T0,Tk,Tz,alsred,Ns,koeftep,Ref,tau0)
jv=-koeftep*(Tk-T0)/y0; jv=abs(real(jv));
cons=abs((Tk-T0)/y0);
sigma=5.668e-8; 
stch=1-Ref;
psi12=1; psi21=psi12; %psi12=2*g; psi21=psi12;
tem=Tz*tem; %alsred=IntegrDlVol;eps=(0.6e-8)*4184/3600; stch=eps/sigma;temv = [20 50 100 250 500]; temv = temv+273.15; tepv = [0.045 0.047 0.051 0.065 0.096];tepv=tepv*4184/3600; kotepv = koef(temv,tepv,length(tepv))';epssh=(9+12)/2; ro=2.6e3; vol=(0.5e-3)^3; hi=-0.5e-3; pi=3.1415926535897932;mu=1+4*pi*hi*ro*vol; epssh=mu*epssh; alsred=IntegrDlVol;lam0=(3000+1862)/2; lam0=0.01/lam0; kap=(lam0*alsred)/(4*pi);nkao=(epssh+kap^2)^0.5; epsdsh=2*nkao*kap; Rs=(70+87)/200;Dis=(4*((1+Rs)^2)-4*((1-Rs)^2)*(1+kap^2))^0.5; n1=(2*(1+Rs)-Dis)/(2*(1-Rs)); n2=(2*(1+Rs)+Dis)/(2*(1-Rs));nkao=(n1+nkao)/2; nv=(1.538+1.542)/2; Rfl=((nv-1)^2+kap^2)/((nv+1)^2+kap^2);Ns=round(y0/tol); Tr=1; Ab=1-Tr;for k=1:Ns Tr=Tr*exp(-alsred*tol); Tr=Tr*(1-Rfl); end; Tr=Tr*exp(-alsred*y0);
R0=stch*sigma*(T0^4)+psi21*Ref*sigma*(Tk^4);
Rtau0=stch*sigma*(Tk^4)+psi12*Ref*sigma*(T0^4);
ppi=denTherFlux(nv,tol,tem,Tk,alsred,Ns,tau0,stch,sigma,tau,R0,Rtau0); ppi=abs(real(ppi));
kotet=abs(jv-ppi); kotet=kotet/cons;
koefte=kotet;
end

function qtp = denTherFlux(nv,tol,tem,Tk,alsred,Ns,tau0,stch,sigma,tau,R0,Rtau0)
ep=1e-18;
%------
int3extau0=vychExpFu(3,tau0,ep);
int4extau0=vychExpFu(4,tau0,ep);
qss=2*(R0*((1-stch)*int3extau0+(1/tau0)*int4extau0-1/(3*tau0)));
qss=qss+2*(Rtau0*(-int4extau0/tau0-0.5+(1/(3*tau0))));
qss=qss+stch*sigma*(Tk^4);
%-------
q=1; inte=0; integ=0; tm=0; y=0; tr=0; tu=0;
%-------
for k=1:Ns+1
taus=alsred*y;
a=tau0-taus;
tm=vychExpFu(2,a,ep);
tr=vychExpFu(3,a,ep);
tu=vychExpFu(3,taus,ep);
inte=(nv^2)*((1-stch)*tm+(tr-tu)/tau0)*sigma*(tem(q)^4);
integ(q)=inte;
y=y+tol;
q=q+1;
end;
for k=1:Ns+1
    integ(k)=proveAdec(integ(k));
end
%disp(integ);
integr=2*trapz(tau,integ);
qtp=abs(integr);
end

function vief = vychExpFu(n,a,er)
g=exp(-a)/a;
if (g<er) 
    g=0; 
else
    g=integroexpon(n,a);
end;
if (a==0)
    g=1/(n-1); 
end
vief=g;
end

function mk = proveAdec(m)
if (isnan(m))   
    m=0; 
end; 
if (isinf(m))  
    m=0;  
end; 
mk=m;
end