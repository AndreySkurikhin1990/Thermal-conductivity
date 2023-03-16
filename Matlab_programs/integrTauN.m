function intgr = integrTauN(y,npp,teta,y0,Tz,Nt,alsr,nom,koeftep)
tau0=alfm(Ntt)*koo(Ntt,Nto+1);
int3taus=intexpfunc3(0,ep,taum,Nto);
int3tau0mins=intexpfunc3(tau0,ep,taum,Nto);
sigma=5.668e-8; 
%temv = [20 50 100 250 500]; %te0=273.15; %temv = temv+te0; %tepv = [0.045 0.047 0.051 0.065 0.096];%tepv = tepv*4184/3600; %kotepv = koef(temv,tepv,length(tepv))'; %koeftep = polyval(kotepv,Tz); 
%------
N=(koeftep*alsr)/(4*sigma*(Tz^4));
tau0=alsr*y0; 
htaush=tau0/Nt; 
tau=alsr*y; 
E3=0; 
inte=0; 
ko=0; 
q=1; 
taush=0; 
tm=0; 
tr=0; 
integri=0;
%-------
for k=1:Nt+1
taush=(k-1)*htaush;
if (nom == k) 
    E3=1/2; 
else
    E3=-(integroexpon(3,abs(tau-taush))); 
end;
tm=int3taus(k); 
E3=E3+tm; 
tr=tm; 
tm=int3tau0Mins(k);
E3=E3+(tau/tau0)*(tm-tr); 
E3=(npp^2)*E3*(teta(q)^4); 
inte(q)=real(E3);
ko(q)=taush; 
q=q+1;
E3=0;
end;
integri=0; 
su=0; 
q=length(inte);
for k=1:q-1
su=su+(ko(k+1)-ko(k))*(inte(k)+inte(k+1))/2;
end
integri=su;
integri=integri/(2*N);
%disp(integri);
intgr=integri;
end