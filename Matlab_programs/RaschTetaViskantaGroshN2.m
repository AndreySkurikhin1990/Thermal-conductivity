format long g; ro=112; cp=1588; alfs=0; alfs=1e6*SredGraf; knuSha=0; knuSha=preobMas(); y0=3e4; y0l=1e3; 
npp=0; npp=Kramers_n(); dl=0; dl=dlvoVer53101; p=length(dl); te0=273.15; ep=1e-20; %vl=299792458; 
temvc=[109 235 109 238]; temvc=temvc+te0; temvh=[585 1000 603 1000]; temvh=temvh+te0; tepo=[3.9596 8.6377 3.5614925 7.123]; qv=(1e4/13.85)*tepo;
for k=1:4
tepv(k)=qv(k)*(y0*1e-6)/(temvh(k)-temvc(k));
end
temvs=(1/2)*(temvc+temvh); tepv1=(tepv(1)+tepv(3))/2; tepv2=(tepv(2)+tepv(4))/2; tem1=(temvs(1)+temvs(3))/2; tem2=(temvs(2)+temvs(4))/2; k1=(tepv2-tepv1)/(tem2-tem1); k2=tepv2-k1*tem2; a=tepv/(ro*cp); 
for k=1:p  
    dl(k)=1e-2/dl(k); %fr(k)=vl/dl(k); 
end;
%-------------
Refl=Refl_Sha(alfs,dl,npp,p,knuSha);
%Refl=0;
ktko=0; for w=1:4
    T0=temvc(w); Tk=temvh(w); koeftep=tepv(w);
for l=1e2:1e2
    Ntt=ceil(y0/y0l); Nto=ceil(y0l/l); Nit=1e3; n0=1; Tz=(T0+Tk)/2; %temv = [20 50 100 250 500]; temv = temv+te0; tepv = [0.045 0.047 0.051 0.065 0.096]; tepv = tepv*4184/3600; kotepv = koef(temv,tepv,length(tepv))'; koeftep = polyval(kotepv,Tz);no=0; al=0; dlv=0; no=BBIn(Tz-te0); j=1; tet=tempvtst(T0-te0,Tk-te0,y0l,Ntt); for k=no(1):no(2) al(j)=alfs(k); dlv(j)=dl(k); j=j+1; end;
%koo=vychKoor(y0,y0l,l,Nto,Ntt);
%-------------
np=real(nsreddv(dl,npp,Tz)); 
alf=1/sprko(Tz); %alf=knusreddv(alfs,npp,dl,Tz);tau=alf*koo; 
tau0=alf*y0;
%-------------
Tna=T0; Tko=Tk; Tsr=Tz;
tetam=(1/Tsr)*(nachPribTem(T0-te0,l,a(w),y0,y0l)+te0);
%-------------
koom=0:l:y0; taum=alf*koom; 
%q=1; j=round(Ntt);
%for j=1:Ntt
%for k=1:Nto+1
%    tetam(q)=teta(j,k);
%    koom(q)=koo(j,k); 
%    q=q+1; 
%end;end;
Ntt=1; Nto=ceil(y0/l); n1=Nto+1;
%------------- 
int3taus=intexpfunc3(0,ep,taum,Nto);
int3tau0mins=intexpfunc3(tau0,ep,taum,Nto);
%-------------
ktn=0; ktk=0;
for j=1:Nit    
    q=1; 
    ktn=ktk;
%-------------
for k=1:Nto+1
%for p=1:Nto+1 tetam(p)=teta(1,p); taum(p)=tau(1,p); koom(p)=koo(1,p); end;
m=0; r=0; Tna=tetam(n0)*Tz; Tko=tetam(n1)*Tz; Tsr=(Tna+Tko)/2;
m=Gfunc(koom(k),Tz,Tna/Tz,Tko/Tz,y0,Tna,Tko,np,l,alf,Refl,koeftep,Nto);
r=integrTau(koom(k),np,tetam,y0,Tz,Nto,alf,int3taus,int3tau0mins,k);
m=m+r; m=ProvAdek(m); 
if (m==0) 
else tetam(q)=m; q=q+1; 
end; 
end;
%-------------
for k=2:Nto+1 
    tetam(k)=(tetam(k)+tetam(k-1))/2; 
end;
%-------------
Tna=tetam(n0)*Tz; Tko=tetam(n1)*Tz; Tsr=(Tna+Tko)/2;
ktk=RaschTepPot(np,taum,y0,l,tetam,Tna,Tko,Tz,alf,Nto,koeftep,Refl,tau0);
disp(j); 
if (abs(ktk-ktn)/ktn<1e-3)    
    break; 
end; 
end;
end;
ktko(w)=ktk;
end;
%sig=5.67051e-8; qe=koeftep*(Tna-Tko)/(((Nto-2*n0)/Nto)*y0l); qe=qe/(sig*(Tz^4)); qr=1-(Tko/Tz)^4; qr=qr/(0.75*alf*y0l+1/(1-Refl)+1/(1-Refl)-1); qe=qe-qr; disp(qe*y0l*sig*(Tz^3)/(1-Tko/Tz));
%q=length(tetam); 
%tetam=Tz*tetam-te0;
%for j=1:q koom(j)=koom(j)/1e3; end;
ktko(1)=(ktko(1)*temvs(1)+ktko(3)*temvs(3))/(temvs(1)+temvs(3)); ktko(2)=(ktko(2)*temvs(2)+ktko(4)*temvs(4))/(temvs(2)+temvs(4)); temvs(1)=(temvs(1)+temvs(3))/2; temvs(2)=(temvs(2)+temvs(4))/2;
pl=plot(temvs(1:2)-te0,ktko(1:2),'-b');
set(pl,'LineWidth',2); 
hold on; grid on; 
xlabel({'Температура, °С'}); 
ylabel({'Коэффициент кондуктивной теплопроводности, Вт/(м*К)'}); 
title({'График зависимости кондуктивного КТП от температуры'});