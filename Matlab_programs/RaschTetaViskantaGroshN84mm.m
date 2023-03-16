format long g; ro=112; cp=1533; y0=3e4; y0l=5e3; 
teta=0; koo=0; tetam=0; koom=0; tau=0; taum=0;
alfs=0; alfs=1e6*SredGraf(); npp=0; npp=Kramers_n(); dl=0; dl=dlvoVer53101(); p=length(dl); te0=273.15; ep=1e-20; %vl=299792458; 
temvh=arrTemHigh(); temvh=temvh+te0; temvc=arrTemCold(); temvc=temvc+te0; tepo=arrTepPot84(); qv=(1/13.85e-4)*tepo; ts=0; pt=length(temvc); 
for k=1:pt
    ts=ts+temvh(k)+temvc(k);
end
ts=ts/2/pt
for k=1:p  
    dl(k)=1e-2/dl(k); %fr(k)=vl/dl(k); 
end;
%-------------
Refl=Refl_Sha(dl,npp,p,ts)
%Refl=0;
for k=1:pt
tepv(k)=qv(k)*(y0*1e-6)/(temvh(k)-temvc(k));
end
temvs=(1/2)*(temvc+temvh); a=tepv/cp/ro;
ktkoh=0; ktkoc=0; wh=1; wc=1; tevho=0; tevco=0;
for w=1:pt
    f=1; ot=0; T0=temvc(w); Tk=temvh(w); koeftep=tepv(w);
for l=5e2:5e2
    Ntt=ceil(y0/y0l); Nto=ceil(y0l/l); Nit=2e2; n0=1; Tz=(T0+Tk)/2 %Tz=T0 %temv = [20 50 100 250 500]; temv = temv+te0; tepv = [0.045 0.047 0.051 0.065 0.096]; tepv = tepv*4184/3600; kotepv = koef(temv,tepv,length(tepv))'; koeftep = polyval(kotepv,Tz);no=0; al=0; dlv=0; no=BBIn(Tz-te0); j=1; tet=tempvtst(T0-te0,Tk-te0,y0l,Ntt); for k=no(1):no(2) al(j)=alfs(k); dlv(j)=dl(k); j=j+1; end;
%koo=vychKoor(y0,y0l,l,Nto,Ntt);
%-------------
np=real(nsreddv(dl,npp,Tz))
%alf=1/sprko(Tz);
alf=knusreddv(alfs,npp,dl,Tz)*1e-6%tau=alf*koo; 
tau0=alf*y0;
%-------------
Tna=T0; Tko=Tk; Tsr=Tz;
tetam=(1/Tsr)*nachPribTemN(Tna,Tko,l,y0,np,koeftep,alf);
%tetam=(1/Tsr)*(nachPribTem(T0-te0,l,a(w),y0,y0l)+te0);
%-------------
koom=0:l:y0; taum=alf*koom; 
%q=1; j=round(Ntt);for j=1:Ntt for k=1:Nto+1    tetam(q)=teta(j,k);    koom(q)=koo(j,k);     q=q+1; end;end;
Ntt=1; Nto=ceil(y0/l);
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
n1=length(tetam); m=0; r=0; Tna=tetam(n0)*Tz; Tko=tetam(n1)*Tz; Tsr=(Tna+Tko)/2;
m=Gfunc(koom(k),Tz,Tna/Tz,Tko/Tz,y0,Tna,Tko,np,l,alf,Refl,ktn,Nto);
r=integrTau(koom(k),np,tetam,y0,Tz,Nto,alf,int3taus,int3tau0mins,k);
m=m+r; m=ProvAdek(m); 
if (m==0) 
else tetam(q)=m; q=q+1; 
end; end;
%-------------
for k=2:Nto+1 
    tetam(k)=(tetam(k)+tetam(k-1))/2; 
end;
%-------------
ktk=abs(RaschTepPotN(np,taum,y0,l,tetam,Tna,Tko,Tz,alf,Nto,koeftep,Refl,tau0));
disp(j);
ot=abs(ktk-ktn)
if (ot<1e-4) 
    break; end; 
if (j>1e2) 
    f=0; break; end;
if (ktk<0) 
    ktk=1; continue; end;
end;
end;
if (f==1) 
    if (ktk>0)
        if (ktk<koeftep)
    if (temvs(w)>ts)
        tevho(wh)=temvs(w);
        ktkoh(wh)=ktk;
        wh=wh+1;
    disp(ktk);
    else
        tevco(wc)=temvs(w);
        ktkoc(wc)=ktk;
        wc=wc+1;
        disp(ktk);
    end;
        end;
    end;
end;
end;
disp(tevho);
disp(ktkoh);
disp(tevco);
disp(ktkoc);
%sig=5.67051e-8; qe=koeftep*(Tna-Tko)/(((Nto-2*n0)/Nto)*y0l); qe=qe/(sig*(Tz^4)); qr=1-(Tko/Tz)^4; qr=qr/(0.75*alf*y0l+1/(1-Refl)+1/(1-Refl)-1); qe=qe-qr; disp(qe*y0l*sig*(Tz^3)/(1-Tko/Tz));
%q=length(tetam); tetam=Tz*tetam-te0;for j=1:q koom(j)=koom(j)/1e3; end;
p=length(tevho); tem1=0; tem2=0; tepv1=0; tepv2=0;
for k=1:p
tem1=tem1+tevho(k); tepv1=tepv1+ktkoh(k)*tevho(k);
end
tepv1=tepv1/tem1
tem1=tem1/p
p=length(tevco);
for k=1:p
tem2=tem2+tevco(k); tepv2=tepv2+ktkoc(k)*tevco(k);
end
disp(temvs-te0);
tepv2=tepv2/tem2
tem2=tem2/p
kn=(tepv2-tepv1)/(tem2-tem1);
sc=tepv1-kn*tem1;
ktep=0; te=3e2:1:1e3; te=te+te0;
for k=1:length(te)
    ktep(k)=kn*te(k)+sc;
end
pl=plot(te,ktep,'-b');
set(pl,'LineWidth',2); 
hold on; grid on; 
xlabel({'Температура, К'}); 
ylabel({'Коэффициент кондуктивной теплопроводности, Вт/(м*К)'}); 
title({'График зависимости кондуктивного КТП от температуры'});