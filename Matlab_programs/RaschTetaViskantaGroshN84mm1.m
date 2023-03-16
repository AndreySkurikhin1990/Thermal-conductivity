format long g; ro=112; cp=1533; y0=3e4*1e-6; y0l=5e3*1e-6; tetan=0; teta=0; koo=0; tetam=0; koom=0; tau=0; taum=0; alfm=0;
alfs=0; alfs=1e6*SredGraf(); npp=0; npp=Kramers_n(); dl=0; dl=dlvoVer53101(); p=length(dl); te0=273.15; ep=1e-20; %vl=299792458; 
temvh=arrTemHigh(); temvh=temvh+te0; temvc=arrTemCold(); temvc=temvc+te0; 
tepo=arrTepPot84(); qv=(1/13.85e-4)*tepo; ts=0; pt=length(temvc); 
for k=1:pt
    ts=ts+temvh(k)+temvc(k);
end
ts=ts/2/pt;
for k=1:p  
    dl(k)=1e-2/dl(k); %fr(k)=vl/dl(k); 
end;
%-------------
Refl=Refl_Sha(dl,npp,p,ts);
%Refl=0;
for k=1:pt
tepv(k)=qv(k)*y0/(temvh(k)-temvc(k));
end
temvs=(1/2)*(temvc+temvh); a=tepv/cp/ro;
ktkm=0; ktkoh=0; ktkoc=0; wh=1; wc=1; tevho=0; tevco=0;
for w=1:pt
    f=1; ot=0; Tk=temvc(w); T0=temvh(w); koeftep=tepv(w);
for le=1e2:1e2
    l=le*1e-6; Ntt=ceil(y0/y0l);Nto=ceil(y0l*1e6/le);
    Nit=3e2; n0=1; %Tz=(T0+Tk)/2; 
    Tz=T0; %temv = [20 50 100 250 500]; temv = temv+te0; tepv = [0.045 0.047 0.051 0.065 0.096]; tepv = tepv*4184/3600; kotepv = koef(temv,tepv,length(tepv))'; koeftep = polyval(kotepv,Tz);no=0; al=0; dlv=0; no=BBIn(Tz-te0); j=1; tet=tempvtst(T0-te0,Tk-te0,y0l,Ntt); for k=no(1):no(2) al(j)=alfs(k); dlv(j)=dl(k); j=j+1; end;
koo=vychKoor(y0,y0l,l,Nto,Ntt);
%-------------
np=real(nsreddv(dl,npp,Tz));
alf=1e6/sprko(Tz);
%alf=knusreddv(alfs,npp,dl,Tz);
tau=alf*koo;
%-------------
Tna=T0; Tko=Tk; Tsr=Tz;
tetan=(1/Tsr)*nachPribTemN(Tna,Tko,l,y0,np,koeftep,alf);
%tetan=(1/Tsr)*(nachPribTem(T0-te0,l,a(w),y0,y0l)+te0);
%-------------
q=1; for j=1:Ntt 
    Tsr=0;
    for k=1:Nto 
        teta(j,k)=tetan(q);    
        q=q+1; end; 
    Tsr=teta(j,1)*Tz;
    %alf=koefpoglPatch(koo(j,1),Tsr); 
    alf=knusreddv(alfs,npp,dl,Tsr);
    %alf=1e6/sprko(Tsr);
    alfm(j)=alf;
end;
teta(Ntt,Nto+1)=T0/Tz;

for k=1:Ntt-1
    teta(k,Nto+1)=teta(k+1,1);
end;

%-------------
ktn=koeftep; ktk=ktn;
for j=1:Nit    
    q=1; 
    ktn=ktk;
%-------------
for yc=1:Ntt
    
for ri=1:Nto+1 
    tetam(ri)=teta(yc,ri); 
    taum(ri)=koo(yc,ri)*alfm(yc); 
    koom(ri)=koo(yc,ri); 
end;
tau0=alfm(Ntt)*koo(Ntt,Nto+1);

int3taus=intexpfunc3(0,ep,taum,Nto);
int3tau0mins=intexpfunc3(tau0,ep,taum,Nto);
%------------- 

for k=1:Nto+1
n1=length(tetam); m=0; r=0; Tna=tetam(n0)*Tz; Tko=tetam(n1)*Tz; Tsr=Tna; %Tsr=(Tna+Tko)/2;
m=Gfunc(koom(k),Tsr,Tna/Tsr,Tko/Tsr,y0,Tna,Tko,np,l,alfm(yc),Refl,ktn,Nto);
r=integrTau(koom(k),np,tetam,y0,Tsr,Nto,alfm(yc),int3taus,int3tau0mins,k,ktn);
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
%disp(j);
ot=abs(ktk-ktn);
%disp(ot);
if (ot<1e-3) 
    break; end; 
if (j>2e2) 
    f=0; break; end;
if (ktk<0) 
    ktk=1; continue; end;

for k=1:Ntt
    teta(yc,k)=tetam(k);
end;

end;
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
ktkm(w)=ktk;
end;
koo=1*koo
teta=1*teta
tau=1*tau
ktkm=1*ktkm
tepv=1*tepv
%sig=5.67051e-8; qe=koeftep*(Tna-Tko)/(((Nto-2*n0)/Nto)*y0l); qe=qe/(sig*(Tz^4)); qr=1-(Tko/Tz)^4; qr=qr/(0.75*alf*y0l+1/(1-Refl)+1/(1-Refl)-1); qe=qe-qr; disp(qe*y0l*sig*(Tz^3)/(1-Tko/Tz));
%q=length(tetam); tetam=Tz*tetam-te0;for j=1:q koom(j)=koom(j)/1e3; end;
wh=1; wc=1; for k=1:length(ktkm)
    if (ktkm(k)<tepv(k))
        if (rem(k,2)==0)
            ktkoh(wh)=ktkm(k);
            tevho(wh)=temvs(k);
            wh=wh+1;
        else
            ktkoc(wc)=ktkm(k);
            tevco(wc)=temvs(k);
            wc=wc+1;
        end
    end
end
ktkoh=1*ktkoh
ktkoc=1*ktkoc
tevho=1*tevho
tevco=1*tevco
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
if isnan(tepv1)
else
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
end