format long g; lm = [3 5 10 20 25 30 50 60 75 100]; nlm = length(lm);
ro=112; cp=650/3; alfs=0; alfs=1e6*SredGraf; knuSha=0; knuSha=preobMas(); 
npp=0; npp=Kramers_n(); dl=0; dl=dlvoVer53101; p=length(dl); te0=273.15; ep=1e-20; %vl=299792458; 
koeftep=0.156; a=koeftep/ro/cp; T0=603+te0; Tk=109+te0;
for k=1:p  
    dl(k)=1e-2/dl(k); %fr(k)=vl/dl(k); 
end;
%-------------
Refl=Refl_Sha(alfs,dl,npp,p,knuSha)
%Refl=0;
le=1; u=1; 
for w=1:nlm
    l=lm(w); y0=3e4; y0l=1e3; Ntt=ceil(y0/y0l); Nto=ceil(y0l/l); Nit=1e3; n0=1; Tz=(T0+Tk)/2; %temv = [20 50 100 250 500]; temv = temv+te0; tepv = [0.045 0.047 0.051 0.065 0.096]; tepv = tepv*4184/3600; kotepv = koef(temv,tepv,length(tepv))'; koeftep = polyval(kotepv,Tz);no=0; al=0; dlv=0; no=BBIn(Tz-te0); j=1; tet=tempvtst(T0-te0,Tk-te0,y0l,Ntt); for k=no(1):no(2) al(j)=alfs(k); dlv(j)=dl(k); j=j+1; end;
%koo=vychKoor(y0,y0l,l,Nto,Ntt);
%-------------
np=real(nsreddv(dl,npp,Tz)); 
alf=1/sprko(Tz); %alf=knusreddv(alfs,npp,dl,Tz);
%tau=alf*koo; 
tau0=alf*y0;
%-------------
Tna=T0; Tko=Tk; Tsr=Tz;
tetam=(1/Tsr)*(nachPribTem(T0-te0,l,a,y0,y0l)+te0);
%-------------
koom=0:l:y0; taum=alf*koom; q=1; 
%j=round(Ntt);
%for j=1:Ntt
%for k=1:Nto+1
%    tetam(q)=teta(j,k);
%    koom(q)=koo(j,k); 
%    q=q+1; 
%end;
%end;
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
disp(ktk);
end;
%disp(n0);
%disp(n1);
%disp(Tz*tetam-te0);
%Tna=tetam(n0)*Tz;
%Tko=tetam(n1)*Tz;
%disp(Tna-te0);
%disp(Tko-te0);
%sig=5.67051e-8; qe=koeftep*(Tna-Tko)/(((Nto-2*n0)/Nto)*y0l); qe=qe/(sig*(Tz^4)); qr=1-(Tko/Tz)^4; qr=qr/(0.75*alf*y0l+1/(1-Refl)+1/(1-Refl)-1); qe=qe-qr; disp(qe*y0l*sig*(Tz^3)/(1-Tko/Tz));
%q=length(tetam); 
%tetam=Tz*tetam-te0;
%for j=1:q koom(j)=koom(j)/1e3; end;
%p=plot(1e-3*koom,tetam,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Координата, мм'}); 
%ylabel({'Температура, °С'}); 
%title({'График зависимости температуры от координаты'});