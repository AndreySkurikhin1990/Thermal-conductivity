format long g; alfs=0; alfs=1e6*SredGraf; knuSha=0; knuSha=preobMas();  npp=0; npp=Kramers_n(); dl=0; dl=dlvoVer53101; p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k); end;
%ronu=0; nsh=0; nsha=0; nsh=Kramers_n_Sha(knuSha,dl); nv=0; hi=0; hiSha=0; np=0; nsham=0;
%for k=1:p  hiSha(k)=dl(k)*knuSha(k)/(4*pi); hi(k)=dl(k)*alfs(k)/(4*pi); nsha(k)=((nsh(k)^2+hiSha(k)^2)/(npp(k)^2+hi(k)^2))^0.5;ronu(k)=0.5+((nsha(k)-1)*(3*nsha(k)+1))/(6*(nsha(k)+1)^2)-(2*(nsha(k)^3)*(nsha(k)^2+2*nsha(k)-1))/((nsha(k)^2+1)*(nsha(k)^4-1)); ronu(k)=ronu(k)+(8*(nsha(k)^4)*((nsha(k)^4)+1)*log(nsha(k)))/((nsha(k)^2+1)*((nsha(k)^4-1)^2));ronu(k)=ronu(k)+(nsha(k)^2)*((nsha(k)^2-1)^2)*log((nsha(k)-1)/(nsha(k)+1))/((nsha(k)^2+1)^3);end;%p=length(dl); su=0; for j=1:p-1    su=su+(dl(j+1)-dl(j))*(nsha(j)+nsha(j+1))/2; end; %su=0; for j=1:p-1 su=su+(dl(j+1)-dl(j))*(ronu(j)+ronu(j+1))/2; end; Refl=su/(dl(p)-dl(1)); Refl=abs(real(Refl)); %su=0; for j=1:p-1   su=su+(dl(j+1)-dl(j))*(npp(j)+npp(j+1))/2; end; %su=0; for j=1:p-1 su=su+(dl(j+1)-dl(j))*(alfs(j)+alfs(j+1))/2; end; alf=su/(dl(p)-dl(1)); alf=real(alf);
%-------------
y0=3e4; l=1e0; y0l=1e3; 
for l=1:20
Ntt=0; Ntt=round(y0/y0l); Nto=round(y0l/l);Nit=5e1; n0=round(Nto/8); te0=273.15; 
koeftep=0.245; T0=1000+te0; Tk=235+te0; Tz=T0; %Tz=(T0+Tk)/2; 
%temv = [20 50 100 250 500]; temv = temv+te0; tepv = [0.045 0.047 0.051 0.065 0.096]; tepv = tepv*4184/3600; kotepv = koef(temv,tepv,length(tepv))'; koeftep = polyval(kotepv,Tz);
%no=0; al=0; dlv=0; no=BBIn(Tz-te0); j=1; tet=tempvtst(T0-te0,Tk-te0,y0l,Ntt); for k=no(1):no(2) al(j)=alfs(k); dlv(j)=dl(k); j=j+1; end;p=length(al); su=0; for j=1:p-1 su=su+(dlv(j+1)-dlv(j))*(al(j)+al(j+1))/2; end; %alf=su/(dlv(p)-dlv(1));
koo=0; ko=0; teta=0; tem=0; 
%-------------
for k=1:Ntt 
    koo(k,1)=ko; ko=ko+y0l; end;
for k=1:Ntt-1 
    koo(k,Nto+1)=koo(k+1,1); end; ko=0;
for k=1:Ntt 
    for j=2:Nto 
        ko=ko+l; koo(k,j)=ko; end; ko=koo(k,1); end; 
koo(Ntt,Nto+1)=y0; ep=1e-20;
Tna=T0; Tko=Tk; Tsr=Tz; teta(1,1)=T0/Tz; teta(Ntt,Nto)=Tk/Tz;
for k=2:Ntt 
    teta(k,1)=((T0^4+(Tk^4-T0^4)*(k-1)/Ntt)^0.25)/Tz; end;
for k=1:Ntt-1 
    teta(k,Nto+1)=teta(k+1,1); end;
for j=1:Ntt 
    Tna=teta(j,1)*Tz; Tko=teta(j,Nto+1)*Tz; Tsr=Tz; %Tsr=(Tna+Tko)/2;
    for k=2:Nto
    teta(j,k)=((Tna^4+(Tko^4-Tna^4)*(k-1)/Nto)^0.25)/Tsr; end; end;
%-------------
tetam=0; taum=0; koom=0; 
for k=1:Nto+1 
tetam(k)=teta(1,k); koom(k)=koo(1,k); end; %nsham=real(trapz(dl,nsha)/(dl(p)-dl(1)));
np=real(nsreddv(dl,npp,Tz)); Refl=0; 
alf=1/sprko(Tz);
taum=alf*koom; tau=alf*koo;
%-------------
tau0=alf*y0l; 
int3taus=0; int3taus(1)=1/2; q=1; for k=1:Nto+1 
    a=taum(k);  y=abs(exp(-a)/a); if (y<ep)  
        tm=0; else  
        tm=integroexpon(3,a); end; int3taus(q)=tm; q=q+1; end; %int4tau=0; int4tau(1)=1/3; q=1; for k=1:Nto+1      a=taum(k); y=abs(exp(-a)/a); if (y<ep)         tm=0; else         tm=(exp(-a)-a*int3taus(q))/3; end; int4tau(q)=tm; q=q+1;end; int4tau=0; int4tau(1)=1/3; q=1; for k=1:Nto+1      a=taum(k); y=abs(exp(-a)/a); if (y<ep)         tm=0; else tm=integroexpon(4,a); end; int4tau(q)=tm; q=q+1;end; %disp(int4tau);
int3tau0mins=0; q=1; for k=1:Nto+1 
    a=taum(k); a=tau0-a; y=abs(exp(-a)/a); if (y<ep) 
        tm=0; else tm=integroexpon(3,a); end; int3tau0mins(q)=tm; q=q+1; end; %int4tau0m=0; int4tau0m(Nto+1)=1/3; q=1; for k=1:Nto+1      a=taum(k); a=tau0-a; y=abs(exp(-a)/a); if (y<ep)         tm=0;  else tm=(exp(-a)-int3tau0mins(q)*a)/3; end; int4tau0m(q)=tm; q=q+1; end; %int4tau0m=0; int4tau0m(Nto+1)=1/3; q=1; for k=1:Nto+1  a=taum(k); a=tau0-a; y=abs(exp(-a)/a); if (y<ep)  tm=0;  else tm=integroexpon(4,a); end; int4tau0m(q)=tm; q=q+1; end; disp(int4tau0m);
%-------------
na=0; ko=0; ktn=0; ktk=0; for j=1:Nit    
    q=1; ktn=ktk;
for k=1:Nto+1
%for p=1:Nto+1 tetam(p)=teta(1,p); taum(p)=tau(1,p); koom(p)=koo(1,p); end;
m=0; r=0; Tna=tetam(n0)*Tz; Tko=teta(Nto-n0)*Tz; Tsr=(Tna+Tko)/2;
m=Gfunc(koom(k),Tz,Tna/Tz,Tko/Tz,y0l,Tna,Tko,np,l,alf,Refl, koeftep,Nto);
r=integrTau(koom(k),np,tetam,y0l,Tz,Nto,alf,int3taus,int3tau0mins,k);
m = m + r;
if (isnan(m))   
    m=0; continue; end; if (isinf(m))  
    m=0; continue; end; if (abs(m)>25e-1)  
    m=0; continue; end; if (abs(m)<4e-1)
    m=0; continue; end; %teta(1,q)=real(m); 
tetam(q)=real(m);q=q+1; end; 
for g=1:Nto  %teta(1,g)=(teta(1,g)+teta(1,g+1))/2; 
tetam(g)=(tetam(g)+tetam(g+1))/2; end; 
Tna=tetam(n0)*Tz; Tko=tetam(Nto-n0)*Tz; Tsr=(Tna+Tko)/2;
jvk=RaschTepPot(np,taum,y0l,l,tetam,Tna,Tko,Tz,alf, Nto,koeftep,Refl,y0l*alf);
ktk=jvk*((koom(Nto-n0)-koom(n0))*(1e-6)); ktk=ktk/(tetam(Nto-n0)*Tz-tetam(n0)*Tz); ktk=abs(real(ktk));
disp(j); for g=1:Nto+1 
    if isnan(tetam(g)) 
    else
        if isinf(tetam(g)) disp(g); end; end; end;
if (abs(ktk-ktn)/ktn<1e-3)
    break; end;
na=ko;
ko=trapz(koom,tetam);
if (abs((na-ko)/ko)<1e-4) 
    break; end;
end;
end
Tna=tetam(n0)*Tz; Tko=tetam(Nto-n0)*Tz; 
%sig=5.67051e-8; qe=koeftep*(Tna-Tko)/(((Nto-2*n0)/Nto)*y0l); qe=qe/(sig*(Tz^4)); qr=1-(Tko/Tz)^4; qr=qr/(0.75*alf*y0l+1/(1-Refl)+1/(1-Refl)-1); qe=qe-qr; disp(qe*y0l*sig*(Tz^3)/(1-Tko/Tz));
jv=RaschTepPot(np,taum,y0l,l,tetam,Tna,Tko,Tz,alf,Nto,koeftep,Refl,y0l*alf)
q=length(tetam); 
tetam=Tz*tetam-te0;
kt=jv*((koom(2*n0)-koom(n0))*(1e-6)); kt=kt/(tetam(2*n0)-tetam(n0));kt=abs(real(kt))
%for j=1:q koom(j)=koom(j)/1e3; end;
p=plot(1e-3*koom(n0+1:q-n0-1),tetam(n0+1:q-n0-1),'-b');set(p,'LineWidth',2); hold on; grid on; xlabel({'Координата, мм'}); ylabel({'Температура, °С'}); title({'График зависимости температуры от координаты'});