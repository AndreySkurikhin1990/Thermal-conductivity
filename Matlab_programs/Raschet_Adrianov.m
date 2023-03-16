format long g; alfs=0; alfs=SredGraf; dl=0; dl=dlvoVer53101; p3=length(dl); for k=1:p3 dl(k)=1e4/dl(k); end; als=0; als=trapz(dl,alfs)/(dl(p3)-dl(1));
knuSha=0; knuSha=preobMas(); npp=0; npp=Kramers_n(); ronu=0; nsh=0; nsha=0; nsh=Kramers_n_Sha(knuSha,dl); hi=0; hiSha=0; for k=1:p3    hiSha(k)=dl(k)*knuSha(k)/(4*pi); hi(k)=dl(k)*alfs(k)/(4*pi);nsha(k)=((nsh(k)^2+hiSha(k)^2)/(npp(k)^2+hi(k)^2))^0.5; ronu(k)=0.5+((nsha(k)-1)*(3*nsha(k)+1))/(6*(nsha(k)+1)^2)-(2*(nsha(k)^3)*(nsha(k)^2+2*nsha(k)-1))/((nsha(k)^2+1)*(nsha(k)^4-1)); ronu(k)=ronu(k)+(8*(nsha(k)^4)*((nsha(k)^4)+1)*log(nsha(k)))/((nsha(k)^2+1)*((nsha(k)^4-1)^2));ronu(k)=ronu(k)+(nsha(k)^2)*((nsha(k)^2-1)^2)*log((nsha(k)-1)/(nsha(k)+1))/((nsha(k)^2+1)^3);end; 
Refl1=real(trapz(dl,ronu)/(dl(p3)-dl(1))); Refl2=Refl1; a1=1-Refl1; a2=1-Refl2; ep=1e-20; n0=round(Nto/10);
x=0; y0=3e4; y0l=8e3; l=2e0; Nto=round(y0l/l); intex30=integroexpon(3,als*y0l);
if (abs(intex30)<ep) else if (isnan(intex30)) intex30=0; end; end;
sig=5.668e-8; te0=273.15; T0=585+te0; Tk=109+te0; q1=0; q2=0; q3=0; q4=0; qr=0; qob=0; koo=0;
te1=0; Tsr=(T0+Tk)/2; for k=1:Nto+1    te1(k)=(T0^4+(Tk^4-T0^4)*(k-1)/Nto)^0.25; end;
temv = [20 50 100 250 500]; tepv = [0.045 0.047 0.051 0.065 0.096]; p4=length(temv); for j=1:p4 temv(j) = temv(j)+te0; tepv(j)=tepv(j)*4184/3600; end; kotepv = koef(temv,tepv,p4)'; qm=0;
for k=1:Nto+1 koo(k)=x; x=x+l; end;
for k=1:Nto+1
x=koo(k);
intex3(k)=integroexpon(3,als*x);
intex3m(k)=integroexpon(3,als*(y0l-x));
if (abs(intex3(k))<ep) else if (isnan(intex3(k))) intex3(k)=0; end; end; if (abs(intex3m(k))<ep) else if (isnan(intex3m(k))) intex3m(k)=0; end; end;
psi1(k)=2*intex3(k)-4*intex30*intex3m(k)*Refl2;
psi1(k)=psi1(k)/(1-Refl1*Refl2*4*intex30^2);
psi2(k)=2*intex3m(k)-4*intex30*intex3(k)*Refl1;
psi2(k)=psi2(k)/(1-Refl1*Refl2*4*intex30^2);
ks=0; F1=0; F2=0; ko1=0; ko2=0; m=1; 
for j=1:Nto+1
if (j<k)
    F1(j)=integroexpon(2,als*(x-ks))+Refl1*psi1(k)*integroexpon(2,als*ks)-Refl2*psi2(k)*integroexpon(2,als*(y0l-ks)); ko1(j)=ks; if (abs(F1(j))<ep) else if (isnan(F1(j))) F1(j)=0; end; end;
else
    F2(m)=integroexpon(2,als*(ks-x))+Refl2*psi2(k)*integroexpon(2,als*(y0l-ks))-Refl1*psi1(k)*integroexpon(2,als*ks); ko2(m)=ks; if (abs(F2(m))<ep) else if (isnan(F2(m))) F2(m)=0; end; end; m=m+1;
end;
ks=ks+l;
end
disp(k);
p1=length(ko1); p2=length(ko2); pf1=0; for j=1:p1   pf1(j)=F1(j)*te1(j)^4; end; for j=1:p2  pf2(j)=F2(j)*te1(j)^4; end;
s1=0; for j=1:p1-1   s1=s1+(ko1(j+1)-ko1(j))*(pf1(j)+pf1(j+1))/2; end; s2=0; for j=1:p2-1   s2=s2+(ko2(j+1)-ko2(j))*(pf2(j)+pf2(j+1))/2; end;
q1(k)=psi1(k)*a1*sig*T0^4; q2(k)=-psi2(k)*a2*sig*Tk^4; q3(k)=2*als*sig*s1; q4(k)=-2*als*sig*s2; qr(k)=q1(k)+q2(k)+q3(k)+q4(k);
if (k<1e2) disp(q1(k)); disp(q2(k)); disp(q3(k)); disp(q4(k)); end;
if (k<Nto+1) koeftep = polyval(kotepv,te1(k)); qm(k)=-koeftep*(te1(k+1)-te1(k))/(koo(k+1)-koo(k)); else qm(k)=qm(k-1); end; qob(k)=qm(k)+qr(k);
end
Ref=Refl1*Refl2; tau=als*y0l; az1=Refl1*(1-exp(-2*tau))+(1-Ref*exp(-2*tau)); az1=az1/2/(1-Ref*exp(-4*tau)); az2=Refl2*(1-exp(-2*tau))+(1-Ref*exp(-2*tau)); az2=az2*exp(-2*tau)/2/(1-Ref*exp(-4*tau)); 
Aef=4*tau*(1-Ref*exp(-4*tau))-(Refl1+Refl2)*(1-exp(-2*tau))^2-2*(1-Ref*exp(-2*tau))*(1-exp(-2*tau)); Aef=Aef/4/tau^2/(1-Ref*exp(-4*tau));
te2=0; koeftep=polyval(kotepv,Tsr); qobs=(T0-Tk)*koeftep/y0l+Aef*sig*(T0^4-Tk^4);
for k=1:Nto+1
    x=koo(k); te2(k)=T0+(sig*(T0^4-Tk^4)/tau-qobs)*x/koeftep+sig*(T0^4-Tk^4)*(az1*exp(-2*als*x)-az2*exp(2*als*x))/2/als^2/y0l/koeftep;
end
p=plot(koo(n0:Nto-n0),qob(n0:Nto-n0),'-b');set(p,'LineWidth',2); hold on; grid on; xlabel({'Координата, мкм'}); ylabel({'Плотность теплового потока, Вт/(м*К)'}); title({'График зависимости q_r(x)'});