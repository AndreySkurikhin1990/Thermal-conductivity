epssh=(9+12)/2; ro=2.6e3; vol=(0.5e-3)^3; hi=-0.5e-3; pi=3.1415926535897932;
mu=1+4*pi*hi*ro*vol; epssh=mu*epssh; alsred=IntegrDlVol; lam0=(3000+1862)/2;
lam0=0.01/lam0; kap=(lam0*alsred)/(4*pi); nkao=(epssh+kap^2)^0.5;
epsdsh=2*nkao*kap; Rs=(70+87)/200; Dis=4*((1+Rs)^2)-4*((1-Rs)^2)*(1+kap^2);
Dis=Dis^0.5; n1=(2*(1+Rs)-Dis)/(2*(1-Rs)); n2=(2*(1+Rs)+Dis)/(2*(1-Rs));
nmus=(1.552+1.542+1.527)/3; nkao=(n1+nkao+nmus)/3;
nv=(1.538+1.542)/2; sig=5.668e-8;
%Rfl=((nv-1)^2+kap^2)/((nv+1)^2+kap^2);
tol=1e-5; y0=5e-4; Ns=round(y0/tol); Tr=1;
%for k=1:Ns Tr=Tr*exp(-alsred*tol); 
%Tr=Tr*(1-Rfl); end
%alsh=(-log(Tr)/y0); 
Tr=Tr*exp(-alsred*y0); Rflk=((nv-nkao)^2)/((nv+nkao)^2); Abs=1-Tr; 
Tr=(1-Rflk)*Tr; Ab=alsred*y0; Abst=Tr;
%Abl=alsh*y0; 
Tw1=480+273.15; Tw2=470+273.15; Tz=(Tw1+Tw2)/2;
temv = [20 50 100 250 500]; temv = temv+273.15;
tepv = [0.045 0.047 0.051 0.065 0.096]; tepv=tepv*4184/3600;
kotepv = koef(temv,tepv,length(tepv))';
lamb1=polyval(kotepv,Tw1); lamb2=polyval(kotepv,Tw2); lamb=(lamb1+lamb2)/2;
%qrez=sig*(Tw1^4-Tw2^4)/(2/Abst-1+3*alsred*y0/4);
%ko=0; tem=0; Nt=100; hy=y0/Nt;
%for k=1:Nt
%    ko(k)=k*hy; tem(k)=(Tw1^4-qrez*0.75*alsred*ko(k)/sig)^0.25-273.15; ko(k)=ko(k)*1e6;
%end
y=0; tem=0; ko=0;
for k=1:Ns
    tem(k)=((Tw1^4+(Tw2^4-Tw1^4)*y/y0)^0.25); koo(k)=y; y=y+tol;
end
r1=Rflk; r2=r1;
Aef=4*Ab*(1-exp(-4*Ab)*r1*r2)-(r1+r2)*(1-exp(-2*Ab))^2-2*(1-r1*r2**exp(-2*Ab))*(1-exp(-2*Ab));
Aef=Aef/(4*(Ab^2)*(1-r1*r2*exp(-4*Ab)));
qob=(Tw1-Tw2)*lamb/y0+Aef*sig*(Tw1^4-Tw2^4);
az1=(r1*(1-exp(-2*Ab))+(1-r1*r2*exp(-2*Ab)))/(2*(1-r1*r2*exp(-4*Ab)));
az2=(r2*(1-exp(-2*Ab))+(1-r1*r2*exp(-2*Ab)))*exp(-2*Ab)/(2*(1-r1*r2*exp(-4*Ab)));
%al=alsh;
al=alsred; y=0;
for k=1:Nt
ko(k)=y;
tem(k)=Tw1+(sig*(Tw1^4-Tw2^4)/(Ab*lamb)-qob/lamb)*ko(k)+(sig*(Tw1^4-Tw2^4))*(az1*exp(-2*al*ko(k))-az2*exp(2*al*ko(k)))/(2*(al^2)*lamb*y0)-273.15;
y=y+tol;
end
p=plot(ko,tem,'-b');
set(p,'LineWidth',3);
hold on;
grid on;
xlabel({'Координата, мкм'});
ylabel({'Температура, град С'});
title({'График зависимости температуры от координаты'});