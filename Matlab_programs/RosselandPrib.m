epssh=(9+12)/2; 
ro=2.6e3; 
vol=(5e-3)^3; 
hi=-0.5e-3; 
pi=3.1415926535897932;
mu=1+4*pi*hi*ro*vol; 
epssh=mu*epssh; 
alsred=IntegrDlVol; 
lam0=(3000+1862)/2;
lam0=0.01/lam0; 
kap=(lam0*alsred)/(4*pi); 
nkao=(epssh+kap^2)^0.5;
epsdsh=2*nkao*kap; 
Rs=(70+87)/200; 
Dis=4*((1+Rs)^2)-4*((1-Rs)^2)*(1+kap^2);
Dis=Dis^0.5; 
n1=(2*(1+Rs)-Dis)/(2*(1-Rs)); 
n2=(2*(1+Rs)+Dis)/(2*(1-Rs));
nmus=(1.552+1.542+1.527)/3; 
nkao=(n1+nkao+nmus)/3;
nv=(1.538+1.542)/2; sig=5.668e-8;
%Rfl=((nv-1)^2+kap^2)/((nv+1)^2+kap^2);
tol=1e-7;
hy=tol;
y0=1e-3; 
Ns=round(y0/tol); Tr=1;
Nt=Ns;
Nit=1;
%for k=1:Ns Tr=Tr*exp(-alsred*tol); %Tr=Tr*(1-Rfl); end %alsh=(-log(Tr)/y0); 
Tr=Tr*exp(-alsred*y0); 
Rflk=((nv-nkao)^2)/((nv+nkao)^2); Abs=1-Tr; 
Tr=(1-Rflk)*Tr; Ab=alsred*y0; Abst=Tr; %Abl=alsh*y0; 
Tw1=630+273.15; 
Tw2=610+273.15; 
temv = [20 50 100 250 500]; temv = temv+273.15;
tepv = [0.045 0.047 0.051 0.065 0.096]; tepv=tepv*4184/3600;
kotepv = koef(temv,tepv,length(tepv))';
qrez=(lamb/y0)*(Tw1-Tw2)+sig*(Tw1^4-Tw2^4)/(2/Abst-1+0.75*alsred*y0)
lmb1=qrez*y0/(Tw1-Tw2)
qrez=(lamb/y0)*(Tw1-Tw2)+sig*(Tw1^4-Tw2^4)/(2/Abst-1)
lmb2=qrez*y0/(Tw1-Tw2)
tem=0; hy=y0/Nt; ko=0; q=1; lambd=0; fl=0;
for j=1:Nit
y=0;
for k=1:Ns+1
    ko(k)=y; 
%    tem(k)=(Tw1^4-qrez*0.75*alsred*y/sig)^0.25; 
    ko(k)=ko(k)*1e6; y=y+hy;
end
%for g=1:length(tem)
%if (abs(imag(tem(g)))>0) || (tem(g)<0) fl=1; break; end
%end
%if (fl==1) 
%    break; 
%end
%Tw2=tem(length(tem));
%Tw1=tem(1);
Tz=(Tw1+Tw2)/2;
lamb=polyval(kotepv,Tz);
qrez=(lamb/y0)*(Tw1-Tw2)+(4*sig)*(Tw1^4-Tw2^4)/(3*alsred*y0)
%if (abs(imag(qrez))>0) break; end
dely=(ko(length(ko))-ko(1))*(1e-6)
lambd(q)=qrez*dely/(Tw1-Tw2)
q=q+1;
end
%j=length(ko);
%dely=(ko(j)-ko(1))*(1e-6)
disp(Tw1-Tw2);
%disp(lambd(q-1));
%disp(qrez*dely/(Tw1-Tw2));
%p=plot(ko,tem-273.15,'-b');
%disp(trapz(lambd)/(length(lambd)-1));
%p=plot(1:length(lambd),lambd,'-b');
%set(p,'LineWidth',3);
%hold on;
%grid on;
%xlabel({'Координата, мкм'});
%ylabel({'Относительная температура'});
%title({'График зависимости температуры от координаты'});