format long g; y0=3e4*1e-6; %temvh=[558 950 560 950 545 927 530 925 600 1000 587 1000 590 1000 540 920 540 920 580 1000 580 1000 560 900]; temvh=temvh+te0;%temvc=[149 350 149 357 168 369 168 369 107 242 111 240 109 233 166 375 171 384 123 294 127 291 184 396]; temvc=temvc+te0;%tepo=[4.7766 11.5016 3.55732 11.4997 3.9805 9.5532 3.6447 9.4779 3.5144 8.6593 3.352 9.218 3.0313 7.7946 3.4023 10.2068 3.92812 11.17333 2.977 8.0448 3.0671 6.1342 4.4624 11.6021]; 
%alfs=0; alfs=RasMasKoAbs(); ep=1e-20; npp=0; npp=Kramers_n(); 
te0=273.15; slu=3; nac=slu; kon=slu;
temvh=arrTemHigh(); temvh=temvh+te0; 
temvc=arrTemCold(); temvc=temvc+te0; 
tepo=arrTepPot84(); qv=(1e4/13.85)*tepo; pt=length(temvc); ts=sredtemp(temvh,temvc,pt); 
dl=RasshDiapDlinVoln(); tepv=koeftepv(qv,y0,temvh,temvc,pt);
temvs=(1/2)*(temvc+temvh); %tepv1=0; tepv2=0; tem1=0; tem2=0; for k=1:p if (rem(k,2)==1) tem1=tem1+temvs(k); tepv1=tepv1+tepv(k); else tem2=tem2+temvs(k); tepv2=tepv2+tepv(k); end; end; tem1=tem1/(p/2);tem2=tem2/(p/2); tepv1=tepv1/(p/2); tepv2=tepv2/(p/2); 
for slu=nac:kon
switch (slu)
    case 1
        Tz=[353+352,621+647]/2+te0
        km=danPoTemTepl(temvs,tepv); kt1=km(1); kt2=km(2); %Засыпка плоско-параллельная, исходный
        kc=danPoTemC(temvs,temvc); ktc1=kc(1); ktc2=kc(2);
        kh=danPoTemH(temvs,temvh); kth1=kh(1); kth2=kh(2);
    case 2
        Tz=[372+354+353,648+653+648]/3+te0
        km=danPoTemTepl2(temvs,tepv); kt1=km(1); kt2=km(2); %Засыпка вертикальная, исходный
        kc=danPoTemC2(temvs,temvc); ktc1=kc(1); ktc2=kc(2);
        kh=danPoTemH2(temvs,temvh); kth1=kh(1); kth2=kh(2);
    case 3
        Tz=[624,900];
        km=danPoTemTepl3(temvs,tepv); kt1=km(1); kt2=km(2); %Засыпка плоско-параллельная, повторные измерения
        kc=danPoTemC3(temvs,temvc); ktc1=kc(1); ktc2=kc(2);
        kh=danPoTemH3(temvs,temvh); kth1=kh(1); kth2=kh(2);
    case 4
        Tz=[626,921];
        km=danPoTemTepl4(temvs,tepv); kt1=km(1); kt2=km(2); %Засыпка вертикальная, после обжига при 1000 °С
        kc=danPoTemC4(temvs,temvc); ktc1=kc(1); ktc2=kc(2);
        kh=danPoTemH4(temvs,temvh); kth1=kh(1); kth2=kh(2);
    case 5
        Tz=[628,924];
        km=danPoTemTepl5(temvs,tepv); kt1=km(1); kt2=km(2); %Засыпка вертикальная, повторы
        kc=danPoTemC5(temvs,temvc); ktc1=kc(1); ktc2=kc(2);
        kh=danPoTemH5(temvs,temvh); kth1=kh(1); kth2=kh(2);
    case 6
        Tz=[625,920];
        km=danPoTemTepl6(temvs,tepv); kt1=km(1); kt2=km(2); %Засыпка плоско-параллельная, после обжига при 1000 °С
        kc=danPoTemC6(temvs,temvc); ktc1=kc(1); ktc2=kc(2);
        kh=danPoTemH6(temvs,temvh); kth1=kh(1); kth2=kh(2);    
end
%qvm=danPoPlotTeplPot(temvs,qv); qk1=qvm(1);qk2=qvm(2);%tem1=temvc(1); tem2=temvh(2); for k=1:p     if (tem2<temvh(k)) tem2=temvh(k);     end;    if (tem1>temvc(k))        tem1=temvc(k);    end; end;
%Tz= 1e2:5e1:12e2; Tz=Tz+te0;
p=length(Tz); 
%koeftep=0; ptpo=0; teco=0; teho=0;GrT=0; npT=0;
for k=1:p
    %teco(k)=polyval(kc,Tz(k));
    %teho(k)=polyval(kh,Tz(k));
    koeftep(k)=abs(polyval(km,Tz(k)));
    %ptpo(k)=koeftep(k)*abs(teho(k)-teco(k))/y0;
    %GrT(k)=ptpo(k)/koeftep(k);
    %npT(k)=1;
    %npT(k)=nsreddvVer(Tz(k));
    %alf(k)=1e-3;
    %alf(k)=1e6/sredRosSieg(Tz(k),npp,alfs,dl);
    %disp(k);
end
%GnpT=0; %npT=1*npT;
%for k=2:p
    %GnpT(k)=(npT(k)-npT(k-1))/(Tz(k)-Tz(k-1));
%end
%sigma=5.67e-8; lamizl=0; lamte=0;
%for k=2:p
    %t=abs(8*npT(k)*sigma*(Tz(k)^3)/(3*alf(k)));
    %t=t*(2*npT(k)+Tz(k)*GnpT(k));
    %lamizl(k)=t;
    %lamte(k)=koeftep(k)-t;
%end
slu=slu'
Tz=Tz'
%lamizl=lamizl'
koeftep=koeftep'
end
%p=plot(Tz(2:p),lamizl(2:p),'-b',Tz(2:p),lamte(2:p),'-k',Tz(2:p),koeftep(2:p),'-c');
%set(p,'LineWidth',3); hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Составляющие коэффициента теплопроводности, Вт/(м*К)'}); 
%title({'График зависимости составляющих коэффициента теплопроводности от температуры'});
%legend('Излучательная компонента', 'Кондуктивная часть КТП', 'Общий КТП','location','best');