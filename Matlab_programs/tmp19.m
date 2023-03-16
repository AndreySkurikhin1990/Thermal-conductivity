%Фракция 8-4 мм
%Фракция 2-0,7 мм
y0=3e4*1e-6;
te0=273.15;
%temvh=arrTemHigh(); 
temvh=arrTemHigh207();
temvh=temvh+te0; 
%temvc=arrTemCold(); 
temvc=arrTemCold207(); 
temvc=temvc+te0; 
temvs=(temvh+temvc)/2;
%tepo=arrTepPot84(); 
tepo=arrTepPot207();
qv=(1e4/13.85)*tepo; 
tepv=koeftepv(qv,y0,temvh,temvc,length(temvc));
%---------------------------
%Засыпка плоско-параллельная, исходный
%tem1=temvs(15)+temvs(17);
%tem2=temvs(16)+temvs(18);
%tepv1=(tepv(15)*temvs(15)+tepv(17)*temvs(17))/tem1;
%tepv2=(tepv(16)*temvs(16)+tepv(18)*temvs(18))/tem2;
%tem1=tem1/2; 
%tem2=tem2/2;
%---------------------------
%Засыпка хаотичная, исходный, фракция 2-0,7 мм
%tem1=temvs(1);
%tem2=temvs(2);
%tepv1=(temvh(1)*temvs(1))/tem1;
%tepv2=(temvh(2)*temvs(2))/tem2;
%tem1=tem1/1
%tem2=tem2/1
%---------------------------
%Фракция 2-0,7 мм (повторные измерения)
%tem1=temvs(3)+temvs(5);
%tem2=temvs(4)+temvs(6);
%tepv1=(tepv(3)*temvs(3)+tepv(5)*temvs(5))/tem1;
%tepv2=(tepv(4)*temvs(4)+tepv(6)*temvs(6))/tem2;
%tem1=tem1/2
%tem2=tem2/2
%---------------------------
%Фракция 2-0,7 мм, после обжига при 1000 °С
tem1=temvs(7)+temvs(9);
tem2=temvs(8)+temvs(10);
tepv1=(tepv(7)*temvs(7)+tepv(9)*temvs(9))/tem1;
tepv2=(tepv(8)*temvs(8)+tepv(10)*temvs(10))/tem2;
tem1=tem1/2
tem2=tem2/2
np=Kramers_n();
dl=RasshDiapDlinVoln();
eps=0;
for k=1:length(dl)
eps(k)=epsilnu(n(k));
end
es=epssred(dl,eps,tem1,np)
es=epssred(dl,eps,tem2,np)