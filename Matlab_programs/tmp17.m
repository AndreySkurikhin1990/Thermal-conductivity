%Фракция 8-4 мм
y0=3e4*1e-6;
te0=273.15;
temvh=arrTemHigh(); 
%temvh=arrTemHigh207();
temvh=temvh+te0; 
temvc=arrTemCold(); 
%temvc=arrTemCold207(); 
temvc=temvc+te0; 
temvs=(temvh+temvc)/2;
tepo=arrTepPot84(); 
%tepo=arrTepPot207();
qv=(1/13.85e-4)*tepo; 
tepv=koeftepv(qv,y0,temvh,temvc,length(temvc));
%---------------------------
%Засыпка плоско-параллельная, исходный, фракция 8-4 мм
tem1=temvs(15)+temvs(17);
tem2=temvs(16)+temvs(18);
tepv1=(tepv(15)*temvs(15)+tepv(17)*temvs(17))/tem1;
tepv2=(tepv(16)*temvs(16)+tepv(18)*temvs(18))/tem2;
tem1=tem1/2; 
tem2=tem2/2;
%tem1=tmp18(tem1,tem2,tepv1,tepv2);
%---------------------------
%Засыпка вертикальная, исходный, фракция 8-4 мм
tepv1=0; tepv2=0; tem1=0; tem2=0;
tem1=temvs(5)+temvs(7)+temvs(11);
tem2=temvs(6)+temvs(8)+temvs(12);
tepv1=(tepv(5)*temvs(5)+tepv(7)*temvs(7)+tepv(11)*temvs(11))/tem1;
tepv2=(tepv(6)*temvs(6)+tepv(8)*temvs(8)+tepv(12)*temvs(12))/tem2;
tem1=tem1/3; tem2=tem2/3;
%tem1=tmp18(tem1,tem2,tepv1,tepv2);
%---------------------------
%Засыпка плоско-параллельная, повторные измерения, фракция 8-4 мм
tepv1=0; tepv2=0; tem1=0; tem2=0;
tem1=temvs(19)+temvs(21)+temvs(23);
tem2=temvs(20)+temvs(22)+temvs(24);
tepv1=(tepv(19)*temvs(19)+tepv(21)*temvs(21)+tepv(23)*temvs(23))/tem1;
tepv2=(tepv(20)*temvs(20)+tepv(22)*temvs(22)+tepv(24)*temvs(24))/tem2;
tem1=tem1/3; tem2=tem2/3;
%tem1=tmp18(tem1,tem2,tepv1,tepv2);
%---------------------------
%Засыпка вертикальная, после обжига при 1000 °С, фракция 8-4 мм
tepv1=0; tepv2=0; tem1=0; tem2=0;
tem1=temvs(1)+temvs(3);
tem2=temvs(2)+temvs(4);
tepv1=(tepv(1)*temvs(1)+tepv(3)*temvs(3))/tem1;
tepv2=(tepv(2)*temvs(2)+tepv(4)*temvs(4))/tem2;
tem1=tem1/2; tem2=tem2/2;
%tem1=tmp18(tem1,tem2,tepv1,tepv2);
%---------------------------
%Засыпка вертикальная, повторы, фракция 8-4 мм
tepv1=0; tepv2=0; tem1=0; tem2=0;
tem1=temvs(9)+temvs(13);
tem2=temvs(10)+temvs(14);
tepv1=(tepv(9)*temvs(9)+tepv(13)*temvs(13))/tem1;
tepv2=(tepv(10)*temvs(10)+tepv(14)*temvs(14))/tem2;
tem1=tem1/2; tem2=tem2/2;
%tem1=tmp18(tem1,tem2,tepv1,tepv2);
%---------------------------
%Засыпка плоско-параллельная, после обжига при 1000 °С, фракция 8-4 мм
tepv1=0; tepv2=0; tem1=0; tem2=0;
tem1=temvs(25);
tem2=temvs(26);
tepv1=(tepv(25)*temvs(25))/tem1;
tepv2=(tepv(26)*temvs(26))/tem2;
%tem1=tmp18(tem1,tem2,tepv1,tepv2);
%---------------------------
%Засыпка хаотичная, исходный, фракция 2-0,7 мм
%tem1=temvs(1);
%tem2=temvs(2);
%tepv1=(temvh(1)*temvs(1))/tem1;
%tepv2=(temvh(2)*temvs(2))/tem2;
%tem1=tem1/1; tem2=tem2/1;
%tem1=tmp18(tem1,tem2,tepv1,tepv2);
%---------------------------
%Фракция 2-0,7 мм (повторные измерения), фракция 2-0,7 мм
%tem1=temvs(3)+temvs(5);
%tem2=temvs(4)+temvs(6);
%tepv1=(tepv(3)*temvs(3)+tepv(5)*temvs(5))/tem1;
%tepv2=(tepv(4)*temvs(4)+tepv(6)*temvs(6))/tem2;
%tem1=tem1/2; tem2=tem2/2;
%tem1=tmp18(tem1,tem2,tepv1,tepv2);
%---------------------------
%Фракция 2-0,7 мм, после обжига при 1000 °С, фракция 2-0,7 мм
%tem1=temvs(7)+temvs(9);
%tem2=temvs(8)+temvs(10);
%tepv1=(tepv(7)*temvs(7)+tepv(9)*temvs(9))/tem1;
%tepv2=(tepv(8)*temvs(8)+tepv(10)*temvs(10))/tem2;
%tem1=tem1/2; tem2=tem2/2;
%tem1=tmp18(tem1,tem2,tepv1,tepv2);