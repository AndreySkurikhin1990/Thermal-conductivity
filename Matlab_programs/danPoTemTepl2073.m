function [ isdan ] = danPoTemTepl2073(temvs,tepv)
%Фракция 2-0,7 мм, после обжига при 1000 °С
tepv1=0; tepv2=0; tem1=0; tem2=0;
tem1=temvs(7)+temvs(9);
tem2=temvs(8)+temvs(10);
tepv1=(tepv(7)*temvs(7)+tepv(9)*temvs(9))/tem1;
tepv2=(tepv(8)*temvs(8)+tepv(10)*temvs(10))/tem2;
tem1=tem1/2; tem2=tem2/2;
k1=(tepv2-tepv1)/(tem2-tem1); k2=tepv2-k1*tem2; 
isdan=[k1 k2];
end