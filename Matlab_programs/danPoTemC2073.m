function [ isdan ] = danPoTemC2073(temvs,temvc)
%Фракция 2-0,7 мм, после обжига при 1000 °С
tc1=0; tc2=0; tem1=0; tem2=0;
tem1=temvs(7)+temvs(9);
tem2=temvs(8)+temvs(10);
tc1=(temvc(7)*temvs(7)+temvc(9)*temvs(9))/tem1;
tc2=(temvc(8)*temvs(8)+temvc(10)*temvs(10))/tem2;
tem1=tem1/2; tem2=tem2/2;
k1=(tc2-tc1)/(tem2-tem1); k2=tc2-k1*tem2; 
isdan=[k1 k2];
end