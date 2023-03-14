function [ isdan ] = danPoTemC2072(temvs,temvc)
%Фракция 2-0,7 мм (повторные измерения)
tc1=0; tc2=0; tem1=0; tem2=0;
tem1=temvs(3)+temvs(5);
tem2=temvs(4)+temvs(6);
tc1=(temvc(3)*temvs(3)+temvc(5)*temvs(5))/tem1;
tc2=(temvc(4)*temvs(4)+temvc(6)*temvs(6))/tem2;
tem1=tem1/2; tem2=tem2/2;
k1=(tc2-tc1)/(tem2-tem1); k2=tc2-k1*tem2; 
isdan=[k1 k2];
end