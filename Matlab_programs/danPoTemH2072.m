function [ isdan ] = danPoTemH2072(temvs,temvh)
%Фракция 2-0,7 мм (повторные измерения)
tc1=0; tc2=0; tem1=0; tem2=0;
tem1=temvs(3)+temvs(5);
tem2=temvs(4)+temvs(6);
tc1=(temvh(3)*temvs(3)+temvh(5)*temvs(5))/tem1;
tc2=(temvh(4)*temvs(4)+temvh(6)*temvs(6))/tem2;
tem1=tem1/2; tem2=tem2/2;
k1=(tc2-tc1)/(tem2-tem1); k2=tc2-k1*tem2; 
isdan=[k1 k2];
end