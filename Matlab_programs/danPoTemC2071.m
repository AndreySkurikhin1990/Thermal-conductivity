function [ isdan ] = danPoTemC2071(temvs,temvc)
%Исходный, фракция 2-0,7 мм
tc1=0; tc2=0; tem1=0; tem2=0;
tem1=temvs(1);
tem2=temvs(2);
tc1=(temvc(1)*temvs(1))/tem1;
tc2=(temvc(2)*temvs(2))/tem2;
tem1=tem1/1; tem2=tem2/1;
k1=(tc2-tc1)/(tem2-tem1); k2=tc2-k1*tem2; 
isdan=[k1 k2];
end