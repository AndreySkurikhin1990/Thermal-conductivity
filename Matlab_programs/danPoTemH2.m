function [ isdan ] = danPoTemH2(temvs,temvh)
%Засыпка вертикальная, исходный
tc1=0; tc2=0; tem1=0; tem2=0;
tem1=temvs(5)+temvs(7)+temvs(11);
tem2=temvs(6)+temvs(8)+temvs(12);
tc1=(temvh(5)*temvs(5)+temvh(7)*temvs(7)+temvh(11)*temvs(11))/tem1;
tc2=(temvh(6)*temvs(6)+temvh(8)*temvs(8)+temvh(12)*temvs(12))/tem2;
tem1=tem1/3; tem2=tem2/3;
k1=(tc2-tc1)/(tem2-tem1); k2=tc2-k1*tem2; 
isdan=[k1 k2];
end