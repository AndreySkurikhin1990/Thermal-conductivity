function [ isdan ] = danPoTemH5(temvs,temvh)
%Засыпка вертикальная, повторы
tc1=0; tc2=0; tem1=0; tem2=0;
tem1=temvs(9)+temvs(13);
tem2=temvs(10)+temvs(14);
tc1=(temvh(9)*temvs(9)+temvh(13)*temvs(13))/tem1;
tc2=(temvh(10)*temvs(10)+temvh(14)*temvs(14))/tem2;
tem1=tem1/2; tem2=tem2/2;
k1=(tc2-tc1)/(tem2-tem1); k2=tc2-k1*tem2; 
isdan=[k1 k2];
end