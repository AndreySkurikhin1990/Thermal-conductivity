function [ isdan ] = danPoTemH4(temvs,temvh)
%Засыпка вертикальная, после обжига при 1000 °С
tc1=0; tc2=0; tem1=0; tem2=0;
tem1=temvs(1)+temvs(3);
tem2=temvs(2)+temvs(4);
tc1=(temvh(1)*temvs(1)+temvh(3)*temvs(3))/tem1;
tc2=(temvh(2)*temvs(2)+temvh(4)*temvs(4))/tem2;
tem1=tem1/2; tem2=tem2/2;
k1=(tc2-tc1)/(tem2-tem1); k2=tc2-k1*tem2; 
isdan=[k1 k2];
end