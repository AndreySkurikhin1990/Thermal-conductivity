function [ isdan ] = danPoTemC4(temvs,temvc)
%Засыпка вертикальная, после обжига при 1000 °С
tc1=0; tc2=0; tem1=0; tem2=0;
tem1=temvs(1)+temvs(3);
tem2=temvs(2)+temvs(4);
tc1=(temvc(1)*temvs(1)+temvc(3)*temvs(3))/tem1;
tc2=(temvc(2)*temvs(2)+temvc(4)*temvs(4))/tem2;
tem1=tem1/2; tem2=tem2/2;
k1=(tc2-tc1)/(tem2-tem1); k2=tc2-k1*tem2; 
isdan=[k1 k2];
end