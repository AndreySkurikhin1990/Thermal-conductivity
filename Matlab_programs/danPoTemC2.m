function [ isdan ] = danPoTemC2(temvs,temvc)
%Засыпка вертикальная, исходный
tc1=0; tc2=0; tem1=0; tem2=0;
tem1=temvs(5)+temvs(7)+temvs(11);
tem2=temvs(6)+temvs(8)+temvs(12);
tc1=(temvc(5)*temvs(5)+temvc(7)*temvs(7)+temvc(11)*temvs(11))/tem1;
tc2=(temvc(6)*temvs(6)+temvc(8)*temvs(8)+temvc(12)*temvs(12))/tem2;
tem1=tem1/3; tem2=tem2/3;
k1=(tc2-tc1)/(tem2-tem1); k2=tc2-k1*tem2; 
isdan=[k1 k2];
end