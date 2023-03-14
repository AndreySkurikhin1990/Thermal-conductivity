function [ isdan ] = danPoTemC(temvs,temvc)
%Засыпка плоско-параллельная, исходный
tc1=0; tc2=0; tem1=0; tem2=0;
tem1=temvs(15)+temvs(17);
tem2=temvs(16)+temvs(18);
tc1=(temvc(15)*temvs(15)+temvc(17)*temvs(17))/tem1;
tc2=(temvc(16)*temvs(16)+temvc(18)*temvs(18))/tem2;
tem1=tem1/2; tem2=tem2/2;
k1=(tc2-tc1)/(tem2-tem1); k2=tc2-k1*tem2; 
isdan=[k1 k2];
end