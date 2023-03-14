function [ isdan ] = danPoTemH6(temvs,temvh)
%Засыпка плоско-параллельная, после обжига при 1000 °С
tc1=0; tc2=0; tem1=0; tem2=0;
tem1=temvs(25);
tem2=temvs(26);
tc1=(temvh(25)*temvs(25))/tem1;
tc2=(temvh(26)*temvs(26))/tem2;
tem1=tem1/1; tem2=tem2/1;
k1=(tc2-tc1)/(tem2-tem1); k2=tc2-k1*tem2; 
isdan=[k1 k2];
end