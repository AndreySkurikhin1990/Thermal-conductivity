%Определение КТП ТК шамота, ИТОМ, КВИ и вермикулита
function [ t ] = tmp55()
format longg; 
vybves=3; 
te0=273.15; tem = 2e2:1e2:12e2; tem=tem+te0;
y0=3e1*1e-3; %vmivmf=1; vfv=1; vyuv=2; vysv=2;
vysh=[2]; 
%vybmaritom=[0, 1, 2, 3];
%vybmarkvi=[4, 5, 6, 7, 8, 9, 10];
vybmaritom=[0];
vmkvi=[10];
%nvyfv=[0, 1]; %фракции вермикулита
%nvysv=[0, 1, 2]; %состояния вермикулита
%nvmivmf=[0, 1, 2]; %стационарные методы измерений - 2019 и 2020
%nvyuv=[1, 2]; %укладка вермикулита
nvyfv=[0]; %фракции вермикулита
nvysv=[2]; %состояния вермикулита
nvmivmf=[2]; %стационарные методы измерений - 2019 и 2020
nvyuv=[2]; %укладка вермикулита
b=length(tem); eps=1e-8;
for k=1:b
    s(k)=0;
    c(k)=0;
end
%-----
%-----
switch (vybves)
    case (1)
for kvfv=1:length(nvyfv)
for kvysv=1:length(nvysv)
for kvyuv=1:length(nvyuv)
for kvmivmf=1:length(nvmivmf)
vyfv=nvyfv(kvfv);
vysv=nvysv(kvysv);
vmivmf=nvmivmf(kvmivmf);
vyuv=nvyuv(kvyuv);
ktptkv=opredKTPTverKarkVerm(vyfv, tem, vysv, vyuv, vmivmf, y0);
for k=1:length(ktptkv)
    s(k)=s(k)+ktptkv(k);
    if (ktptkv(k)>eps)
    c(k)=c(k)+1;
    end
end
if (vyfv==1)
    break;
end
if (vysv>0)
    break;
end
end
if (vyfv==0)
    break;
end
end
end
end
%-----Шамот
    case (0)
        for vysha=1:length(vysh)
            vysham=vysh(vysha)
    ktptks=opredKTPTverKarkSha(tem, vysham)
    %for k=1:length(ktptks)
    %s(k)=s(k)+ktptks(k);
    %if (ktptks(k)>eps)
    %c(k)=c(k)+1;
    %end
    %end
        end
%-----ИТОМ
    case (2)
        disp('ITOM');
            for k=1:length(vybmaritom)
                vm=vybmaritom(k)
                ktptkitom=opredKTPTverKarkitom(vm, tem, y0)
            end
    case (3)
%for vmkvi=4:10
vmkvi=vmkvi'
ktptkk=opredKTPTverKarkkvi(vmkvi, tem);
%end
%t=ModelirTeploperItom(ktptks,ktptkv,tem);
%vmitom=0; t=RazPorItom(vmitom); 
end
c=c;
for k=1:length(s)
    if (c(k)>eps)
        s(k)=s(k)/c(k);
    end
end
t=s;
end