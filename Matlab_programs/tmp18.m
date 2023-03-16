%поиск ЛКТП при заданных температурах
function n = tmp18(tem1,tem2,tepv1,tepv2)
tem1=1*tem1
tem2=1*tem2
alf1=1e6/sredRosSieg(tem1)
alf2=1e6/sredRosSieg(tem2)
npT1=nsreddvVer(tem1)
npT2=nsreddvVer(tem1-1);
GnT=npT1-npT2
sigma=5.67e-8;
t=abs(8*npT1*sigma*(tem1^3)/(3*alf1));
lamizl=t*(2*npT1+tem1*GnT)
ktp=tepv1
lamte=ktp-lamizl
npT1=nsreddvVer(tem2)
npT2=nsreddvVer(tem2-1);
GnT=npT1-npT2
t=abs(8*npT1*sigma*(tem2^3)/(3*alf2));
lamizl=t*(2*npT1+tem2*GnT)
ktp=tepv2
lamte=ktp-lamizl
n=0;
end