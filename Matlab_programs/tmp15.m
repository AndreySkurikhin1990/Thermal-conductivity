%определение отражения при заданном массиве температур
te=1e2:1e-1:1e3;
te0=273.15;
te=te+te0;
np=Kramers_n();
dl=RasshDiapDlinVoln();
eps=0;
ronu=0;
dv=0;
p=length(dl);
for k=1:p
eps(k)=epsilnu(np(k));
ronu(k)=dependRefln(np(k));
dv(k)=dl(k)/np(k);
end
es=0;
ro=0;
p=length(te);
for k=1:p
es(k)=epssred(dv,eps,te(k),np);
ro(k)=Refsred(ronu,np,dv,te(k));
end
te=1*te'
%es=1*es'
ro=1*ro'
%disp(Refsred(ronu,np,dl,626));
%disp(Refsred(ronu,np,dl,642));
%disp(Refsred(ronu,np,dl,630));
%disp(Refsred(ronu,np,dl,773));
%disp(Refsred(ronu,np,dl,910));
%disp(Refsred(ronu,np,dl,904));