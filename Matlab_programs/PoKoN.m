function kor = PoKoN(sig, ro, teta1, teta2, alp, m, d, lam)
dt=0;
ta=teta1-dt; %tetad
tb=teta2+dt; %tetad
dete=1e-2;
eps=1e3;
h=1;
hk=1e5;
e=1e-4;
while ((eps>e) && (h<hk))
tc=(ta+tb)/2;
ksi=(1+ro)/(1-ro);
ksi1=(1-exp(-m*d));
ksi1=ksi1/(1+exp(-m*d));
teta14=teta1^4;
teta24=teta2^4;
sig1=sig*(teta14-teta24);
sig2=alp*d/2;
ksi2=sig1/(ksi+sig2); %H11
ksi5=(alp*lam/2);
ksi3=ksi5/(ksi+sig2); %H12
ksi4=ksi1*alp/m; %H22
teta04a=(teta1^4+teta2^4-ta^4);
teta04b=(teta1^4+teta2^4-tb^4);
teta04c=(teta1^4+teta2^4-tc^4);
teta0a=(teta04a)^(1/4);
teta0b=(teta04b)^(1/4);
teta0c=(teta04c)^(1/4);
Ha=(ksi2+ksi3*(teta0a-ta));
Hb=(ksi2+ksi3*(teta0b-tb));
Hc=(ksi2+ksi3*(teta0c-tc));
Aa=Ha*ksi4+2*sig*teta04a+2*ksi5*teta0a;
Ab=Hb*ksi4+2*sig*teta04b+2*ksi5*teta0b;
Ac=Hc*ksi4+2*sig*teta04c+2*ksi5*teta0c;
ta4=ta^4;
tb4=tb^4;
tc4=tc^4;
fa=2*sig*(ta4)+2*sig2*Ha+2*ksi5*ta-ksi4*Ha-Aa;
fb=2*sig*(tb4)+2*sig2*Hb+2*ksi5*tb-ksi4*Hb-Ab;
fc=2*sig*(tc4)+2*sig2*Hc+2*ksi5*tc-ksi4*Hc-Ac;
if (((fa*fc)<0) && ((fb*fc)>0))
    tb=tc; 
elseif (((fa*fc)>0) && ((fb*fc)<0))
    ta=tc;
%elseif (((fa*fc)>0) && ((fb*fc)>0))
    %ta=ta+dete;
    %tb=tb-dete;
end
eps=abs(ta-tb);
h=h+1;
end
kor=teta0c;
end