function kor = PoiKor(sig, Aef, lamvo, tn, dlv, jv)
ta=tn-50;
tb=tn+50;
eps=1;
h=1;
while (eps>1e-7)
tc=(ta+tb)/2;    
fc=sig*Aef*(tn^4-tc^4)+lamvo*(tn-tc)/dlv-jv;
fa=sig*Aef*(tn^4-ta^4)+lamvo*(tn-ta)/dlv-jv;
fb=sig*Aef*(tn^4-tb^4)+lamvo*(tn-tb)/dlv-jv;
if ((fa*fc)<0)
elseif ((fb*fc)>0)
    tb=tc; 
end
if ((fa*fc)>0) 
elseif ((fb*fc)<0) 
    ta=tc; 
end
eps=abs(fa-fb);
h=h+1;
if (h>1e2) 
    break; end
end
kor=tc;
end