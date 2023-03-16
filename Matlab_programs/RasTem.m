function delT = RasTem(ns, rat, kotepvo, jv, t2, l, lv, eps, sig);
dlv=rat*lv/((rat+1)*ns);
dlg=(lv-dlv)/ns;
epsil=eps/sig;
epu=1e-15;
tepf=0.24;
Aef=2/epsil-1;
for k=1:ns t(k)=t2; end
nc=1;
c=2;
for k=2:2:ns
t(c)=t(nc)-jv*dlg/tepf;
lamvo=polyval(kotepvo,t(c));
nc=nc+2;
t(nc)=PoiKor(sig, Aef, lamvo, t(c), dlv, jv);
c=c+2;
end
%disp(t(nc));
delT=-t(ns)+t(1);
end