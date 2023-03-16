function eps = tmp34()
te=[626
921
633
923
628
924
626
907
624
900
625
920];
npp=Kramers_n();
for k=1:length(npp)
    epsi(k) = epsilnu(npp(k));
end
dv = RasshDiapDlinVoln();
for k=1:length(te)
ep(k) = epssred(dv,epsi,te(k),npp);
end
eps=ep';
end