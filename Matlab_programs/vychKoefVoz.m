function magn = vychKoefVoz(ko,T)
te0=273.15;
tem=arrTempAir();
tem=tem+te0;
p=length(tem);
t1=1; f=1;
for k=1:p
        if ((tem(k)<T) && (f>0))
            t1=t1+1;
            f=0;
        else
            break;
        end
end
    if ((f>0) || (t1==p))
        t1=p-1;
    end
    if ((t1==1) && (f==0))
        t1=2;
    end
    k=t1;
    T1=tem(k-1);
    T2=tem(k);
    temvo=[T1,T2];
if (ko==1)
    Pr=PrAir();
    Pr1=Pr(k-1);
    Pr2=Pr(k);
    arPr=[Pr1,Pr2];
    arKoPr=polyfit(tem,Pr,length(arPr))';
    Pran=polyval(arKoPr,T);
    magn=Pran;
end
    if (ko==2)
        rovoz=densAir();
    ro1=rovoz(k-1);
    ro2=rovoz(k);
    rovo=[ro1,ro2];
    arKoPlotVoz=polyfit(tem,rovoz,length(rovo))';
    PlotVoz=polyval(arKoPlotVoz,T);
    magn=PlotVoz;
    end
if (ko==3)
        muvoz=koefDinVyazAir();
    mu1=muvoz(k-1);
    mu2=muvoz(k);
    muvo=[mu1,mu2];
    arKoDinVyazVoz=polyfit(tem,muvoz,length(muvo))';
    koDinVisAir=polyval(arKoDinVyazVoz,T);
    magn=koDinVisAir;
end
if (ko==4)
    lam=koefTeploprovAir();
    lam1=lam(k-1);
    lam2=lam(k);
    arlam=[lam1,lam2];
    arKoTep=polyfit(tem,lam,length(arlam))';
    koTepVoz=polyval(arKoTep,T);
    magn=koTepVoz;
end
if (ko==5)
        cp=teploemPostDavAir();
        cp1=cp(k-1);
        cp2=cp(k);
        cpvo=[cp1,cp2];
        arUdTeplPostDavl=polyfit(tem,cp,length(cpvo))';
        UdTepPostDavl=polyval(arUdTeplPostDavl,T);
        magn=UdTepPostDavl;
end
end