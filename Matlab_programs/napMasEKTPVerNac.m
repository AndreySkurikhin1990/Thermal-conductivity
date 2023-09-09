function [ npv ] = napMasEKTPVerNac(vyfv, vysv, vybvykhmas, ys0, vyuv, vmivmf, ete) %isrp=0; %выбор метода измерений для фракции 2-0,7 мм: 0 - нестационарный, 1 - стационарный %выбор фракции: 0 - фракция 2-0,7 мм, 1 - фракция 8-4 мм, 2 - фракция 1,6-0,35 мм %выбор укладки: 1 - плоскопараллельная, 2 - вертикальная %выбор состояния: 0 - исходное, 1 - после повторных измерений, 2 - после прокаливания при 1000 град °С
format longg; vybvesch=1; %cem=length(ete); dete=1e2; te0=273.15; tnac=2e2+te0; cem=11; ete=tnac:dete:tnac+(cem-1)*dete;
    ident=0; thol=napMasEKTP(vyfv, vysv, vmivmf, ete, ident, ys0, vyuv);
    ident=ident+1; tgor=napMasEKTP(vyfv, vysv, vmivmf, ete, ident, ys0, vyuv);
    ident=ident+1; qob=napMasEKTP(vyfv, vysv, vmivmf, ete, ident, ys0, vyuv);
    ident=ident+1; ektpv=napMasEKTP(vyfv, vysv, vmivmf, ete, ident, ys0, vyuv);
    ident=ident+1; ete=napMasEKTP(vyfv, vysv, vmivmf, ete, ident, ys0, vyuv);
    thol=thol;
    ete=ete;
    f=1;
    k=0; tholn=vydelPol(thol, tgor, qob, ektpv, ete, k, f, vybvesch);
    k=k+1; tgorn=vydelPol(thol, tgor, qob, ektpv, ete, k, f, vybvesch);
    k=k+1; qobn=vydelPol(thol, tgor, qob, ektpv, ete, k, f, vybvesch);
    k=k+1; ektpvn=vydelPol(thol, tgor, qob, ektpv, ete, k, f, vybvesch);
    k=k+1; eten=vydelPol(thol, tgor, qob, ektpv, ete, k, f, vybvesch);
    eten=eten;
    tholn=tholn;
    switch (vybvykhmas)
        case (0)
            npv=tholn;
        case (1)
            npv=tgorn;
        case (2)
            npv=qobn;
        case (3)
            npv=ektpvn;
        case (4)
            npv=eten;
        otherwise
            npv=[0];
    end
end
           
function [ vyma ] = PoiskZavVelTem(v, efte, h, f)
	hf=1e0; e=1e-6; vybvesch=1;
    n=length(efte);
    nvyfv=[0,1]; nnvyfv=length(nvyfv); %фракции вермикулита
	nvysv=[0,1,2]; nnvysv=length(nvysv); %состояния вермикулита
	nvmivmf=[1,2]; nnvmivmf=length(nvmivmf); %стационарные методы измерений - 2019 и 2020
	nvyuv=[1,2]; nnvyuv=length(nvyuv); %укладка вермикулита
    m=length(nvyfv)*length(nvysv)*(length(nvmivmf)+length(nvyuv));
temvs=zeros(m,n);
temhq=zeros(m,n);
temcq=zeros(m,n);
ktpq=zeros(m,n);
tepv=zeros(m,n);
    nomer=1;
    vyfrve=0;
    vyukve=1;
    vysove=0;
	for kvf=1:nnvyfv
		vyfrve=nvyfv(kvf);
		for jvsv=1:nnvysv
			vysove=nvysv(jvsv);
			if ((vyfrve==0) || (vyfrve==2))
				for qvmi=1:nnvmivmf
					vymivmf=nvmivmf(qvmi); 
                    ident=0; temvct=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve);
                    ident=ident+1; temvht=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve);
                    ident=ident+1; tepvt=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve);
                    ident=ident+1; ktp=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve);
                    ident=ident+1; ts=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve);
                    fl=-1;
                    k=0; temvctn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, fl, vybvesch);
                    k=k+1; temvhtn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, fl, vybvesch);
                    k=k+1; tepvtn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, fl, vybvesch);
                    k=k+1; ktpn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, fl, vybvesch);
                    k=k+1; tsn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, fl, vybvesch);
                    if (f>0)
                    fl=1;
                    end
                    ident=0; temcq=rasMas(efte, temvhtn, temvctn, ktpn, tepvtn, tsn, temvs, temhq, temcq, ktpq, tepv, nomer, vyfrve, vysove, vymivmf, vyukve, fl, ident);
                    fl=-1;
                    ident=ident+1; temhq=rasMas(efte, temvhtn, temvctn, ktpn, tepvtn, tsn, temvs, temhq, temcq, ktpq, tepv, nomer, vyfrve, vysove, vymivmf, vyukve, fl, ident);
                    ident=ident+1; tepv=rasMas(efte, temvhtn, temvctn, ktpn, tepvtn, tsn, temvs, temhq, temcq, ktpq, tepv, nomer, vyfrve, vysove, vymivmf, vyukve, fl, ident);
                    ident=ident+1; ktpq=rasMas(efte, temvhtn, temvctn, ktpn, tepvtn, tsn, temvs, temhq, temcq, ktpq, tepv, nomer, vyfrve, vysove, vymivmf, vyukve, fl, ident);
                    ident=ident+1; temvs=rasMas(efte, temvhtn, temvctn, ktpn, tepvtn, tsn, temvs, temhq, temcq, ktpq, tepv, nomer, vyfrve, vysove, vymivmf, vyukve, fl, ident);
                    nomer=nomer+1;
                    if (vysove>0) 
                        break;
                    end
                end
            elseif (vyfrve==1)
				for qvuv=1:nnvyuv
					vyukve=nvyuv(qvuv);
                    ident=0; temvct=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve); 
                    ident=ident+1; temvht=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve); 
                    ident=ident+1; tepvt=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve); 
                    ident=ident+1; ktp=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve); 
                    ident=ident+1; ts=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve); 
                    j=-1;
                    k=0; temvctn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, j, vybvesch); 
                    k=k+1; temvhtn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, j, vybvesch); 
                    k=k+1; tepvtn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, j, vybvesch); 
                    k=k+1; ktpn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, j, vybvesch); 
                    k=k+1; tsn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, j, vybvesch); 
                    if (f>0)
                    fl=1;
                    end
                    ident=0; temcq=rasMas(efte, temvhtn, temvctn, ktpn, tepvtn, tsn, temvs, temhq, temcq, ktpq, tepv, nomer, vyfrve, vysove, vymivmf, vyukve, fl, ident);
                    fl=-1;
                    ident=ident+1; temhq=rasMas(efte, temvhtn, temvctn, ktpn, tepvtn, tsn, temvs, temhq, temcq, ktpq, tepv, nomer, vyfrve, vysove, vymivmf, vyukve, fl, ident);
                    ident=ident+1; tepv=rasMas(efte, temvhtn, temvctn, ktpn, tepvtn, tsn, temvs, temhq, temcq, ktpq, tepv, nomer, vyfrve, vysove, vymivmf, vyukve, fl, ident);
                    ident=ident+1; ktpq=rasMas(efte, temvhtn, temvctn, ktpn, tepvtn, tsn, temvs, temhq, temcq, ktpq, tepv, nomer, vyfrve, vysove, vymivmf, vyukve, fl, ident);
                    ident=ident+1; temvs=rasMas(efte, temvhtn, temvctn, ktpn, tepvtn, tsn, temvs, temhq, temcq, ktpq, tepv, nomer, vyfrve, vysove, vymivmf, vyukve, fl, ident);
                    nomer=nomer+1;
                end
            end
        end
    end
    n=length(efte);
    m=nomer;
    mv=zeros(m,n);
    vymas=zeros(1,n);
    for k=1:n
        for j=1:m
            for i=1:n
                if (abs(temvs(j,i)-efte(k))<=hf)
                    mv(j,i)=hf;
                end
            end
        end
        vymas(k)=usrmas(mv, temvs, temhq, temcq, ktpq, tepv, efte, v);
        for j=1:m
            for i=1:n
                mv(j,i)=0;
            end
        end
    end
    v=v';
switch (v)
    case (0)
        temc=vymas;
    case (1)
        temh=vymas;
    case (2)
        tepv=vymas;
    case (3)
        ktpq=vymas;
    case (4)
        temvs=vymas;
end
    vyma=vymas;
end
%-----
function [ vyma ] = rasMas(efte, temvhtn, temvctn, ktpn, tepvtn, tsn, temvs, temhq, temcq, ktpq, tepv, nomer, vyfrve, vysove, vymivmf, vyukve, f, v)
hf=1;
e=1e-1;
    for k=1:length(efte)
            tf=efte(k);
                    for w=1:length(tsn)
                        temr=tsn(w);
                        if (abs(tf-temr)<=hf)
                            temvs(nomer,w)=temr;
                            temhq(nomer,w)=temr*temvhtn(w); 
                            temcq(nomer,w)=temr*temvctn(w); 
                            ktpq(nomer,w)=temr*ktpn(w);
                            tepv(nomer,w)=tepvtn(w)*temr;
                        break; 
                        end
                    end
    end
    if ((f>0) && (vyfrve==2))
                        vyfv=vyfrve';
                        vysv=vysove';
                        if (vyfv<e)
                        vymivmf=vymivmf';
                        else
                        vyuv=vyukve';
                        end
                        temvctn=temvctn;
                        temvhtn=temvhtn;
                        tepvtn=tepvtn;
                        ktpn=ktpn;
                        tsn=tsn;
    end
switch (v)
    case (0)
        vymas=temcq;
    case (1)
        vymas=temhq;
    case (2)
        vymas=tepv;
    case (3)
        vymas=ktpq;
    case (4)
        vymas=temvs;
end
vyma=vymas;
end
%-----
function u = usrmas(mv, temvs, temhq, temcq, ktpq, tepv, efte, v)
hf=1;
n=length(efte);
s=0;
t=0;
e=1e-6;
for i=1:n
    for j=1:n
        if (mv(i,j)>e)
            switch (v)
    case (0)
        s=s+temcq(i,j);
    case (1)
        s=s+temhq(i,j);
    case (2)
        s=s+tepv(i,j);
    case (3)
        s=s+ktpq(i,j);
    case (4)
        s=s+temvs(i,j);
            end
        end
    end
end
for i=1:n
    for j=1:n
        if (mv(i,j)>e)
            if (v<4)
            t=t+temvs(i,j);
            else
                t=t+hf;
            end
        end
    end
end
if (t>e) 
    u=s/t;
else
    u=0;
end
end
%------
function [ oa ] = opredTempHolGor(ektp, ete, h0, v, efte) %моделирование процесса теплообмена в образце %n - длина массива ektp, l - длина массивов temvs, qob, ktp, temvh, temvc, ni - длина efte, qob - плотность теплового потока, которую может создать лабораторная установка 
    ektp=ektp;
    ete=ete;
	hf=1e0; te0=273.15; tn=22.0+te0;
    f=1;
	k=2; qon=PoiskZavVelTem(k, efte, h0, f);
    f=-1;
	k=4; tena=PoiskZavVelTem(k, efte, h0, f);
	ni=length(efte); e=1e-10; 
	koeq=koefPribSha(qon, tena);
    dmkoef=length(koeq);
    temvc=zeros(1,ni);
    qob=zeros(1,ni);
    temvh=zeros(1,ni);
    temvs=zeros(1,ni);
    ktp=zeros(1,ni);
		for k=1:ni
			ts=efte(k); g=0; p=g;
			for j=1:dmkoef
                g=g+(ts^p)*koeq(j); 
                p=p+hf;
            end
			qob(k)=g;
            temvc(k)=0; temvh(k)=0; temvs(k)=0; ktp(k)=0;
        end
        ete=ete;
        ektp=ektp;
        for k=1:ni
            qo=qob(k); laef=0.0;
			if (qo>e)
                etem=efte(k); laef=opredKTPTKTochSha(ektp,ete,etem); 
                dt=qo*h0/laef;
                thol=etem-dt/2e0;
                tgor=thol+dt;
            else
                thol=0.0; tgor=0.0; qob(k)=0.0; ktp(k)=0.0; 
            end
             temvs(k)=(tgor+thol)/2e0; temvc(k)=thol; temvh(k)=tgor; ktp(k)=laef; qob(k)=qo;
        end
        ktp=ktp;
        temvs=temvs;
    vybvesch=1; j=-1; 
    k=0; temvcn=vydelPol(temvc, temvh, qob, ktp, temvs, k, j, vybvesch);
    k=k+1; temvhn=vydelPol(temvc, temvh, qob, ktp, temvs, k, j, vybvesch);
    k=k+1; qobn=vydelPol(temvc, temvh, qob, ktp, temvs, k, j, vybvesch);
    k=k+1; ktpn=vydelPol(temvc, temvh, qob, ktp, temvs, k, j, vybvesch);
    k=k+1; temvsn=vydelPol(temvc, temvh, qob, ktp, temvs, k, j, vybvesch);
        switch (v)
            case (0)
                vm=temvcn;
            case (1)
                vm=temvhn;
            case (2)
                vm=qobn;
            case (3)
                vm=ktpn;
            case (4)
                vm=temvsn;
        end
        vm=vm;
	oa=vm;
end
%----------------
function [ oa ] = napMasEKTP(vyfrve, vysove, vymeizvemafr, te, ident, h, vyukve)
	ko=1e-4; F=13.85*ko; hf=1e0;
	if ((vyfrve==0) || (vyfrve==2)) %фракция 2-0,7 мм или фракция 1,6-0,35 мм
		th207=arrTemHigh207(); temvh=th207;
		tc207=arrTemCold207(); temvc=tc207;
		tp207=arrTepPot207(); tepv=tp207/F;
		n207=length(th207); temvs=(temvh+temvc)/2e0;
		for k=1:n207
		ktp(k)=abs(temvh(k)-temvc(k))/h; 
        ktp(k)=tepv(k)/ktp(k);
        end
        if (vysove==0) %если выбран исходный вермикулит
		if (vymeizvemafr==0) %нестационарный метод - установка Netzsch
		k=0; koefc=oprkoefKTPiskhchao(k, te, h);
        k=1; koefh=oprkoefKTPiskhchao(k, te, h); 
        k=2; koefq=oprkoefKTPiskhchao(k, te, h);
		k=3; koeft=oprkoefKTPiskhchao(k, te, h);
        k=4; koefs=oprkoefKTPiskhchao(k, te, h);
        elseif (vymeizvemafr==1) %данные 2020 года - ГОСТ 12170 - стационарный метод
		ktp=arrKTP_2020();
        temvh=arrTem1_2020();
		temvc=arrTem2_2020();
		temvs=arrTem3_2020();
        tst=(temvh+temvc)/2e0; ts=tst; temvs=ts;
		ts=danIskh207(ts, tst);
        ktp=danIskh207(ktp, tst);
        temvh=danIskh207(temvh, tst);
		temvc=danIskh207(temvc, tst);
        temvs=danIskh207(temvs, tst);
		koeft=koefPribSha(ktp, ts); 
		koefh=koefPribSha(temvh, ts);
		koefc=koefPribSha(temvc, ts);
		koefs=koefPribSha(temvs, ts);
        dmk=length(koeft); q=length(ts); tepv=0;
		for k=1:q
			s=0; r=0; 
            for j=1:dmk
                s=s+koeft(j)*(ts(k)^r); 
                r=r+hf;
            end
			tepv(k)=s*abs(temvh(k)-temvc(k))/h;
        end
			koefq=koefPribSha(tepv, ts); 
        elseif (vymeizvemafr==2) %данные 2019 года - ГОСТ 12170
		koefq=danPoTemTepl2071(temvs, tepv); 
        koefh=danPoTemTepl2071(temvs, temvh);
		koefc=danPoTemTepl2071(temvs, temvc);
        koeft=danPoTemTepl2071(temvs, ktp);
		koefs=danPoTemTepl2071(temvs, temvs);
        end
        elseif (vysove==1) %после повторных измерений
			koefc=danPoTemTepl2072(temvs, temvc); 
            koefh=danPoTemTepl2072(temvs, temvh); 
			koefq=danPoTemTepl2072(temvs, tepv); 
            koeft=danPoTemTepl2072(temvs, ktp); 
			koefs=danPoTemTepl2072(temvs, temvs);
        elseif (vysove == 2) %после обжига при 1000 °С
			koefc=danPoTemTepl2073(temvs, temvc); 
            koefh=danPoTemTepl2073(temvs, temvh); 
			koefq=danPoTemTepl2073(temvs, tepv); 
            koeft=danPoTemTepl2073(temvs, ktp); 
			koefs=danPoTemTepl2073(temvs, temvs);
        elseif (vysove == 3) %после повторного обжига при 1000 °С
			koefc=danPoTemTepl2074(temvs, temvc); 
            koefh=danPoTemTepl2074(temvs, temvh); 
			koefq=danPoTemTepl2074(temvs, tepv); 
            koeft=danPoTemTepl2074(temvs, ktp); 
			koefs=danPoTemTepl2074(temvs, temvs); 
        end
    elseif (vyfrve==1) %фракция 8-4 мм
		tc84=arrTemCold84(); temvc=tc84;
		th84=arrTemHigh84(); temvh=th84;
        tp84=arrTepPot84(); tepv=tp84/F; 
		nf84=length(tc84); n84=nf84; nn=nf84; 
        temvs=(th84+tc84)/2e0; 
		for k=1:n84
		ktp(k)=abs(th84(k)-tc84(k))/h; 
		ktp(k)=tepv(k)/ktp(k);
        end
        if (vyukve==1) %плоско-параллельная засыпка
			if (vysove==0) %исходный
				koefc=danPoTemTepl840(temvs, temvc); 
                koefh=danPoTemTepl840(temvs, temvh); 
				koefq=danPoTemTepl840(temvs, tepv); 
                koeft=danPoTemTepl840(temvs, ktp); 
				koefs=danPoTemTepl840(temvs, temvs);
            elseif (vysove==1) %после повторных измерений
				koefc=danPoTemTepl842(temvs, temvc); 
                koefh=danPoTemTepl842(temvs, temvh); 
				koefq=danPoTemTepl842(temvs, tepv); 
                koeft=danPoTemTepl842(temvs, ktp); 
				koefs=danPoTemTepl842(temvs, temvs);
            elseif (vysove==2) %после обжига
				koefc=danPoTemTepl845(temvs, temvc); 
                koefh=danPoTemTepl845(temvs, temvh); 
				koefq=danPoTemTepl845(temvs, tepv); 
                koeft=danPoTemTepl845(temvs, ktp); 
				koefs=danPoTemTepl845(temvs, temvs); 
            end
        elseif (vyukve==2) %вертикальная засыпка
			if (vysove==0) %исходный
				koefc=danPoTemTepl841(temvs, temvc); 
                koefh=danPoTemTepl841(temvs, temvh); 
				koefq=danPoTemTepl841(temvs, tepv); 
                koeft=danPoTemTepl841(temvs, ktp); 
				koefs=danPoTemTepl841(temvs, temvs);
            elseif (vysove == 1) %после повторных измерений
				koefc = danPoTemTepl844(temvs, temvc); 
                koefh = danPoTemTepl844(temvs, temvh); 
				koefq = danPoTemTepl844(temvs, tepv); 
                koeft = danPoTemTepl844(temvs, ktp); 
				koefs = danPoTemTepl844(temvs, temvs); 
            elseif (vysove == 2) %после обжига
				koefc = danPoTemTepl843(temvs, temvc);
                koefh = danPoTemTepl843(temvs, temvh);
				koefq = danPoTemTepl843(temvs, tepv);
                koeft = danPoTemTepl843(temvs, ktp);
				koefs = danPoTemTepl843(temvs, temvs);
            end
        end
    end
    switch (ident)
        case (0)
            vm=koefc;
        case (1)
            vm=koefh;
        case (2)
            vm=koefq;
        case (3)
            vm=koeft;
        case (4)
            vm=koefs;
    end
    nt=length(te); nvm=length(vm); 
    ma=zeros(1,nt);
    for k=1:nt
        t=te(k); s=0; r=s;
	for j=1:nvm
        s=s+vm(j)*(t^r);
        r=r+hf;
    end
    ma(k)=s;
	end
	oa=ma;
end
%----------------
function [ oa ] = oprkoefKTPiskhchao(v, efte, h) %vmiv - выбор метода измерений %0 - установка Netzsch - нестационарный метод
        mt=arrTem_Netzsch();
		ktpv=arrKTP_Netzsch();
		k=0; thvn=opredTempHolGor(ktpv, mt, h, k, efte);
        k=k+1; tgvn=opredTempHolGor(ktpv, mt, h, k, efte); 
		k=k+1; qovn=opredTempHolGor(ktpv, mt, h, k, efte);
        k=k+1; ktpvn=opredTempHolGor(ktpv, mt, h, k, efte);
        k=k+1; tsvn=opredTempHolGor(ktpv, mt, h, k, efte);
            f=-1;
            vybvesch=1;
    k=0; thvnn=vydelPol(thvn, tgvn, qovn, ktpvn, tsvn, k, f, vybvesch);
    k=k+1; tgvnn=vydelPol(thvn, tgvn, qovn, ktpvn, tsvn, k, f, vybvesch);
    k=k+1; qovnn=vydelPol(thvn, tgvn, qovn, ktpvn, tsvn, k, f, vybvesch);
    k=k+1; ktpvnn=vydelPol(thvn, tgvn, qovn, ktpvn, tsvn, k, f, vybvesch);
    k=k+1; tsvnn=vydelPol(thvn, tgvn, qovn, ktpvn, tsvn, k, f, vybvesch);
        thvnn=thvnn;
        tsvnn=tsvnn;
		koefc=koefPribSha(thvnn, tsvnn);
        %le=length(thvnn);
        %koefc=koefPribVer(thvn, tsvn);
        koefh=koefPribSha(tgvnn, tsvnn);
		koefq=koefPribSha(qovnn, tsvnn);
        koeft=koefPribSha(ktpvnn, tsvnn); 
		koefs=koefPribSha(tsvnn, tsvnn); 
    switch (v)
        case (0)
            vm=koefc;
        case (1)
            vm=koefh;
        case (2)
            vm=koefq;
        case (3)
            vm=koeft;
        case (4)
            vm=koefs;
    end
    oa=vm;
end
%----------------
function [ oa ] = danIskh207(ma, x)
		s=0; t=0; q=3;
		for k=1:q
            s=s+ma(k)*x(k); 
            t=t+x(k);
        end
        s=s/t;
		k=1; m(k)=s;
		k=2; p=4; m(k)=ma(p);
		k=3; p=5; q=6; m(k)=(ma(p)*x(p)+ma(q)*x(q))/(x(p)+x(q)); %600, 800 и 1000 °C
oa=m;
end
%---------
function [ tem ] = arrTem_Netzsch()
te0=273.15;
tem=[27.967, 93.833, 192.5, 341.667, 491.467, 641.2, 790.933, 991.133]+te0;
end
function [ ktp ] = arrKTP_Netzsch()	
ktp=[70, 80, 103, 150, 207, 283, 373, 477]*1e-3; 
end
function [ a ] = arrKTP_2020()
a=[0.175566644058715, 0.176801537812368, 0.179324717653617, 0.211768953068592, 0.237194543621728, 0.237231989760775];
end
function [ a ] = arrTem1_2020()
te0=273.15;
a=[585, 600, 585, 800, 1000, 1000]+te0;
end
function [ a ] = arrTem2_2020()
te0=273.15;
a=[119.75, 138, 129.5, 200, 273, 261]+te0;
end
function [ a ] = arrTem3_2020()
te0=273.15;
a=[377, 396, 383.5, 548, 703, 697.25]+te0;
end
%---------
function [ isdan ] = danPoTemTepl840(temvs, temvh) %Засыпка плоско-параллельная, исходный
	e=1e-6;
	p=15; q=17; tem1=temvs(p)+temvs(q);
	p=16; q=18; tem2=temvs(p)+temvs(q);
	p=15; q=17; tho1=(temvh(p)*temvs(p)+temvh(q)*temvs(q))/tem1;
	p=16; q=18; tho2=(temvh(p)*temvs(p)+temvh(q)*temvs(q))/tem2;
	tem1=tem1/2e0; tem2=tem2/2e0;
	k1=0;
	if (abs(tem2-tem1)>e) 
        k1=(tho2-tho1)/(tem2-tem1); 
    end
	k2=tho2-k1*tem2;
	isdan=[k2, k1]; 
end
function [ isdan ] = danPoTemTepl841(temvs, temvh) %Засыпка вертикальная, исходный
	tc1=0; tc2=0; e=1e-6;
	p=5; q=7; r=11; tem1=temvs(p)+temvs(q)+temvs(r);
	p=6; q=8; r=12; tem2=temvs(p)+temvs(q)+temvs(r);
	p=5; q=7; r=11; 
    if (abs(tem1)>e) 
        tc1=(temvh(p)*temvs(p)+temvh(q)*temvs(q)+temvh(r)*temvs(r))/tem1;
    end
	p=6; q=8; r=12; 
    if (abs(tem2)>e) 
        tc2=(temvh(p)*temvs(p)+temvh(q)*temvs(q)+temvh(r)*temvs(r))/tem2;
    end
	tem1=tem1/3e0; tem2=tem2/3e0;
	k1=0;
    if (abs(tem2-tem1)>e) 
        k1=(tc2-tc1)/(tem2-tem1); 
    end
    k2=tc2-k1*tem2;
	isdan=[k2, k1]; 
end
function [ isdan ] = danPoTemTepl842(temvs, temvh) %Засыпка плоско-параллельная, повторы
	e=1e-6;
	p=19; q=21; r=23; tem1=temvs(p)+temvs(q)+temvs(r);
	p=20; q=22; r=24; tem2=temvs(p)+temvs(q)+temvs(r);
	p=19; q=21; r=23; 
    tc1=0;
    if (abs(tem1)>0) 
        tc1=(temvh(p)*temvs(p)+temvh(q)*temvs(q)+temvh(r)*temvs(r))/tem1; 
    end
	p=20; q=22; r=24; tc2=0;
    if (abs(tem2)>e) 
        tc2=(temvh(p)*temvs(p)+temvh(q)*temvs(q)+temvh(r)*temvs(r))/tem2; 
    end
	tem1=tem1/3e0; tem2=tem2/3e0;
	k1=0;
    if (abs(tem2-tem1)>e) 
        k1=(tc2-tc1)/(tem2-tem1); 
    end
    k2=tc2-k1*tem2;
    isdan=[k2,k1];
end
function [ isdan ] = danPoTemTepl843(temvs, temvh) %Засыпка вертикальная, после обжига при 1000 °С
	e=1e-6;
	p=1; q=3; tem1=temvs(p)+temvs(q);
	p=2; q=4; tem2=temvs(p)+temvs(q);
	p=1; q=3; tc1=0;
    if (abs(tem1)>e) 
        tc1=(temvh(p)*temvs(p)+temvh(q)*temvs(q))/tem1;
    end
	p=2; q=4; tc2=0; 
    if (abs(tem2)>e) 
        tc2=(temvh(p)*temvs(p)+temvh(q)*temvs(q))/tem2;
    end
	tem1=tem1/2e0; 
    tem2=tem2/2e0;
	k1=0; k2=0; 
    if (abs(tem2-tem1)>e) 
        k1=(tc2-tc1)/(tem2-tem1); 
    end
    k2=tc2-k1*tem2;
	isdan=[k2,k1]; 
end
function [ isdan ] = danPoTemTepl844(temvs, temvh) %Засыпка вертикальная, повторы
	e=1e-6;
	p=9; q=13; tem1=temvs(p)+temvs(q);
	p=10; q=14; tem2=temvs(p)+temvs(q);
	p=9; q=13; tc1=0;
    if (abs(tem1)>e) 
        tc1=(temvh(p)*temvs(p)+temvh(q)*temvs(q))/tem1;
    end
	p=10; q=14; tc2=0;
    if (abs(tem2)>e) 
        tc2=(temvh(p)*temvs(p)+temvh(q)*temvs(q))/tem2;
    end
	tem1=tem1/2e0; tem2=tem2/2e0;
    k1=0;
	if (abs(tem2-tem1)>e) 
        k1=(tc2-tc1)/(tem2-tem1); 
    end
    k2=tc2-k1*tem2;
	isdan=[k2,k1];
end
function [ isdan ] = danPoTemTepl845(temvs, temvh) %Засыпка плоско-параллельная, после обжига при 1000 °С
	e=1e-6;
	k=25; tem1=temvs(k);
	k=26; tem2=temvs(k);
	k=25; tc1=0;
    if (abs(tem1)>e) 
        tc1=(temvh(k)*temvs(k))/tem1;
    end
	k=26; tc2=0;
    if (abs(tem2)>e) 
        tc2=(temvh(k)*temvs(k))/tem2;
    end
	k1=0;
    if (abs(tem2-tem1)>e) 
        k1=(tc2-tc1)/(tem2-tem1);
    end
    k2=tc2-k1*tem2;
	isdan=[k2,k1]; 
end
%---------
function [ isdan ] = danPoTemTepl2071(temvs, tepv) %Засыпка исходная, фракция 2-0,7 мм
    tepv1=0; tepv2=0; e=1e-6;
	k=1; tem1=temvs(k);
	k=2; tem2=temvs(k);
	k=1; 
    if (abs(tem1)>e) 
        tepv1=(tepv(k)*temvs(k))/tem1;
    end
	k=2; 
    if (abs(tem2)>e) 
        tepv2=(tepv(k)*temvs(k))/tem2;
    end
	k1=0;
    if (abs(tem2-tem1)>e) 
        k1=(tepv2-tepv1)/(tem2-tem1); 
    end
    k2=tepv2-k1*tem2;
	isdan=[k2,k1];
end
function [ isdan ] = danPoTemTepl2072(temvs, tepv) %Фракция 2-0,7 мм (повторные измерения)
	e=1e-6;
	p=3; q=5; tem1=temvs(p)+temvs(q);
	p=4; q=6; tem2=temvs(p)+temvs(q);
	p=3; q=5; tepv1=0; 
    if (abs(tem1)>e) 
        tepv1=(tepv(p)*temvs(p)+tepv(q)*temvs(q))/tem1;
    end
	p=4; q=6; tepv2=0;
    if (abs(tem2)>e) 
        tepv2=(tepv(p)*temvs(p)+tepv(q)*temvs(q))/tem2;
    end
	tem1=tem1/2e0; tem2=tem2/2e0;
	k1=0;
    if (abs(tem2-tem1)>e) 
        k1=(tepv2-tepv1)/(tem2-tem1); 
    end
    k2=tepv2-k1*tem2;
	isdan=[k2,k1]; 
end
function [ isdan ] = danPoTemTepl2073(temvs, tepv) %Фракция 2-0,7 мм, после обжига при 1000 °С
	e=1e-6;
	p=7; q=9; tem1=temvs(p)+temvs(q);
	p=8; q=10; tem2=temvs(p)+temvs(q);
	p=7; q=9; tepv1=0;
    if (abs(tem1)>e)
    tepv1=(tepv(p)*temvs(p)+tepv(q)*temvs(q))/tem1;
    end
	p=8; q=10; tepv2=0;
    if (abs(tem2)>e)
    tepv2=(tepv(p)*temvs(p)+tepv(q)*temvs(q))/tem2;
    end
	tem1=tem1/2e0; tem2=tem2/2e0;
	k1=0;
    if (abs(tem2-tem1)>e) 
        k1=(tepv2-tepv1)/(tem2-tem1); 
    end
	k2=tepv2-k1*tem2;
	isdan=[k2,k1]; 
end
function [ isdan ] = danPoTemTepl2074(temvs, tepv) %Фракция 2-0,7 мм, после повторного обжига при 1000 °С
	tepv1=0; tepv2=0; e=1e-6;
	p=3; tem1=temvs(p);
	p=4; tem2=temvs(p);
	p=3; 
    if (abs(tem1)>e)
    tepv1=(tepv(p)*temvs(p))/tem1;
    end
    p=4; 
    if (abs(tem2)>e)
	tepv2=(tepv(p)*temvs(p))/tem2;
    end
	k1=0;
    if (abs(tem2-tem1)>e) 
        k1=(tepv2-tepv1)/(tem2-tem1); 
    end
	k2=tepv2-k1*tem2;
	isdan=[k2,k1]; 
end
%---------
function [ temcol ] = arrTemCold84()
te0=273.15;
temcol=[168.0, 369.0, 168.0, 369.0, 148.6, 356.5, 184.0, 396.0, 148.75 ...
            350.0, 166.0, 375.0, 171.0, 383.5, 106.75, 242.0, 123.0, 294.0 ... 
            111.0, 240.0, 109.0, 232.75, 127.0, 291.0, 221.0, 443.75]+te0;
end
function [ temh ] = arrTemHigh84()
te0=273.15;
temh=[545.0, 927.0, 530.0, 925.0, 560.0, 950.0, 560.0, 900.0, 558.0 ...
           950.0, 540.0, 920.0, 540.0, 920.0, 600.0, 1000.0, 580.0, 1000.0 ...
           587.0, 1000.0, 590.0, 1000.0, 580.0, 1000.0, 483.0, 850.0]+te0;
end
function [ tep ] = arrTepPot84()
tep=[3.9805, 9.5532, 3.6447, 9.4779, 3.55732, 11.4997, 4.4624, 11.6021 ...
    4.7766, 11.5016, 3.4023, 10.2068, 3.92812, 11.17333, 3.5144, 8.6593 ...
    2.977, 8.0448, 3.352, 9.218, 3.0313, 7.7946, 3.0671, 6.1342, 1.73466, 4.32967];
end
function [ temcol ] = arrTemCold207()
te0=273.15;
temcol=[109.0, 235.0, 101.0, 199.0, 108.75, 238.0, 124.0, 266.0, 111.0, 262.0]+te0;
end
function [ temh ] = arrTemHigh207()
te0=273.15;
temh=[585.0, 1000.0, 603.0, 1000.0, 603.0, 1000.0, 571.5, 1000.75, 583.0, 1000.0]+te0;
end
function [ tepot ] = arrTepPot207()
tepot=[3.9596, 8.6377, 2.3003, 5.3674, 3.56149, 7.123, 2.12992, 7.6956, 2.3003, 6.9009];
end
%---------
function [ vyma ] = vydelPol(temvcs, temvhs, qos, ektpvs, temvss, v, f, vybves)
tev0 = 273.15; ko=1e2; n=length(qos); e=1e-6; q=0;
switch (vybves)
    case (0)
matepr=11.5*ko;
    case (1)
matepr = 11.0*ko;
end
matepr=matepr + tev0;
if (f>0)
for k = 1:n
if ((temvcs(k)>e) && (temvhs(k)>e) && (qos(k)>e) && (ektpvs(k)>e) && (temvss(k)>e) && (temvhs(k)<matepr) && (temvcs(k)<matepr) && (temvss(k)<matepr)) 
        q=q+1; 
end
end
qn = q;
if (qn==0)
    temvcs=temvcs;
    temvhs=temvhs;
    qos=qos;
    ektpvs=ektpvs;
    temvss=temvss;
end
q = 1; vm=zeros(1, qn);
if (qn>0)
for k=1:n
if ((temvcs(k)>e) && (temvhs(k)>e) && (qos(k)>e) && (ektpvs(k)>e) && (temvss(k)>e) && (temvhs(k)<matepr) && (temvcs(k)<matepr) && (temvss(k)<matepr)) 
                switch (v)
                case (0) 
                    vm(q) = temvcs(k);
                case (1) 
                    vm(q) = temvhs(k);
                case (2) 
                    vm(q) = qos(k);
                case (3) 
                    vm(q) = ektpvs(k);
                case (4) 
                    vm(q) = temvss(k);
                end
            q=q+1;
end
end
else
switch (v)
                case (0) 
                    vm = temvcs;
                case (1) 
                    vm = temvhs;
                case (2) 
                    vm = qos;
                case (3) 
                    vm = ektpvs;
                case (4) 
                    vm = temvss;
end
end
else
for k = 1:n
if ((temvcs(k)>e) && (temvhs(k)>e) && (qos(k)>e) && (ektpvs(k)>e) && (temvss(k)>e)) 
        q=q+1; 
end
end
qn = q;
if (qn==0)
    temvcs=temvcs;
    temvhs=temvhs;
    qos=qos;
    ektpvs=ektpvs;
    temvss=temvss;
end
q = 1; vm=zeros(1, qn);
if (qn>0)
for k=1:n
if ((temvcs(k)>e) && (temvhs(k)>e) && (qos(k)>e) && (ektpvs(k)>e) && (temvss(k)>e)) 
                switch (v)
                case (0) 
                    vm(q) = temvcs(k);
                case (1) 
                    vm(q) = temvhs(k);
                case (2) 
                    vm(q) = qos(k);
                case (3) 
                    vm(q) = ektpvs(k);
                case (4) 
                    vm(q) = temvss(k);
                end
            q=q+1;
end
end
else
switch (v)
                case (0) 
                    vm = temvcs;
                case (1) 
                    vm = temvhs;
                case (2) 
                    vm = qos;
                case (3) 
                    vm = ektpvs;
                case (4) 
                    vm = temvss;
end
end    
end
	vyma=vm;
end
%---------
function [ vm ] = vydelPolN(thol, tgor, qob, ektpv, ete, v, fl, vyve)
te0=273.15; matepr=0.0; ko=1e2; vybmar=0;
    switch (vyve)
        case (0)
            if ((vybmar>=4) && (vybmar<=6)) matepr=11.5*ko;
            elseif ((vybmar>=7) && (vybmar<=8)) matepr=13.0*ko;
            elseif ((vybmar>=9) && (vybmar<=10)) matepr=12.7*ko;
            elseif ((vybmar>=11) && (vybmar<=13)) matepr=13.0*ko;
            end
        case (1)
            matepr=11.0*ko; 
        case (2)
            matepr=11.0*ko; 
        case (3)
            if (vybmar==4) matepr=1e1*ko; 
            elseif (vybmar==5) matepr=10.5*ko; 
            elseif ((vybmar>=6) && (vybmar<=7)) matepr=11.0*ko; 
            elseif ((vybmar>=8) && (vybmar<=10)) matepr=11.5*ko;
            end
    end
    matepr=matepr+te0;
    t=-1; tholn=podvydelPol(thol, fl, t, matepr);
    t=-t; tgorn=podvydelPol(tgor, fl, t, matepr);
    t=-t; qobn=podvydelPol(qob, fl, t, matepr);
    ektpvn=podvydelPol(ektpv, fl, t, matepr);
    eten=podvydelPol(ete, fl, t, matepr);
    nk=length(ete);
    nom=0;
    for k=1:nk
        nom(k)=tholn(k)*tgorn(k)*qobn(k)*ektpvn(k)*eten(k);
    end
    if (length(nom)==1)
        disp('1');
    end
    q=1; tho=0; tgo=0; qo=0; ektp=0; et=0;
    for k=1:nk
        x=nom(k);
        if (x==1)
            tho(q)=thol(k);
            tgo(q)=tgor(k);
            qo(q)=qob(k);
            ektp(q)=ektpv(k);
            et(q)=ete(k);
            q=q+1;
        end
    end
    switch (v)
        case (0)
        vyma=tho;
        case (1)
        vyma=tgo;
        case (2)
        vyma=qo;
        case (3)
        vyma=ektp;
        case (4)
        vyma=et;
    end
    vm=vyma;
end
%---------
function [ vm ] = podvydelPol(po, fl, t, templa)
e=1e-6; nk=length(po);
for k=1:nk
    nom(k)=1;
end
for k=1:nk
	if (po(k)<e)
        nom(k)=0;
    end
end
if ((fl>0) && (t>0))
            for k=1:nk
                if (po(k)>templa)
                    nom(k)=0;
                end
            end
end
vm=nom;
end
function [ ko ] = koefPribSha(ktp, te)
	yx2=0.0; yx=0.0; x4=0.0; x3=0.0; x2=0.0; x=0.0; y=0.0; le=length(ktp); p=le;
	for k=1:le
		yx2=yx2+ktp(k)*(te(k)^2e0); 
		yx=yx+ktp(k)*te(k); 
		y=y+ktp(k);
		x4=x4+(te(k)^4e0); 
		x3=x3+(te(k)^3e0); 
		x2=x2+(te(k)^2e0); 
		x=x+te(k);
    end %применение метода наименьших квадратов
	b=[yx2, yx, y]; 
	A = [x4, x3, x2; x3, x2, x; x2, x, p];
	de = det(A);
    A1 = [yx2, yx, y; x3, x2, x; x2, x, p];
	de1 = det(A1);
    A2 = [x4, x3, x2; yx2, yx, y; x2, x, p];
	de2 = det(A2);
    A3 = [x4, x3, x2; x3, x2, x; yx2, yx, y];
	de3 = det(A3);
	ko=[de3/de, de2/de, de1/de];
end
function [ ko ] = koefPribVer(ktp, te)
	yx2=0; yx=0; x4=0; x3=0; x2=0; x=0; y=0; le=length(ktp); p=le;
	for k=1:le
		yx2=yx2+ktp(k)*(te(k)^2e0); 
		yx=yx+ktp(k)*te(k); 
		y=y+ktp(k);
		x4=x4+(te(k)^4e0); 
		x3=x3+(te(k)^3e0); 
		x2=x2+(te(k)^2e0); 
		x=x+te(k);
    end %применение метода наименьших квадратов
	b=[yx2, yx, y]'; 
    A = [x4, x3, x2; x3, x2, x; x2, x, p];
	ko=fliplr((inv(A)*b)');
end
function ktpo = opredKTPTKTochToch(ktptks, te, temp)
f = 1; p = 0; ep=1e-6; ktp = 0.0; n=length(te); nn=1;
if ((temp>=te(nn)) && (temp<=te(n)))
for k = 1:n
if ((te(k) >= temp) && (f>0))
        p = k; f = 0; break;
end
end
elseif (temp>=te(n))
    p=n; f=0;
elseif (temp<=te(1))
    p=2; f=0;
end
if ((f==0) && (p>1))
    x2=te(p);
    x1=te(p-1);
	dt = x2 - x1;
if (abs(dt) > ep)
    y2=ktptks(p);
    y1=ktptks(p - 1);
    b=y1;
			ko = (y2 - y1) / dt;
            if (p==n)
                b=y2;
            end
			ktp = b + ko*(temp - x1);
else 
    ktp=0;
end
end
ktpo=ktp;
end
%---
function ktpo = opredKTPTKTochSha(ktptks, te, temp)
	ce=length(te); n=ce; f=1; p=0; k=0; nn=1; nk=n; u=-1; d=2;
	e=1e-6; ktp=0.0; ko=0.0; x1=0.0; x2=0.0; y1=0.0; y2=0.0; b=0.0; dt=0.0; 
	if (temp<te(nn)) 
        p=d; f=u;
    elseif (temp>=te(nk)) 
        p=nk; f=u;
    elseif ((temp>=te(nn)) && (temp<te(nk-1)))
	for k = nn:(nk-1)
		if ((te(k+1)>=temp) && (f>0) && (te(k)<temp)) 
            p = k+1; f = u; break;
        end
    end
	elseif ((temp>=te(nk-1)) && (temp<te(nk))) 
        p=nk; f=u;
    end
    f=f';
    temp=temp;
    te=te;
	if ((f==u) && (p>nn))
		x2=te(p);
		x1=te(p - 1);
		dt = x2 - x1;
		if (abs(dt) > e)
			y2=ktptks(p);
			y1=ktptks(p - 1);
			b=y1;
			if (temp>te(nk)) 
                b=y2;
            end
			ko = (y2 - y1) / dt;
			ktp = b + ko*(temp - x1);
        end
    end
	ktpo=ktp;
end