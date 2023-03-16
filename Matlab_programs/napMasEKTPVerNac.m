function [ npv ] = napMasEKTPVerNac(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf, vybvykhmas, ys0, vyuv, vmivmf, ete) %isrp=0; %выбор метода измерений для фракции 2-0,7 мм: 0 - нестационарный, 1 - стационарный %выбор фракции: 0 - фракция 2-0,7 мм, 1 - фракция 8-4 мм, 2 - фракция 1,6-0,35 мм %выбор укладки: 1 - плоскопараллельная, 2 - вертикальная %выбор состояния: 0 - исходное, 1 - после повторных измерений, 2 - после прокаливания при 1000 град °С
format longg; hf = 1e0; 
h=ys0; cem=length(ete); %dete=1e2; te0=273.15; tnac=2e2+te0; cem=11; ete=tnac:dete:tnac+(cem-1)*dete;
switch (vybvesch)
    case (0) %шамот
    vybsha=vybmar; vystsha=0; cems=cem; 
    por=novNapMas(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf);
    kektps=poisMasKoefsha(vybsha); dmkos=length(kektps);
	for k = 1:cems
        t=ete(k); s=0; g=0; 
	for j = 1:dmkos
        s = s + kektps(j)*(t^g); 
        g = g + hf;
    end
    ektps(k) = s; 
    end
    f=1; fl=1; m=1; k=0; thol=opredTempHolGor(ektps, ete, h, k, ete, m, f, fl);
    f=0; fl=0; m=-m; k=1; tgor=opredTempHolGor(ektps, ete, h, k, ete, m, f, fl);   
    k=2; qob=opredTempHolGor(ektps, ete, h, k, ete, m, f, fl);
    k=3; ektps=opredTempHolGor(ektps, ete, h, k, ete, m, f, fl);
    k=4; ete=opredTempHolGor(ektps, ete, h, k, ete, m, f, fl);
    case (3) %КВИ
        vpk=vybmar; cemk=cem; 
    kektpk=poisMasKoefkvi(vpk); dmkok=length(kektp); k=0;
	por=novNapMas(vybvesch, vpk, vyfv, vysv, vpmf, vpkf);
	for k = 1:cemk
        t=ete(k); s=0; g=0; 
	for j = 1:dmkok
        s = s + kektpk(j)*(t^g); 
        g = g + hf;
    end
    ektpk(k) = s; 
    end
    f=1; fl=1; m=1; k=0; thol=opredTempHolGor(ektpk, ete, h, k, ete, m, f, fl);
    f=0; fl=0; m=-m; k=1; tgor=opredTempHolGor(ektpk, ete, h, k, ete, m, f, fl);   
    k=2; qob=opredTempHolGor(ektpk, ete, h, k, ete, m, f, fl);
    k=3; ektpv=opredTempHolGor(ektpk, ete, h, k, ete, m, f, fl);
    k=4; ete=opredTempHolGor(ektpk, ete, h, k, ete, m, f, fl);        
    case (2) %ИТОМ
    vpi=vybitom; cemi=cem; 
	kektpi=poisMasKoefItom(vpi); k=0;
	por=novNapMas(vybvesch, vpi, vyfv, vysv, vpmf, vpkf);
	for k = 1:cemi
        t=ete(k); s=0; g=0; 
	for j = 1:dmkoi
        s = s + kektpi(j)*(t^g); 
        g = g + hf;
    end
    ektpi(k) = s; 
    end
    f=1; fl=1; m=1; k=0; thol=opredTempHolGor(ektpi, ete, h, k, ete, m, f, fl);
    f=0; fl=0; m=-m; k=1; tgor=opredTempHolGor(ektpi, ete, h, k, ete, m, f, fl);   
    k=2; qob=opredTempHolGor(ektpi, ete, h, k, ete, m, f, fl);
    k=3; ektpv=opredTempHolGor(ektpi, ete, h, k, ete, m, f, fl);
    k=4; ete=opredTempHolGor(ektpi, ete, h, k, ete, m, f, fl);        
    case (1) %вермикулит
	ident=0; f=1; thol=napMasEKTP(vyfv, vysv, vmivmf, ete, ident, ys0, vyuv, f);
	ident=1; f=-f; tgor=napMasEKTP(vyfv, vysv, vmivmf, ete, ident, ys0, vyuv, f);
    ident=2; qob=napMasEKTP(vyfv, vysv, vmivmf, ete, ident, ys0, vyuv, f);
    ident=3; ektpv=napMasEKTP(vyfv, vysv, vmivmf, ete, ident, ys0, vyuv, f);
    ident=4; ete=napMasEKTP(vyfv, vysv, vmivmf, ete, ident, ys0, vyuv, f);
end
    k=0; f=1; tholn=vydelPol(thol, tgor, qob, ektpv, ete, k, f, vybvesch);
    k=1; tgorn=vydelPol(thol, tgor, qob, ektpv, ete, k, f, vybvesch);
    k=2; qobn=vydelPol(thol, tgor, qob, ektpv, ete, k, f, vybvesch);
    k=3; ektpvn=vydelPol(thol, tgor, qob, ektpv, ete, k, f, vybvesch);
    k=4; eten=vydelPol(thol, tgor, qob, ektpv, ete, k, f, vybvesch);
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
           
function [ mv ] = PoiskZavVelTem(v, efte, h, fla, flag, fl1, fl2, fl3)
	hf=1e0; e=1e-6; vyukve=1; n=length(efte); c=0; vybvesch=1;
    for k=1:n
    temvs(k)=0; temhq(k)=0; temcq(k)=0; ktpq(k)=0; tepv(k)=0; cemt(k)=0;
    temvct(k)=0; temvht(k)=0; tepvt(k)=0; ktp(k)=0; ts(k)=0; 
    end
	nvyfv=[0,1]; nnvyfv=length(nvyfv); %фракции вермикулита
	nvysv=[0,1,2]; nnvysv=length(nvysv); %состояния вермикулита
	nvmivmf=[1,2]; nnvmivmf=length(nvmivmf); %стационарные методы измерений - 2019 и 2020
	nvyuv=[1,2]; nnvyuv=length(nvyuv); %укладка вермикулита
	for kvf=1:nnvyfv
		vyfrve=nvyfv(kvf);
		for jvsv=1:nnvysv
			vysove=nvysv(jvsv);
			if ((vyfrve==0) || (vyfrve==2))
				for qvmi=1:nnvmivmf
					vymivmf=nvmivmf(qvmi); 
                    ident=0; temvct=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve);
                    ident=1; temvht=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve); 
                    ident=2; tepvt=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve);
                    ident=3; ktp=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve);
                    ident=4; ts=napMasEKTP(vyfrve, vysove, vymivmf, efte, ident, h, vyukve);
                    k=0; f=-1; temvctn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, f, vybvesch);
                    k=1; temvhtn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, f, vybvesch);
                    k=2; tepvtn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, f, vybvesch);
                    k=3; ktpn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, f, vybvesch);
                    k=4; tsn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, f, vybvesch);
                    if ((fla>0) && (flag>0) && (fl1>0) && (fl2>0) && (fl3>0))
                        temvctn=temvctn;
                        temvhtn=temvhtn;
                        tepvtn=tepvtn;
                        ktpn=ktpn;
                        tsn=tsn;
                    end
                    q=length(tsn); d=1; k=1; cfp=tepvt(k); 
                    for k=(d+1):q
                        cft=tepvt(k); 
                        if ((cft<=cfp) && (d>0))
                            d=-1; break;
                        end
                            cfp=cft;
                    end
						if (d>0)
                            for w=1:q 
                                tf=tsn(w); cf=tepvtn(w)*tf;
							for k=1:n
							if (abs(tf-efte(k))<=hf)
                                temvs(k)=temvs(k)+tf;
								temhq(k)=temhq(k)+tf*temvhtn(w); 
                                temcq(k)=temcq(k)+tf*temvctn(w); 
                                ktpq(k)=ktpq(k)+tf*ktpn(w);
								tepv(k)=tepv(k)+cf; 
                                cemt(k)=cemt(k)+hf; 
								break; 
                            end
                            end
                            end
                        end
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
                    k=0; j=-1; temvctn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, j, vybvesch); 
                    k=1; temvhtn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, j, vybvesch); 
                    k=2; tepvtn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, j, vybvesch); 
                    k=3; ktpn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, j, vybvesch); 
                    k=4; tsn=vydelPol(temvct, temvht, tepvt, ktp, ts, k, j, vybvesch); 
                    if ((fla>0) && (flag>0) && (fl1>0) && (fl2>0) && (fl3>0))
                        temvctn=temvctn;
                        temvhtn=temvhtn;
                        tepvtn=tepvtn;
                        ktpn=ktpn;
                        tsn=tsn;
                    end
                    q=length(tsn); d=1; k=1; cfp=tepvt(k); 
                    for k=(d+1):q
                        cft=tepvt(k); 
                        if ((cft<=cfp) && (d>0))
                            d=-1; break;
                        end
                            cfp=cft;
                    end
						if (d>0)
                            for w=1:q 
                                tf=tsn(w); cf=tepvtn(w)*tf;
							for k=1:n
							if (abs(tf-efte(k))<=hf)
								temvs(k)=temvs(k)+tf;
								temhq(k)=temhq(k)+tf*temvhtn(w); 
                                temcq(k)=temcq(k)+tf*temvctn(w); 
                                ktpq(k)=ktpq(k)+tf*ktpn(w);
								tepv(k)=tepv(k)+cf; 
                                cemt(k)=cemt(k)+hf; 
								break; 
                            end
                            end
                            end
                        end
                end
            end
        end
    end
		for k=1:n
			cf=temvs(k); 
            if (cf>e)
				temhq(k)=temhq(k)/cf; 
                temcq(k)=temcq(k)/cf; 
                ktpq(k)=ktpq(k)/cf; 
                tepv(k)=tepv(k)/cf;
            end
			cf=cemt(k); 
            if (cf>e) 
                temvs(k)=temvs(k)/cf; 
            else
                temvs(k)=0;
            end
        end
        if ((fla>0) && (flag>0) && (fl1>0) && (fl2>0) && (fl3>0))
        temcq;
        temhq;
        tepv;
        ktpq;
        temvs;
        end
switch (v)
    case (0)
        mv=temcq;
    case (1)
        mv=temhq;
    case (2)
        mv=tepv;
    case (3)
        mv=ktpq;
    case (4)
        mv=temvs;
end
end
%----------------
function [ oa ] = opredtemphc(effktp, efftem, tgv, thv, qon, h, qob, tgor, thol, v) %n=3 - длина массива коэффициентов приближающего многочлена, tgv - температура горячей стенки, thv - температура холодной стенки, dlma - длина массива ЭКТП
	k=0; j=0; ts=0; g=0; p=0; hf=1e0; r=0; s=0; tego=0; teho=0;
	tesr=(tgv+thv)/2e0;
	kq=koefPribSha(qon, tesr);
	kgo=koefPribSha(tgv, tesr);
	kho=koefPribSha(thv, tesr);
    dlma=length(effktp);
	for k=1:dlma
		ts=efftem(k);
		g=0; p=0; 
        for j=1:n
            g=g+(ts^p)*kq(j);
            p=p+hf;
        end
if (g<0) 
            g=0; 
end
        qob(k)=g;
		p=effktp(k);
		g=fabs(g*h/p); 
		tego=ts+g/2e0; teho=ts-g/2e0;
if ((tego<0) || (teho<0)) 
if (tego<0)
			r=0; s=0; 
            for j=1:n
                r=r+(ts^s)*kgo(j); 
                s=s+hf;
            end
if (r<0) 
                r=0; tego=r; 
end
end
if (teho<0)
			r=0; s=0; 
            for j=1:n
            r=r+(ts^s)*kho(j); 
            s=s+hf;
            end
			if (r<0) 
                r=0; 
            end
            teho=r; 
end
end
	tgor(k)=tego; thol(k)=teho;
    end
switch (v)
        case (0)
            oa=thol;
        case (1)
            oa=tgor;
        case (2)
            oa=qob;
end
end
%----------------
function [ oa ] = opredTempHolGor(ektp, ete, h0, v, efte, fl, f, fla) %моделирование процесса теплообмена в образце %n - длина массива ektp, l - длина массивов temvs, qob, ktp, temvh, temvc, ni - длина efte, qob - плотность теплового потока, которую может создать лабораторная установка 
	nit=1e10; hf=1e0; ep=1e-3; d=1e-4; te0=273.15; tn=22+te0; Thna=tn; dt=hf;
	k=2; fn=1; fln=1; qon=PoiskZavVelTem(k, efte, h0, fn, fln, f, fl, fla);
	k=4; fn=-fn; fln=-fln; tena=PoiskZavVelTem(k, efte, h0, fn, fln, f, fl, fla);
	nf=length(qon); ni=length(efte);
	koeq=koefPribSha(qon, tena);
    dmkoef=length(koeq);
    n=length(ete);
		for k=1:n
			ts=ete(k); g=0; p=g;
			for j=1:dmkoef
                g=g+(ts^p)*koeq(j); 
                p=p+hf;
            end
			qob(k)=g;
            temvc(k)=0; 
            temvh(k)=0; 
            temvs(k)=0; 
            ktp(k)=0;
        end
        g=0;
        for k=1:n
            laef=ektp(k); 
			if (qob(k)>ep)
				p=0; Thnac=Thna+g*dt; del=hf; etem=ete(k); ktp(k)=laef;
				while ((del>ep) && (p<nit))
					thol=Thnac+p*d; %Tg - массив температур горячих стенок, Th - массив температур холодных стенок
					tgor=thol+qob(k)*h0/laef;
					del=abs(2e0*etem-(thol+tgor));
					p=p+hf;
                end
				g=g+hf;
            else
                thol=0; tgor=0; qob(k)=0; ktp(k)=0; 
            end
			 temvs(k)=(tgor+thol)/2e0; temvc(k)=thol; temvh(k)=tgor;
        end
        vybvesch=1;
    k=0; j=1; temvcn=vydelPol(temvc, temvh, qob, ktp, ete, k, j, vybvesch);
    k=1; temvhn=vydelPol(temvc, temvh, qob, ktp, ete, k, j, vybvesch);
    k=2; qobn=vydelPol(temvc, temvh, qob, ktp, ete, k, j, vybvesch);
    k=3; ktpn=vydelPol(temvc, temvh, qob, ktp, ete, k, j, vybvesch);
    k=4; temvsn=vydelPol(temvc, temvh, qob, ktp, ete, k, j, vybvesch);
        if ((fl>0) && (f>0) && (fla>0))
        qon=qon;
        tena=tena;
		ktp=ktp;
        qob=qob;
        ete=ete;
            temvcn=temvcn;
            temvhn=temvhn;
            qobn=qobn;
            ktpn=ktpn;
            temvsn=temvsn;
        end
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
	oa=vm;
end
%----------------
function [ oa ] = napMasEKTP(vyfrve, vysove, vymeizvemafr, te, ident, h, vyukve, fl)
	F=13.85*1e-4; hf=1e0;
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
		k=0; f=1; koefc=oprkoefKTPiskhchao(vymeizvemafr, k, te, h, f, fl);
        k=1; f=-f; koefh=oprkoefKTPiskhchao(vymeizvemafr, k, te, h, f, fl); 
        k=2; koefq=oprkoefKTPiskhchao(vymeizvemafr, k, te, h, f, fl);
		k=3; koeft=oprkoefKTPiskhchao(vymeizvemafr, k, te, h, f, fl);
        k=4; koefs=oprkoefKTPiskhchao(vymeizvemafr, k, te, h, f, fl);
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
function [ oa ] = oprkoefKTPiskhchao(vmiv, v, efte, h, f, fl) %vmiv - выбор метода измерений
	if (vmiv==0) %0 - установка Netzsch - нестационарный метод
		mt=arrTem_Netzsch();
		ktpv=arrKTP_Netzsch();
		m=1; k=0; thvn=opredTempHolGor(ktpv, mt, h, k, efte, m, f, fl);
        m=-m; k=1; tgvn=opredTempHolGor(ktpv, mt, h, k, efte, m, f, fl);   
		k=2; qovn=opredTempHolGor(ktpv, mt, h, k, efte, m, f, fl);
        k=3; ktpvn=opredTempHolGor(ktpv, mt, h, k, efte, m, f, fl);
        k=4; tsvn=opredTempHolGor(ktpv, mt, h, k, efte, m, f, fl);
		koefc=koefPribSha(thvn, tsvn);
        koefh=koefPribSha(tgvn, tsvn);
		koefq=koefPribSha(qovn, tsvn);
        koeft=koefPribSha(ktpvn, tsvn); 
		koefs=koefPribSha(tsvn, tsvn); 
    end
    switch (v)
        case (0)
            oa=koefc;
        case (1)
            oa=koefh;
        case (2)
            oa=koefq;
        case (3)
            oa=koeft;
        case (4)
            oa=koefs;
    end
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
	e=1e-4;
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
	tc1=0; tc2=0; e=1e-4;
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
	e=1e-4;
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
	e=1e-4;
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
	e=1e-3;
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
	e=1e-3;
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
    tepv1=0; tepv2=0; e=1e-3;
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
	e=1e-3;
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
	e=1e-3;
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
	tepv1=0; tepv2=0; e=1e-3;
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

function [ vm ] = vydelPol(thol, tgor, qob, ektpv, ete, v, fl, vyve)
te0=273.15; 
    switch (vyve)
        case (0)
            templa=175*1e1+te0;     
        case (1)
            templa=134*1e1+te0; 
        case (2)
            templa=1750+te0;
        case (3)
            templa=1750+te0;
    end
    t=-1; tholn=podvydelPol(thol, fl, t, templa);
    t=-t; tgorn=podvydelPol(tgor, fl, t, templa);
    t=-t; qobn=podvydelPol(qob, fl, t, templa);
    ektpvn=podvydelPol(ektpv, fl, t, templa);
    eten=podvydelPol(ete, fl, t, templa);
    nk=length(tholn);
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
        if (x>0)
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
function [ kti ] = podpoisMasKoefItom(ktp1, ktp2, t1, t2)
ktpitom(2) = (ktp2 - ktp1) / (te2 - te1); 
ktpitom(1) = ktp1 - ktpitom(2) * te1;
kti=ktpitom;		
end
function [ vm ] = poisMasKoefItom(no, kti)
	f=length(kti); k=0; tei0=273.15; t1=2e2; t2=38e1; kk=1e-2; km=kk/1e1; te200=t1+tei0; te380=t2+tei0;
    switch (no)
        case (0)
		ktpit = [9e0, 12e0]*kk;
		vm=podpoisMasKoefItom(ktpit(1), ktpit(2), te200, te380);
        case (1)
		ktpit = [120, 139]*km; %из Диссертации
		ktpit = [18, 19]*kk; %Данные 2017 года	
        case (2) 
		ktpit = [18.3,19.4]*kk; %из Диссертации
		ktpit = [26, 37]*kk; %Данные 2017 года
        case (3)
        ktpit = [23, 25]*kk; %из Диссертации
		ktpit = [42, 52]*kk; %Данные 2017 года
        otherwise
	disp('Net takoy marki ITOM!');
    end
	vm=podpoisMasKoefItom(ktpit(1), ktpit(2), te200, te380);
end
function [ vm ] = poisMasKoefkvi(vyb, kktp)
	f=length(kktp); t1=25e0; t2=5e2; dt=t2-t1; kn=0.0;
    switch (vyb)
    case (3)
    kktp=[0.068,0.00015]; %350
    case (4) %400
    kktp = [0.082,0.000125]; %1
	kn=(0.156-0.087)/dt; kktp(2)=kn; kn=kn*t2; kktp(1) = 0.156-kn; %2
	kn=(0.155-0.087)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.155-kn; %3 
	kktp = [0.14, 0.00016]; %4 Spirina выбор
	%kn=(0.146-0.087)/dt; kktp(2)=kn; kn=kn*t2; kktp(1) = 0.146-kn; %к вопросу о стандартизации КВИ
    case (5) %500
    kktp = [0.103, 0.0001]; %из Ахтямова, 1991
	kn=(0.165-0.105)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.165-kn; %выбор
	%kn=(0.178-0.105)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.178-kn; %Двухслойные
	%kn=(0.152-0.105)/dt; kktp(2)=kn; kn=kn*t2; kktp(1) = 0.152-kn; //к вопросу о стандартизации КВИ
    case (6) %600
	kktp = [0.116, 0.00015]; %выбор
	kn=(0.201-0.12)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.201-kn;
	kn=(0.195-0.12)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.195-kn;
	kktp = [0.17, 0.00015]; %Spirina
	kn=(0.196-0.12)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.196-kn; %Двухслойные выбор
	%kn=(0.178-0.12)/dt; kktp[1]=kn; kn=kn*t2; kktp[0] = 0.178-kn; %к вопросу о стандартизации КВИ
    case (7) %700
    kktp = [0.146, 0.00017];
	kn=(0.216-0.15)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.216-kn; 
	kn=(0.235-0.15)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.235-kn; %к вопросу о стандартизации КВИ
	kn=(0.251-0.16)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.251-kn; %Двухслойные //выбор
    case (8) %800
    kktp =[0.156, 0.00018];
	kn=(0.226-0.16)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.226-kn; %выбор
	%kn=(0.25-0.16)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.25-kn;
	%kktp = [0.21, 0.00014]; %Spirina
	%kn=(0.23-0.16)/dt; kktp(2)=kn; kn=kn*t2; kktp(1) = 0.23-kn; %к вопросу о стандартизации КВИ
    case (9) %900
    kktp = [0.185, 0.00019];
	kn=(0.23-0.195)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.23-kn; %выбор
	%kn=(0.29-0.195)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.29-kn;
    case (10) %1000
    kktp = [0.246, 0.00025];
	kn=(0.287-0.25)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.287-kn; %выбор
	%kn=(0.36-0.25)/dt; kktp(2) = kn; kn=kn*t2; kktp(1) = 0.36-kn;
	%kn=(0.35-0.25)/dt; kktp(2)=kn; kn=kn*t2; kktp(1) = 0.35-kn; //к вопросу о стандартизации КВИ
    end
	vm=kktp;
end
function [ vm ] = poisMasKoefsha(vrsh)
sfeosha=0; smgosha=0; salosha=0; ssiosha=0;
ko=1e-2; p20=2e1*ko; p24=24.0*ko; p30=3e1*ko; p33=33.0*ko; p16=16.0*ko; p10=1e1*ko;
	switch (vrsh)
        case (0)
            sfeosha = 1.24*ko; smgosha = 0.29*ko; salosha = 41.9*ko; ssiosha = 54.0*ko; porsha=21.8*ko; %Роучка
        case (1)
            sfeosha = 1.21*ko; smgosha = 0.29*ko; salosha = 42.6*ko; ssiosha = 1e0-salosha-sfeosha-smgosha; porsha=11.0144*ko; %ШПД-41
        case (2)
            sfeosha = 1.64*ko; smgosha = 0.36*ko; salosha = 35.9*ko; ssiosha = 59.1e-2;  porsha=25.2*ko; %ШБ-1 2-1
        case (3)
            sfeosha = 1.66*ko; smgosha = 0.4*ko; salosha = 37.3*ko; ssiosha = 57.4*ko;  porsha=26.5*ko; %ШВ-1 1-1
        case (4)
            sfeosha = 1.24*ko; smgosha = 0.29*ko; salosha = 41.9*ko; ssiosha = 54*ko;  porsha=11.5*ko; %ШПД
        case (5)
            sfeosha = 1.54*ko; smgosha = 0.3*ko; salosha = 38.6*ko; ssiosha = 56.5*ko;  porsha=16.5*ko; %ШКУ-32 3-1
    end
	if ((salosha >= 28e-2) && (salosha <= 38e-2))
		if ((porsha >= p20) && (porsha < p24)) 
            vybsha = 0; vystsha = 0; 
		kektp = [-0.435e-9, 0.685e-6, 0.134e-3, 0.725]; 
        elseif ((porsha >= p24) && (porsha < p30))
                vybsha = 0; vystsha = 1; 
		kektp = [-0.867e-9, 1.77e-6, -0.523e-3, 0.806]; %задание коэффициентов - шамот средней пористости
        elseif ((porsha >= p16) && (porsha < p20))
            vybsha = 1; 
		kektp = [-0.397e-9, 0.71e-6, 0.011e-3, 0.851]; %уплотненный шамот
        elseif ((porsha >= p30) && (porsha <= p33))
            vybsha = 2; 
		kektp = [-0.377e-9, 0.918e-6, -0.338e-3, 0.77]; %низкоплотный шамот
        elseif ((porsha >= p10) && (porsha<p16))
            vybsha = 3; 
		kektp = [0.0, -0.607e-6, 1.14e-3, 0.641]; %повышенной плотности шамот
        end
    end
	if ((salosha>38e-2) && (salosha <= 45e-2))
		if ((porsha >= p20) && (porsha < p24))
            vybsha = 0; vystsha = 0; 
		kektp = [-0.124e-9, 0.215e-6, 0.125e-3, 1.01];
        elseif ((porsha >= p24) && (porsha < p30))
            vybsha = 0; vystsha = 1; 
		kektp = [-0.333e-9, 0.805e-6, -0.289e-3, 0.903]; %задание коэффициентов - шамот средней пористости
        elseif ((porsha >= p16) && (porsha < p20))
            vybsha = 1; 
		kektp = [0.0, -0.154e-6, 0.369e-3, 1.03]; %уплотненный шамот
        elseif ((porsha >= p30) && (porsha < p33))
                vybsha = 2; 
		kektp = [-0.377e-9, 0.918e-6, -0.338e-3, 0.77]; %низкоплотный шамот
        elseif ((porsha >= p10) && (porsha < p16))
            vybsha = 3; 
		kektp = [0, -0.141e-6, 0.437e-3, 1.32]; %повышенной плотности шамот
        end
    end
    vm=fliplr(kektp);
end