function [ vmua ] = RaschAlphaTvKar(por, vr, vss)
	l=2; k=0; d=10; q=0; cfk=l; cei=0; razfra=31;
	j=0; jk=1e6; pj=0; h=0;
	x=0.0; y=0.0; p=0.0; xt=0.0; yp=0.0; hf=1e0; te0=273.15; tc=22.0+te0;
	pors=1e-6; lamtem=1448.0; dvpn = lamtem*pors / tc / 4e0;
	pn = zeros(1,d); alx = zeros(1,d); rcf = zeros(1,cfk); e=1e-15;
	rapo=vybFunRasPorpoRaz(por, vr, vss, 0);
    srra=vybFunRasPorpoRaz(por, vr, vss, 1);
    prgr=vybFunRasPorpoRaz(por, vr, vss, 2);
    legr=vybFunRasPorpoRaz(por, vr, vss, 3);
    pt=vybFunRasPorpoRaz(por, vr, vss, 4);
	k=1; srp=pt(k); marp=max(prgr);
    cei=length(srra); 
	x = 0.0;
    koe=1e-6;
    if (srp<hf)
        ko=1e6;
    else ko=1e6;
    end
    srp=srp*ko; marp=ko*marp; 
    for j = 1:razfra
        x = x + srra(j) * rapo(j); 
    end
    pors = x*por / srp;
	pr = 0.0; rmf = 0.0; prf = 0.0; po = 0.0; alsf = zeros(1,cfk);
	po=jk;
	for j = 1:d 
        pn(j) = 0.0; alx(j) = 0.0;
    end
	for j = 1:cfk
        rcf(j) = 0.0; alsf(j) = 0.0;
    end
	for k = 1:l
		x = RasFracXeffSha60(k); %размер частицы
		rcf(k) = x*ko; y = x*pors; %размер поры
		for j = 1:d 
            pn(j) = 0.0;
        end
		for j = 1:jk
			pr = rand(); 
			yp = y*pr; xt = yp*(hf - por) / por;
			p = x / (xt + yp);
			pr = 0.0; 
            for h = 1:d
                xt = pr + hf; 
                if ((p >= pr) && (p < xt)) 
                    pn(h) = pn(h) + hf; pr = xt;
                end
            end
		pr = 0.0; 
        for j = 1:d
            pn(j) = pn(j) / po; 
            pr = pr + pn(j);
        end
		for j = 1:d
            pn(j) = pn(j) / pr;
		for j = 2:d
			p = 0.0; 
            for h = 1:j
                p = p + hf;
            end
			yp = pors*x*koe / (p - hf); xt = (hf - pors)*x*koe / p; 
			if (yp>dvpn)
				lamtem = F0_lamT(yp*tc); 
                if (lamtem<0.0) 
                    lamtem = 0.0; 
                end
                    if (lamtem>hf) 
                        lamtem = hf;
                    end
				alx(j) = BolTochRasAlpha(0, j, yp, xt, tc, ktpvosha, vtesha, etesha, cem, dmkv, snms, Tsh, Rsh, 2 * N, kttk, Ash, 1, Rash, Tash, Aash, Rtsh, Ttsh, Atsh)*lamtem;
            end
			else alx(j) = 0.0;
            end
		pt = RaschRTA(d, xt, 0.0, 0.0, yp, tc, 1, 0, 0, 0); k=1; pp = pt(k); altc = pp(k); 
		k=1; alx(k) = 0.0; k=k+1; alx(k) = altc / (hf - por);
		for j = 1:d
			p = 0.0; 
            for h = 1:j
                p = p + hf;
			yp = pors*x*koe / (p - hf); 
            if (j>1) 
                xt = oprProcSoderpoPoris(rapo, legr, yp, cei); 
            else
                xt = 1e0;
            end
			alsf(k) = alsf(k) + pn(j) * alx(j) * xt;
            end
        end
	x = 0.0; yp = 0.0; 
    for j = 1:cfk
        x = x + alsf(j); 
        yp = yp + hf; 
    end
    x = x / yp; 	
	x = (x / altc) %ослабление  ѕ за счет пористой структуры вермикулита
	return x;
        end