function [ vyma ] = RaschVneshIzluchSha(Tra, Ref, prot, hlr1, hrl1, kost, hlr0, hrl0, fl)
	k = 3; q = 2; e=1e-15;
    khrl = zeros(1,kost); khlr = zeros(1,kost); ri = zeros(1,q);
	hrl1 = zeros(1,kost); hlr1 = zeros(1,kost);
	for k = 1:kost
        khrl(k) = 0.0; 
        khlr(k) = 0.0; 
        hrl1(k) = 0.0; 
        hlr1(k) = 0.0;
    end
    for k = 1:q
        ri(k) = 0.0;
    end
	k=1; khlr(k) = 1e0; k=kost; khrl(k)=1e0; j=k;
	for k = 2:kost
		khlr(k) = khlr(k - 1) * Tra(k - 1);
		khrl(j - 1) = khrl(j) * Tra(j);
		j=j-1;
    end 
    khrl
    khlr
	q = 1; slr = 0.0; srl = 0.0; 
    for j = 1:kost %узнаем, какая доля внешнего излучения дошла до j-ой стенки - учитываем вклад от всех стенок
		hlr1(j) = 0.0; hrl1(j) = 0.0; %внешнее излучение падает с обеих сторон
		tta = 0.0; ttb = 0.0; m = q;
		for k = 1:kost
			w0 = prot(q); 
            w1 = prot(q + 1); 
            w2 = prot(q + 2); 
            w3 = prot(q + 3); 
            wk = w2 + w3; %для внешнего излучения, распространяющегося слева направо
			tt = Ref(k) * w2 + Tra(k) * w3; %доля энергии, отраженная от k-ой стенки, и прошедшая через нее, пришла справа к j-ой стенке
			if (abs(wk)>e)
				tta = tta + tt*khlr(k) * w0 / wk; %доля энергии, отраженная от k-ой стенки, и прошедшая через нее, пришла слева к j-ой стенке
				ttb = ttb + tt*khlr(k) * w1 / wk;
            end %для учета при подсчете лучистой энергии, распространяющейся справа налево
			q = q + 4;
        end
        q = m + 4 * kost; m = q; %считаем, что стенки, до которых дошло излучение, как бы "сами" излучают, из-за того, что на них попало внешнее излучение
		for k = kost:1:(-1)
			w0 = prot(q - 3); 
            w1 = prot(q - 2); 
            w2 = prot(q - 1); 
            w3 = prot(q); 
            wk = w2 + w3; 
            q = q - 4; %для внешнего излучения, распространяющегося справа налево
			tt = Ref(k) * w3 + Tra(k) * w2; %деление на два потока
			if (fabs(wk) > 0.0)
				tta = tta + tt*khrl(k) * w0 / wk; %для учета лучистой энергии, распространяющейся слева направо, падает слева от j-ой стенки
				ttb = ttb + tt*khrl(k) * w1 / wk;
            end %то, что попадает справа налево на стенку, справа от j-ой стенки
        end
        q = m + 1;
		hlr1(j) = tta*hlr0; %узнаем величину в абсолютных единицах этой лучистой энергии - слева от пластинки, падает слева направо
		hrl1(j) = ttb*hrl0; %справа от пластинки, падает справа налево
    end
    slr = slr + Tra(kost) * hlr1(kost) + hrl0*Ref(kost); qq=1; ri(qq) = slr; %то, что вышло из последней стенки вправо
	srl = srl + Tra(qq) * hrl1(qq) + Ref(qq) * hlr1(qq); ri(qq+1) = srl; %то, что вышло из первой стенки влево
	hlr1(qq) = hlr0; 
    hrl1(kost) = hrl0; 
    %for (k=0; k<ks; k++) cout << "hrl " << k << " = " << hrl1[k] << "\thlr " << k << " = " << hlr1[k] << endl;
	disp('Ras Vn Iz Sha k);
    switch (fl)
        case (0)
            vyma=hlr1;
        case (1)
            vyma=hrl1;
        case (2)
            vyma=ri;
    end
end