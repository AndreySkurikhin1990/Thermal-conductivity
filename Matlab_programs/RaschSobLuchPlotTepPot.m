function [ vyma ] = RaschSobLuchPlotTepPot(kost, prot, Ts, Tss, Tra, Ref, Ab, slr, srl, ao, sao, fl)
	q = 3;
	sig = 5.67e-8;
	sislr = zeros(1,kost); sisrl = zeros(1,kost);
	q = 1; 
    for j = 1:kost %узнаем, какая энергия собственного излучения дошла до j-ой стенки
		sislr(j) = 0.0; sisrl(j) = 0.0; tta = 0.0; ttb = 0.0;
		for k = 1:kost
			w0 = prot(q); 
            w1 = prot(q + 1); 
            w2 = prot(q + 2); 
            w3 = prot(q + 3); 
            wk = w2 + w3;
			tt = (Ts(k)^4e0)*w2 + (Tss(k)^4e0)*w3;
			if (abs(wk)>0.0)
				ttc = Ab(k) * w0*tt / wk; 
                tta = tta + ttc; %то, что попадает слева направо
				tt = Ab(k) * w1*tt / wk; 
                ttb = ttb + tt;
            end %то, что попадает справа налево
			q = q + 4;
        end
		sislr(j) = sig*tta; %собственное излучение стенок слева направо и справа налево
		sisrl(j) = sig*ttb;
    end
    q=1;
	slr = slr + Tra(kost) * sislr(kost);
	slr = slr + Ab(kost) * sig*(Tss(kost)^4e0);
	srl = srl + Tra(q) * sisrl(q);
	srl = srl + Ab(q) * sig*(Ts(q)^4e0);
	sislr(q) = 0.0; 
    sisrl(kost) = 0.0; %стенок вне многослойной стенки нет, излучение ниоткуда не падает, есть только внешнее излучение %излучение за крайней правой стенкой и за крайней левой
	p=1;
    ao(p) = slr; 
    p=p+1; 
    ao(p) = srl;
    switch (fl)
        case (0)
            vyma=ao;
        case (1)
            vyma=sislr;
        case (2)
            vyma=sisrl;
    end
    end