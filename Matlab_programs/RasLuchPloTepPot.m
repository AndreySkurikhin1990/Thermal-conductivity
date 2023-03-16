function [ vyma ] = RasLuchPloTepPot(kost, hrl1, hlr1, ao, Tra, Ref, Ab, sao, sislr, sisrl, Tsr,fl)
	k = 3; saot = sao;
	Erlr = zeros(1,kost); Elrr = zeros(1,kost); Erl = zeros(1,kost); Elr = zeros(1,kost); Etem = zeros(1,kost);
	Eres = zeros(1,kost); Eresl = zeros(1,kost); Eresr = zeros(1,kost);	
	for j = 1:kost
        Erlr(j) = 0.0; 
        Elrr(j) = 0.0; 
        Erl(j) = 0.0; 
        Elr(j) = 0.0; 
        Eres(j) = 0.0; 
        Eresl(j) = 0.0; 
        Eresr(j) = 0.0;
    end
	srl = 0.0; srl = 0.0; 
    tth = 0.0; tat = 0.0;
	for k = 1:kost-1
		srl = hrl1(k) * Tra(k); %слева от k-ой стенки
		slr = hlr1(k); 
        tta = slr + srl*Ref(k - 1); 
        Elr(k) = Elr(k) + tta; %в воздушном промежутке между (k-1)-ой и k-ой стенками - излучение, идущее слева направо
		ttb = srl + slr*Ref(k); 
        Erl(k) = Erl(k) + ttb; %излучение, идущее справа налево
		Eresl(k) = Eresl(k) + tta - ttb; %вправо - влево
		slr = sislr(k); 
        srl = sisrl(k) * Tra(k); 
        tta = slr + srl*Ref(k - 1); %слева - направо
		Elr(k) = Elr(k) + tta; 
        ttb = srl + slr*Ref(k); %справа - налево
		Erl(k) = Erl(k) + ttb; 
        Eresl(k) = Eresl(k) + tta - ttb; %результирующее излучение слева от k-ой стенки
    end
    for k = 1:kost 
		slr = hlr1(k) * Tra(k); %справа от k-ой стенки
		srl = hrl1(k + 1); 
        tta = srl + slr*Ref(k + 1); %влево
		Erlr(k) = Erlr(k) + tta; %справа от стенки, излучение идет справа налево
		ttb = srl*Ref(k) + slr; 
        Elrr(k) = Elrr(k) + ttb; %вправо
		Eresr(k) = Eresr(k) + ttb - tta; %результирующее излучение справа от k-ой стенки
		slr = sislr(k) * Tra(k); 
        srl = sisrl(k + 1); 
        tta = slr + srl*Ref(k); %вправо
		Elrr(k) = Elrr(k) + tta; 
        ttb = srl + slr*Ref(k + 1); %влево
		Erlr(k) = Erlr(k) + ttb; 
        Eresr(k) = Eresr(k) + tta - ttb;
    end
	Eres = Eresl - Eresr;
	ao(saot) = usrednen(Tsr, kost); saot=saot+1; %2
	ao(saot) = usrednen(Eresl, kost - 1); saot=saot+1; %3
	ao(saot) = usrednen(Eresr, kost - 1); saot=saot+1; %4
	ao(saot) = usrednen(Eres, kost - 1); saot=saot+1; %5
	ao(saot) = usrednen(Erl, kost - 1); saot=saot+1; %6
	ao(saot) = usrednen(Elr, kost - 1); saot=saot+1;
    ao(saot) = ao(saot) - ao(saot - 1); saot=saot+1; %7
	ao(saot) = usrednen(Elrr, kost - 1); saot=saot+1; %8
	ao(saot) = usrednen(Erlr, kost - 1); 
    ao(saot) = ao(saot - 1) - ao(saot); saot=saot+1; %9
	ao(saot) = (ao(saot - 1) + ao(saot - 3)) / 2.0; %10 
	for k=1:kost
        if (k==0) 
            tta=Eresl(k); 
        elseif (k==kost) 
            tta=Eresr(k); 
        else tta=(Eresr(k-1)+Eresl(k))/2e0; 
            Etem(k)=tta;
        end
    end
	switch (fl)
        case (0)
            vyma=ao;
        case (1)
            vyma=Eres;
        case (2)
            vyma=Etem;
    end
    end
        
    function ao = usrednen(usr, kost)
s = 0.0; k = 0.0; minz = 1e-15; hf = 1e0;
	for j=1:kost
        if (abs(usr(j))>minz) 
            s = s + usr(j); k = k + hf; 
        end
    end
	s = s / k; 
    ao = s; %средняя величина
    end