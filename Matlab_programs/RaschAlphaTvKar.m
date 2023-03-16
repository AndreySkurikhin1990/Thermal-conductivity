function [ vmua ] = RaschAlphaTvKar(por, vr, vss)
{
	int l=2, k=0, d=10, q=0, cfk=l, cei=0;
	long j=0, jk=100*100*100, pj=0, h=0;
	double x=0.0, y=0.0, p=0.0, xt=0.0, yp=0.0, hf=1e0, te0=273.15, tc=22.0+te0;
	double altc, *srra=NULL, pors=1e-6, lamtem=1448.0, dvpn = lamtem*pors / tc / 4e0;
	double *pn = new double[d], *alx = new double[d], *rcf = new double[cfk], e=1e-15;
	double *rapo, srp, marp, ce, cet, *pp=NULL, *legr=NULL, **pt=NULL, *prgr=NULL;
	pt=vybFunRasPorpoRazSha(por, vr, vss);
	k=0; rapo=pt[k]; k++; srra=pt[k]; k++; prgr=pt[k]; k++;
	legr=pt[k]; k++; pp=pt[k]; j=0; srp=pp[j]; if (pp) delete[]pp; k++;
	pp=pt[k]; marp=pp[j]; if (pp) delete[]pp; if (pt) delete[]pt;
	cet = ce; j = 0; while (cet > e) { cet = cet - hf; j++; } cei = j; //cout << "cem_srp = " << cei << "\tsrp = " << srp << "\t"; cout << endl; for (j=0; j<cei; j++) if (j<10) cout << "j = " << j << "\trpr = " << rapo[j] << "\tlegr = " << legr[j] << "\t"; cout << endl;
	x = 0.0; for (j = 0; j < cei; j++) if (j <= 31) x = x + srra[j] * rapo[j]; pors = x*por / srp;
	double pr = 0.0, rmf = 0.0, prf = 0.0, po = 0.0, *alsf = new double[cfk];
	for (j = 0; j < RAND_MAX; j++) rmf = rmf + hf; for (j = 0; j < jk; j++) po = po + hf;
	if ((pn) && (alx) && (rcf) && (alsf)) { for (j = 0; j < d; j++) { pn[j] = 0.0; alx[j] = 0.0; } 
	for (j = 0; j < cfk; j++) { rcf[j] = 0.0; alsf[j] = 0.0; } }
	else { cout << snm << endl; j = getchar(); exit(1); }
	long lt; unsigned int st; lt = time(NULL); st = (unsigned int)(lt - (lt % 2)) / 2; srand(st);
	for (k = 0; k < l; k++) {
		x = RasFracXeffSha60(k); //размер частицы
		rcf[k] = x*1e-6; y = x*pors; //размер поры
		for (j = 0; j < d; j++) pn[j] = 0.0;
		for (j = 0; j < jk; j++) {
			pj = rand(); prf = 0.0; for (h = 0; h < pj; h++) prf = prf + hf; pr = prf / rmf;
			yp = y*pr; xt = yp*(hf - pors) / pors;
			p = x / (xt + yp);
			pr = 0.0; for (h = 0; h < d; h++)	{ xt = pr + hf; if ((p >= pr) && (p < xt)) pn[h] = pn[h] + hf; pr = xt; }
		}
		pr = 0.0; for (j = 0; j < d; j++) { pn[j] = pn[j] / po; pr = pr + pn[j]; } //cout << "Summa = " << pr << endl; for (j=0; j<d; j++) cout << "pn ( " << j << " ) = " << pn[j] << "\t";
		for (j = 0; j < d; j++) pn[j] = pn[j] / pr;
		for (j = 2; j < d; j++) {
			p = 0.0; for (h = 0; h<j; h++) p = p + hf;
			yp = pors*x*1e-6 / (p - hf); xt = (hf - pors)*x*1e-6 / p; //cout << "x = " << x << "\txp = " << yp << "\txt = " << xt << "\t";
			if (yp>dvpn) {
				lamtem = F0_lamT(yp*tc); if (lamtem<0.0) lamtem = 0.0; if (lamtem>hf) lamtem = hf;
				alx[j] = BolTochRasAlpha(0, j, yp, xt, tc, ktpvosha, vtesha, etesha, cem, dmkv, snms, Tsh, Rsh, 2 * N, kttk, Ash, 1, Rash, Tash, Aash, Rtsh, Ttsh, Atsh)*lamtem;
			}
			else alx[j] = 0.0;
		}
		pt = RaschRTA(d, xt, 0.0, 0.0, yp, tc, 1, 0, 0, 0); pp = pt[0]; altc = pp[0]; delete[]pp; delete[]pt;
		alx[0] = 0.0; alx[1] = altc / (hf - pors); //for (j=0; j<d; j++) cout << "j = " << j << "\talx = " << alx[j] << "\t";
		for (j = 1; j<d; j++) {
			p = 0.0; for (h = 0; h <= j; h++) p = p + hf;
			yp = pors*x*1e-6 / (p - hf); if (j>1) xt = oprProcSoderpoPoris(rapo, legr, yp, cei); else xt = 1e0; //cout << "j = " << j << "\txt = " << xt << "\t";
			alsf[k] = alsf[k] + pn[j] * alx[j] * xt;
		}
	}
	x = 0.0; yp = 0.0; for (j = 0; j < cfk; j++) { x = x + alsf[j]; yp = yp + hf; } x = x / yp; //for (j=0; j<cfk; j++) cout << "j = " << j << "\tal_sr = " << alsf[j] << endl; 
	delete[]rapo; delete[]srra; delete[]rcf; delete[]pn; delete[]alx; 
	delete[]alsf; delete[]legr; delete[]prgr;
	x = (x / altc); printf("dkops = %0.10lf\n", x); //ослабление  ѕ за счет пористой структуры вермикулита
	return x;
}