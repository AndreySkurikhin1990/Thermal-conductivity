#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <dos.h>
#include <iostream> 
#include <time.h>
using namespace std;
class KoePog { public: double alp, tem; KoePog *nex; };
//-------
double **RaschRTAkvi(int, double, double, double, int, int, double, int, double, int, int, int);
double **vybFunRasPorpoRazSha(double, int, int);
double F0_lamT(double);
double oprProcSoderpoPoris(double *, double *, double, int);
double *zadrkt(int, int, double, int, double, double, int, int, double, int, int, double *);
double *EffectTols(double *, double *, double *, double, double, int);
double RasFracXeffVer60(int);
double RasFracXeffVer60_100(int);
double RasFracXeffVer100_150(int);
double RasFracXeffVer150_200(int);
double RaschAlphaTvKarSha();
double **chaRTA(int, double *, double *, double *, double *, double *, double *, int);
double **izmRTA(double *, int, double *, double *, double *, double *, double *, double *, int, int);
double **RaschVneshIzluch(double *, double *, double *, double *, double *, int, double, double, char *);
double opredKTPTKToch(double *, double *, double, int);
double *opredKoefOtr(double *, double *, int, double, int, int, double *, double *, int);
double ReflSredSha(double);
double ReflSredVer(double);
double ReflSreditom(double, int);
double ReflSredkvi(double, int);
double *kopo(double, double, int, double, int, double, double, int, double *, char *, char *, char *, double, double, double, int, int);
double *podschchiele(int, int, int, double *, double *);
double *KoefPoglRosselNac(double *, int, int, double, double, double, double *, double *, int, double, double, double *, double *, int, int, int);
double **rasPorpoRazVer(double, int, int, int, int, int);
//-------
double RaschAlphaTvKar(double por, int vr, int vss, char *snm)
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
	x = (x / altc); printf("dkops = %0.10lf\n", x); //ослабление КП за счет пористой структуры вермикулита
	return x;
}
double RasFracXeffSha60(int v)
{
	int l = 2, k; double *mkbr = new double[l], *mv = new double[l], *tol = new double[l];
	k = 0; mkbr[k] = 238.57; k++; mkbr[k] = 227.973; k++;
	k = 0; mv[k] = 1.706; k++; mv[k] = 1.1; k++;
	k = 0; tol[k] = 0.7; k++; tol[k] = 0.68; k++;
	double rokbr = 2.75, rov = 2.69, *xv = EffectTols(mkbr, mv, tol, rov, rokbr, l), xsr = 0.0;
	if (!v) xsr = xv[0]; else if (v == 1) xsr = xv[1];
	delete[]mkbr; delete[]mv; delete[]tol; delete[]xv; return xsr;
}
double BolTochRasAlpha(int vyte, int kost, double hvp, double htk, double tc, double *ktpvo, double *vte, double *ete, int ce, int dmk, char *snm, double *Tra, double *Ref, int nao, double *ktptk, double *Ab, int vyv, double *Refa, double *Traa, double *Aba, double *Reft, double *Trat, double *Abt)
{
	int j, l = nao, k; double *prot = NULL, hlr0 = 1e0, hrl0 = 0.0, alks, ra = fabs(hlr0 - hrl0), tol = 1e0;
	double *ao = new double[l], *Ts = NULL, *Tss = NULL, *Tsr = NULL, *Tna = NULL, **mu, *po, *sislr = NULL, *sisrl = NULL;
	double *reiz, *hrl1 = NULL, *hlr1 = NULL, d = 0.0, slr = 0.0, srl = 0.0, hf = 1e0; for (j = 0; j < kost; j++) d = d + hf;
	prot = zadrkt(1, kost, d, 0, htk, hvp, 0, 0, tc, 1, 0, prot); //vyv - выбор вещества: 0 - шамот, 1 - вермикулит
	if (vyv>=4) { cout << "Net takogo veschestva!" << endl; k = getchar(); exit(1); }
	mu = RaschRTA(kost, htk, 0.0, 0.0, 1, 0, hvp, 0, tc, 0, 0, 1);
	k=0; 
	Refa = mu[k]; k++; Traa = mu[k]; k++; Aba = mu[k];  k++; Ref = mu[k];  k++;
	Tra = mu[k];  k++; Ab = mu[k];   k++; Reft = mu[k]; k++; Trat = mu[k]; k++; Abt = mu[k]; delete[]mu;
	mu = opredTempLPStenSha(Ts, Tss, Tsr, Tna, tc, 2 * kost, ktpvo, vte, ete, ktptk, ce, dmk, htk, hvp, 0.0, kost, snm);
	k = 0; Ts = mu[k]; k++; Tss = mu[k]; k++; Tsr = mu[k]; k++; Tna = mu[k]; delete[]mu;
	mu = RaschVneshIzluchSha(Tra, Ref, prot, hlr1, hrl1, kost, hlr0, hrl0, snm);
	k = 0; hlr1 = mu[k]; k++; hrl1 = mu[k]; k++; reiz = mu[k]; delete[]mu;
	k = 0; slr = reiz[k]; k++; srl = reiz[k];
	delete[]Tsr; delete[]Tna; delete[]hlr1; delete[]hrl1;
	mu = RaschSobLuchPlotTepPot(kost, prot, Ts, Tss, Tra, Ref, Ab, slr, srl, ao, 0, snm);
	k = 0; ao = mu[k]; k++; sislr = mu[k]; k++; sisrl = mu[k]; delete[]mu;
	k = 0; slr = ao[k]; k++; srl = ao[k];
	delete[]ao; delete[]prot; delete[]Ts; delete[]Tss; delete[]reiz;
	slr = fabs(srl - slr); ra = fabs(hlr0 - hrl0); ra = slr / ra;
	tol = d*htk + (d - hf)*hvp;
	alks = -log(ra) / tol;
	return alks;
}
double oprProcSoderpoPoris(double *raspor, double *legr, double rp, int ce)
{
	double *prgr = new double[ce], x = 0.0, hf = 1e-6; int k;
	if (prgr) { for (k = 0; k < (ce - 1); k++) prgr[k] = legr[k + 1]; prgr[ce - 1] = prgr[ce - 2] + hf; }
	else { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < ce; k++) if ((legr[k]<rp) && (prgr[k]>rp)) { x = raspor[k] * (rp - legr[k]) / hf; break; } //cout << "x = " << x << endl;
	for (k = 0; k<ce; k++) if (rp>prgr[k]) x = x + raspor[k]; if (x < 0.0) x = 0.0; //cout << "ra_po = " << rp << "\tx = " << x << "\t"; cout << endl;
	delete[]prgr; return x;
}
double **RaschRTA(int kost, double htm, double kta, double ktb, int izm, int vyte, double hvo, int v, double ti, int c, int w, int u, double tocras, char *snm, double *Ra, double *Rt, double *R, double *Aa, double *At, double *A, double *Ta, double *Tt, double *T, double *alpha, int vybves, int cem, double *mko, double *ete, double dkoal, double dkosp)
{
	int k=kost, j=0, i=0; double *te=NULL, e=tocras, **mauk=NULL, hf=1e0;
	if (!v) {
		if ((vybves>=0) && (vybves<4)) {
		alpha=kopo(kta, ktb, 2, ti, kost, htm, hvo, 0, alpha, snm, sfno, skpt, dkoal, dkosp, hk, cem, cektp);
		te = kopo(kta, ktb, 0, ti, kost, htm, hvo, 0, alpha, snm, sfno, skpt, dkoal, dkosp, hk, cem, cektp); }
		if (!c) { k = kost;
if (Ra) delete[]Ra; Ra=new double[k];
if (R) delete[]R; R=new double[k];
if (Rt) delete[]Rt; Rt=new double[k];
Ra=opredKoefOtr(te, Ra, k, tocras, 1, cem, mko, ete, 0);
for (j = 0; j < k; j++) { Rt[j] = Ra[j]; R[j] = Ra[j]; }
if (Ta) delete[]Ta; Ta=new double[k];
if (Tt) delete[]Tt; Tt=new double[k];
if (T) delete[]T; T=new double[k];
if (Aa) delete[]Aa; Aa=new double[k];
if (At) delete[]Atver; At=new double[k];
if (A) delete[]A; A=new double[k];
	for (j = 0; j < k; j++) {
	Ta[j]=-alphaVer[j]*htm; Ta[j]=exp(Ta[j]); Aa[j]=hf-Taver[j]-Raver[j];
	if (Aaver[j] < 0.0) { Aaver[j] = alphaVer[j] * htm; Taver[j] = 1e0 - Raver[j] - Aaver[j]; }
				Atver[j] = Aaver[j]; Ttver[j] = Taver[j]; Aver[j] = Aaver[j]; Tver[j] = Taver[j];
			}
			k = kost; mauk = izmRTAVer(te, k, 1, Rver, Tver, Aver, Raver, Taver, Aaver, 0);
			i=6; mauk[i] = Rtver; i++; mauk[i] = Ttver; i++; mauk[i] = Atver;
		}
		else if (c == 1) { i=1; mauk = new double*[i]; 
		if (!mauk) { cout << snm << endl; j = getchar(); exit(1); } i=0; mauk[i] = te; }
	}
	else if (v == 1) {
		te = kopoVer(0.0, 0.0, 1, ti, kost, c, htm, hvo, w, alpha); i=1;
		mauk=new double*[i]; if (!mauk) { cout << snm << endl; j = getchar(); exit(1); } i=0; mauk[i] = te;
	}
	if (!u) mauk[0] = te; 
	return mauk;
}
double *kopo(double kta, double ktb, int vyb, double tl, int kost, double htk, double hvp, int w, double *al, char *snm, char *sfno, char *skpt, double dkoal, double dkosp, double hk, int cem, int cektp)
{
	int k=0, q=100, j=0, r=0, f=0;
	double te0=273.15, t0=2e1, te=t0+te0, dte=1e0, t=0.0, *p=NULL, a1=0.0, a2=0.0, t1=0.0, t2=0.0;
	char *s=new char[q]; KoePog *kp=new KoePog, *ne=NULL, *roo=NULL, *pre=NULL;
	if ((!kp) || (!s)) { cout << snm << endl; j=getchar(); exit(1); } for (j=0; j < q; j++) s[j] = '\0';
	ifstream fin; fin.open(skpt); if (!fin.is_open()) { cout << sfno << endl; j=getchar(); exit(1); }
	roo = kp; k = 0; while (!fin.eof()) {
		fin.getline(s, q, '\n');
		ne = new KoePog; if (!ne) { cout << snm << endl; j = getchar(); exit(1); }
		kp->alp = atof(s)*dkoal*dkosp; kp->tem = te; 
		kp->nex = ne; pre = kp; kp = ne; k++; te = te + dte; kp->nex = NULL;
	}
	if (ne) delete[]ne; kp=NULL; pre->nex=kp; if (s) delete[]s; r=k;	fin.close(); //cout << "r = " << r << "\tvyb = " << vyb << endl;
	if ((!vyb) || (vyb == 2)) { //поиск массива коэффициентов поглощения при температурах стенок
		f = 2; double *xi = new double[kost], *teks = new double[kost], *koe = new double[f], knat=0.0;
		if ((!teks) || (!koe) || (!xi)) { cout << snm << endl; j = getchar(); exit(1); }
		teks = oprRasTemNach(cem, cektp, teks, koe, kost, xi, htk, hvp, hk, kta, ktb, w);
		for (k = 0; k < kost; k++) {
			kp = roo->nex; pre = roo; j = 0; t = teks[k];
			while ((kp) && (j<r)) {
				if (kp->tem>t) {
					a1 = pre->alp; a2 = kp->alp; 
					t1 = pre->tem; t2 = kp->tem; 
					knat = (a2 - a1) / (t2 - t1); 
					al[k] = a1 + knat*(t - t1); break; }
				else { pre = kp; kp = kp->nex; } j++; } } 
		if (xi) delete[]xi; if (koe) delete[]koe;
		if (vyb == 2) { p = al; delete[]teks; }
		else p = teks;
	}
	if (vyb == 1) { //поиск конкретного коэффициента поглощения при данной температуре tl
		kp = roo->nex; j = 0; pre = roo; t = 0.0;
		while ((kp) && (j<r)) {
			if (kp->tem>tl)
			{
				a1 = pre->alp; a2 = kp->alp; 
				t1 = pre->tem; t2 = kp->tem; 
				t = a1 + (a2 - a1)*(tl - t1) / (t2 - t1); break; }
			else { pre = kp; kp = kp->nex; } j++; } j = 1; p = new double[j]; j=0;
		if (p) p[j] = t; else { cout << snm << endl; j = getchar(); exit(1); } }
	kp = roo; while (kp) { ne = kp->nex; delete kp; kp = ne; } //удаление списка
	return p;
}
double *oprRasTemNach(int ce, int cee, double *teks, double *koe, int kost, double *xi, double htch, double hvozd, double hko, double a, double b, int w, double tocras, int cemunau, int vt, double *thol, double *tgor, double tol)
{
	int k=0, j=0, f=cemunau; 
	double e=tocras, hkx=0.0, ht=0.0, tsre=0.0, **mu=new double*[f]; 
	if (!w) {
		hkx=hko; k=0; koe[k]=(thol[vt]-tgor[vt])/tol; k++; koe[k]=tgor[vt];
	}
	else if (w==1) { hkx=tol/2e0; k=0; koe[k]=a; k++; koe[k]=b; }
	k=0; xi[k]=hkx+htch/2e0; ht=hvozd+htch; for (k=1; k<kost; k++) xi[k]=xi[k-1]+ht; //массив середин каждой из стенок по толщине
	j=0; for (k=0; k<kost; k++) teks[k]=koe[j]*xi[k]+koe[j+1];
	k=0; j=kost-1; tsre=(teks[k]+teks[j])/2e0; //qobsha = opredKTPTKTochSha(qob, etesha, tsre, ce); //for (k=0; k<kost; k++) cout << "k = " << k << "\txk = " << xi[k] << "\ttex = " << teks[k] << endl;
	if (mu) delete[]mu;
	return teks;
}
double *zadrkt(int zf, int kost, double d, int vyte, double htk, double hvo, int prod, int vy, double tc, int c, int u, double *rta, char *snm)
{
	int j; double **mu = RaschRTA(kost, htk, 0.0, 0.0, 1, vyte, hvo, 0, tc, 1, 0, 0), *tt = mu[0], *te = new double[kost];
	if (!te) { cout << snm << endl; j = getchar(); exit(1); }
	else for (j = 0; j < kost; j++) te[j] = tt[j];
	if (!u) {
		int k = 0, kst, q = 0, m = 4, b;
		double *prs = new double[m*kost*kost], *pr = new double[m], Er;
		if ((!prs) || (!pr)) { cout << snm << endl; k = getchar(); exit(1); }
		for (k = 0; k < (m*kost*kost); k++) prs[k] = 0.0;
		q = 0; for (kst = 1; kst <= kost; kst++)
		{
			for (b = 0; b < m; b++) pr[b] = 0.0;
			for (k = 1; k <= kost; k++) {
				pr = izstN(k, kst, m, kost);
				for (b = 0; b < m; b++) prs[q + b] = pr[b]; q = q + m;
			}
		}
		if (prod == 1) Er = opredLuchSost(prs, te, q, zf, vyte);
		delete[]pr; delete[]te; delete[]tt; delete[]mu; return prs;
	}
	else { delete[]tt; delete[]mu; return te; }
}
double *izstN(int izst, int kst, int l, int ocs, double *R, double *T, int N)
{
	int k; double *o;
	if (abs(izst - kst) <= N) o = podschchiele(izst, kst, ocs, R, T);
	else {
		o = new double[l]; if (!o) { cout << "No memory!"; k = getchar(); exit(1); }
		for (k = 0; k < l; k++) o[k] = 0.0;
	}
	return o;
}
double *FuncRaschIzl(double ta, double tb, double d, double tc, int m, double *prot, double de, double *ao, double rez, double rap, double b, double c, int vybo, double *qobve, double *eteve, int cemve, int vyve, char *snm, double dpct, int vyte, double *ektpve, double y0ve, double *ktrve, double *kttkve, int kost, double *ktpvove, double *vteve, int dmkvove, double *ecktpve, double htch, double hvoz, double *Refe, double *Trae, double *Abse, double *Reft, double *Trat, double *Abst, double *txve, double *qxve, double *dtxve, int sctve, double tocrasve, double poristost, int dmao, double *stchsrve, double *Tpctve, double b0, double c0, double mk, char *navyfa, int fla)
{
	int j, zf = 1, cels = 6, k, nvm = 11, ksu = 2 * kost, nfm = 28;
	double qrr = 0.0, qr = 0.0, gtm = 0.0, dtma = 0.0, t = 0.0, hlr0 = 0.0, hrl0 = 0.0, rhlr0 = 0.0, slr = 0.0, srl = 0.0, ra10, **mu, ra9, ra, ra7;
	double *Ts = NULL, *Tss = NULL, *sislr = NULL, *sisrl = NULL, *Tna = NULL, *reiz = NULL, *hrl1 = NULL, *hlr1 = NULL, *Tsr = NULL, qis;
	double *fm = new double[nfm], ktptk, ecktp, ktpr, dp = (1e0 - dpct), qs, *tevy = NULL, r, tesr, dote, p, *vm = new double[nvm];
	double qobvee, gg, ep = tocrasve, rr, ra0, tcc=(ta+tb)/2e0, *Evsv=NULL;
	if (!vm) { cout << snm << endl; k = getchar(); exit(1); }
	gtm = LuchKTPChudnovsky(Abse, tc, kost, hvoz); //по Чудновскому
	qobvee = opredKTPTKTochSha(qobve, eteve, tc, cemve); if (qobvee <= 0.0) qobvee = opredKTPTKTochSha(qobve, eteve, eteve[vyte], cemve);
	qs = dp*qobvee;
	p = opredKTPTKTochSha(ektpve, eteve, tc, cemve);
	dtma = fabs(tb - ta);
	gtm = dtma / y0ve;
	ktpr = opredKTPTKTochSha(ktrve, eteve, tc, cemve);
	qr = ktpr*gtm; qrr = qr;
	ktptk = opredKTPTKTochSha(kttkve, eteve, tc, cemve);
	ecktp = opredKTPTKTochSha(ecktpve, eteve, tc, cemve);
	if (!m) { t = (dpct*ecktp*b + (1e0 - dp*b)*ktptk); if (fabs(t) > ep) hlr0 = qs*ktpr / t; else hlr0 = 0.0; }
	else hlr0 = rez; hrl0 = 0.0; rhlr0 = hlr0 - hrl0;
	if (qs<rhlr0) { rhlr0 = qs; hlr0 = rhlr0 + hrl0; }
	mu = opredTempLPStenSha(Ts, Tss, Tsr, Tna, tc, ksu, ktpvove, vteve, eteve, kttkve, cemve, dmkvove, htch, hvoz, qs - rhlr0, kost, snm);
	k = 0; Ts = mu[k]; k++; Tss = mu[k]; k++; Tsr = mu[k]; k++; Tna = mu[k]; delete[]mu;
	mu = RaschVneshIzluchSha(Trae, Refe, prot, hlr1, hrl1, ks, hlr0, hrl0, snms);
	k = 0; hlr1 = mu[k]; k++; hrl1 = mu[k]; k++; reiz = mu[k]; slr = reiz[0]; srl = reiz[1]; ra0 = slr - srl; delete[]mu;
	tevy = reshnewtrafs(prot, Tna, hlr1, hrl1, zf, Trae, Abse, kost, kttkve, eteve, cemve, ktpvove, vteve, dmkvove, qs, txve, qxve, sctve, htch, hvoz, tocrasve, ektpve, cels + 2 * kost, tb, ta, y0ve, snms);
	tesr = tevy[2]; txve[sctve] = tesr; qxve[sctve] = tevy[1]; dtxve[sctve] = tevy[0]; gtm = tevy[4];
	if ((m>0) && (fabs(tevy[4]) > 0.0)) r = tevy[1] * b*dp / tevy[4] + opredKTPTKTochSha(kttkve, eteve, tesr, cemve)*(1e0 - b*dp); else r = 0.0; dote = tevy[5]; delete[]tevy;
	if (r<0.0) r = 0.0; if (fabs(p) > ep) gg = fabs((r - p) / p)*1e2; else gg = mk; if (mk < gg) { b0 = b; c0 = c; mk = gg; }
	for (j = 0; j < dmao; j++) ao[j] = 0.0;
	mu = RaschSobLuchPlotTepPot(kost, prot, Ts, Tss, Trae, Refe, Abse, slr, srl, ao, 0, snm); 
	k = 0; ao = mu[k]; k++; sislr = mu[k]; k++; sisrl = mu[k]; delete[]mu;
	mu = RasLuchPloTepPot(kost, hrl1, hlr1, ao, Trae, Refe, Abse, 2, sislr, sisrl, Tsr, snm); 
	k = 0; ao = mu[k]; k++; delete[]reiz; reiz = mu[k]; k++; Evsv=mu[k]; delete[]mu; 
	k = 0; slr = ao[k]; k++; srl = ao[k]; ra = slr - srl;
	qr = RasIzlSerStenNac(Refe, Trae, Abse, prot, Reft, Trat, Abst, tc, hlr0, hrl0, b, vyve, poristost, dpct, htch, hvoz, kost, eteve, kttkve, cemve, vyte, qobve, ktpvove, vteve, snm, dmao, stchsrve, Tpctve, dmkvove); qis = qr; //0 - шамот, 1 - вермикулит, 2 - ИТОМ
	qr = (qr*(d - 1e0) + ra)*dp*b / d;
	k = 7; ra7 = (ra + ao[k] * (d - 1e0)) / d; ra7 = ra7*dp*c*b; //ra7=qrr+ra7; 
	k = 9; ra9 = (ra + ao[k] * (d - 1e0)) / d; ra9 = ra9*dp*c*b; //ra9=ra9+qrr; //c - поправка на разрешающие угловые коэффициенты из-за небесконечности стенок, b - их аспектное соотношение (поправка на их некубичность)
	k = 10; ra10 = (ra + ao[k] * (d - 1e0)) / d; ra10 = ra10*dp*c*b; //ra10=qrr+ra10; 
	rez = ra; t = opredKTPTKTochSha(kttkve, eteve, tc, cemve); 
	if (fabs(p) > ep) gg = fabs((r - p) / p)*1e2; else gg = mk; if (mk>gg) { b0 = b; c0 = c; mk = gg; }
	if (fabs(qr) > ep) rr = ra9/qr; else rr = 0.0;
	k = 0; fm[k] = rr; k++; fm[k] = dote;     k++; fm[k] = qr;      k++; fm[k] = qrr;      k++; fm[k] = r;         k++; fm[k] = gg;    k++;
	fm[k] = b;         k++; fm[k] = c;        k++; fm[k] = qis;     k++; fm[k] = ao[5];    k++; fm[k] = ao[7];     k++; fm[k] = ao[9]; k++;
	fm[k] = ao[10];    k++; fm[k] = qobvee;   k++; fm[k] = rap;     k++; fm[k] = t;        k++; fm[k] = ra9;       k++; fm[k] = ra10;  k++;
	fm[k] = ta;        k++; fm[k] = tb;       k++; fm[k] = ra0;	    k++; fm[k] = ra9 - qr; k++; fm[k] = ra10 - qr; k++; fm[k] = ra7;   k++; 
	fm[k] = ra7 - qr;  k++; fm[k] = ra;       k++; fm[k] = hksha;   k++; fm[k] = tcc;	   k++;
	if (fla) { vyvodfile(fm, nfm, 2, 0.0, navyfa); //cout << "tc = " << tc << "\ttNR = " << dote << "\tqr = " << qr << "\tqrr = " << qrr; cout << "\tktp = " << r << "\topktp = " << gg << "\tb = " << b << "\tc = " << c << "\tqss = " << qis; cout << "\tao5 = " << ao[5] << "\tao7 = " << ao[7] << "\tao9 = " << ao[9] << "\tao10 = " << ao[10] << "\tqob = " << qobvee; cout << "\trap = " << rap << "\tktp_tk = " << t << "\tra9 = " << ra9 << "\tra10 = " << ra10; cout << "\tta = " << ta << "\ttb = " << tb << "\tslr-srl = " << slr - srl << "\tra9-qr = " << (ra9 - qr) << "\tra10-qr = " << ra10 - qr << endl; 
	} delete[]fm; delete[]sislr; delete[]sisrl; delete[]hlr1; delete[]hrl1; 
	delete[]reiz; delete[]Tsr; delete[]Ts; delete[]Tss; delete[]Tna; delete []Evsv;
	k = 0; vm[k] = ra; k++; vm[k] = qr; k++; vm[k] = qrr; k++; vm[k] = ra9 - qr; k++; vm[k] = dote; k++; 
	vm[k] = b0; k++; vm[k] = c0; k++; vm[k] = mk; k++; vm[k] = rr; k++; vm[k] = gg; k++; vm[k] = r; return vm;
}
double oprEffDoliTepPeren(double ko, double d, double por, double ektp, int kost, double qobsh, double e, double *ktpvo, double *vte, int cem, char *snm, double *ete, double kttk, int dmkvo, double *ete, double tsr)
{
	int k=0, j=0, f=0, ksu=2*kost; double hvozd=ko, htch=0.0, hf=1e0, *ecktp=NULL;
	double *ecktp=new double[kost], *tepo=new double[ksu], r, p, sa, sb, sc, t, fa, fb, fc;
	if ((!tepo) || (!ecktp)) { cout << snm << endl; k=getchar(); exit(1); }
	for (j=0; j<ksu; j++) tepo[j]=0.0;
	hvozd=ko; htch=hvozd*(hf-por)/por; //cout << "hvoz = " << hvozd << "\tht = " << htch << endl;
	tepo=opredTempStenShaFragm(tepo, ksu, ktpvo, vte, ete, kttk, cem, dmkvo, htch, hvozd, qobsh, ete[j], -hf); //в середине слоя
	r=0.0; for (k=0; k<ksu; k++) for (f=0; f<ksu; f++) { p=tepo[k]-tepo[f]; if (p>r) r=p; } //максимальная разность температур - на краях стенок
	p = hvozd*(d - hf) + htch*d; r = r / p; t = qobsh / r; ecktp=t; //средний КТП многослойной стенки
	f = 10000; sa = 0.0; sb = hf; k = 0;
		do {
			sc = (sa + sb) / 2e0;
			fa = kttk * sa + ecktp * (hf - sa) - ektp; //эффективные КТП многослойной стенки и  перемычки должны сравняться
			fb = kttk * sb + ecktp * (hf - sb) - ektp; //чтобы найти относительные доли площадей сечения переноса общей ПТП
			fc = kttk * sc + ecktp * (hf - sc) - ektp;
			if ((fc*fb>0.0) && (fa*fc<0.0)) sb = sc;
			if ((fc*fa>0.0) && (fb*fc<0.0)) sa = sc;
			r = fabs(sa - sb); k++;	} 
		while ((r>e) && (k<f));
	if (tepo) delete[]tepo; if (ecktp) delete[]ecktp;
	return sc;
}
double *KorrZnachVozdPros(double hps, double ksf, double por, int vy, double e)
{ 	int j=0, k=1000, q=3; double pa = 1e-2, pb = 1e0, *po=NULL, pc=0.0, hf = 1e0, ra = hf;
	double fa, fb, fc, ta, tb, tc, tca, tcb, tcc, ka = hps*pa / por, kb = pb*hps / por, kc, *pkc=new double[q]; //cout << "hps = " << hps << "\tksuf = " << ksf << "\tpor = " << por << endl;
	j=0; pa=1e-2; pb = 1e0; ra = hf; ka = hps*pa / por; kb = pb*hps / por;
	while ((ra>e) && (j < k)) { //подтягиваем пористость к значению, которое задали изначально, во время подстройки ЭКТП
		pc = (pa + pb) / 2e0; //эффективная пористость
		kc = hps*pc / por;
		tc = oprEffDoliTepPeren(kc, ksf, pc); //доля площади, переносимая чистой теплопроводностью твердого каркаса
		tcc = kc*(hf - pc) / pc; //толщина твердой части слоя
		fc = (hf - tc)*kc / (kc + tcc) - por;
		ka = hps*pa / por;
		ta = oprEffDoliTepPeren(ka, ksf, pa); //определяем долю площади сечения перемычки
		tca = ka*(hf - pa) / pa;
		fa = (hf - ta)*ka / (ka + tca) - por;
		kb = hps*pb / por;
		tb = oprEffDoliTepPeren(kb, ksf, pb); //через перемычку тепло распространяется чистой теплопроводностью
		tcb = kb*(hf - pb) / pb;
		fb = (hf - tb)*kb / (kb + tcb) - por;
		if (((fc*fb) > 0.0) && ((fa*fc)<0.0)) pb = pc; if ((fc*fa>0.0) && (fb*fc < 0.0)) pa = pc;
		j++; ra = fabs(pa - pb); }
	cout << "Tolschina vozdushnoy prosloyki = " << kc << "\tDolya Ploschadi CTP = " << tc << "\tEffectivnaya poristost = " << pc << endl;
	j=0; pkc[j]=kc; j++; pkc[j]=tc; j++; //скорректированное значение толщины воздушной прослойки (размер поры), когда ввели перемычку
	pkc[j]=KorrZnachVozdProsUtoch(hps, ksf, por, tc, e, tem)
	return pkc;
} //доля площади, через которую происходит перенос тепла чистой теплопроводностью
double KorrZnachVozdProsUtoch(double hps, double ksf, double por, double dpctv, double e, double tj)
{
	int j = 0, k = 1000, q=0; double pa=1e-3, hf=1e0, pb=hf, pc=0.0, ra=hf;
	double fa, fb, fc, tca, tcb, tcc, ka = hps*pa / por, kb = pb*hps / por, kc; //cout << "hps = " << hps << "\tksuf = " << ksf << "\tpor = " << por << endl;
	while ((ra>e) && (j < k)) { //подтягиваем пористость к значению, которое задали изначально, во время подстройки ЭКТП
		pc = (pa + pb) / 2e0;
		kc = hps*pc / por; //один период стенки
		tcc = kc*(hf-pc) / pc; //толщина твердой части
		fc = (hf-dpctv)*kc / (kc + tcc) - por;
		ka = hps*pa / por;
		tca = ka*(hf-pa) / pa;
		fa = (hf-dpctv)*ka / (ka + tca) - por;
		kb = hps*pb / por;
		tcb = kb*(hf-pb) / pb;
		fb = (hf-dpctv)*kb / (kb + tcb) - por;
		if ((fc*fb > 0.0) && (fa*fc<0.0)) pb = pc;
		if ((fc*fa>0.0) && (fb*fc < 0.0)) pa = pc;
		j++; ra = fabs(pa - pb);
	} //скорректированное значение толщины воздушной прослойки (размер поры), когда ввели перемычку
	return kc;
}
double opredLuchSost(double *prot, double *Tsr, int le, int zf, int vyte, int N, int kost, double bn, double bk, double cn, double ck, char *snm, double qob, int vt, int vtn, int sctx, int nxt, int nnxt, double y0, double *dkosct, double *dkoscm, int dkoscl, double dkosp, double dkoal, double *tkusc, double *kusc, int dmkoosc, double *thol, double *tgor, double *ete, int cem, double smgo, double ssio, double salo)
{
	int k=0, j=0, q=0, m=0, l=3, mma=3000, mm=0, mt=0, vy=9, idr=1;
	double **mauk=NULL, *pt=NULL, s=0.0, p=0.0, rez=0.0, rap=0.0, *temras=NULL;
	double *ao=new double[2*N], *po=NULL, d=0.0, hb=5e-1, hb0=hb, hf=1e0, *ktpr=NULL;
	d=0.0; for (j = 0; j < kost; j++) d = d + hf;
	double de=hf, ta=0.0, tb=ta, tc=ta, fa=ta, fb=tb, fc=tc;
	double ba=bn, bb=bk, bc=bn, xk=0.0, dxk=1e-3, co=0.0, bo=co;
	double ca=cn, cb=ck, cc=cb, hc=1e-1, hc0=hc, tt=0.0, raz=1e2;
	double mk=raz, mkm=mk, c00=0.0, b0=bc, c0=cc, kdc=5e0, kdb=5e0, tmin=0.0, tmax=0.0;
	double *R=NULL, *T=NULL, *A=NULL, *Ra=NULL, *Ta=NULL, *Aa=NULL, *Rt=NULL, *Tt=NULL, *At=NULL;
	if (!ao) { cout << snm << endl; k = getchar(); exit(1); }
	for (k=0; k<(2*N); k++) ao[k]=0.0; //cout << "qo = " << qob << "\n";
	if ((vt==vtn) && (sctx<=2) && (!temras)) {
		j=nxt-nnxt; temras=new double[j]; if (!temras) { cout << snm << endl; k = getchar(); exit(1); }
		ta=(thol[vt]-tgor[vt])/y0; tb=tgor[vt];
		xk=0.0; for (k=0; k<j; k++) { temras[k]=tb+ta*xk; xk=xk+dxk; }	}
	if ((vt==vtn) && (sctx<=2) && (!ktpr)) {
		ktpr=KoefPoglRosselNac(ete, idr, cem, smgo, ssio, salo, dkosct, dkoscm, dkoscl, dkosp, dkoal, tkusc, kusc, dmkoosc, idr, 1);
		for (k = 0; k < cem; k++) cout << "ktr = " << ktpr[k] << "\ttem = " << ete[k] << endl; } 
	l=vt; tmin=thol[l]; tmax=tgor[l];
	ta=tmin; tb=tmax; tc=temras[sctx];
	mauk=RaschRTA(kost, h, 0.0, 0.0, 1, vt, hvo, 0, tc, 0, 0, 1);
	k=0; Ra=mauk[k]; k++; Ta=mauk[k]; k++; Aa=mauk[k]; k++; R=mauk[k]; k++; T=mauk[k]; k++;
	A=mauk[k]; k++; Rt=mauk[k]; k++; Tt=mauk[k]; k++; At=mauk[k]; delete[]mauk;
	ca=cn; cb=ck; m=0; q=0;
	while ((m<mma) && (raz>1e-3) && (!q)) {
		cc=(ca+cb)/2e0;
		pt=FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, cc, vy, qob, ete, cem, 0, snm, dpct, vt, ektp, y0, ktpr, kttk, kost, ktpvo, vte, dmkvo, ecktp, h, hvo, R, T, A, Rt, Tt, At, tx, qx, dtx, sctx, tocras, por, 2*N, stchsr, Tpct, b0, c0, mk, svfdv, 0);
		k=3; fc=pt[k]; delete[]pt;
		if (fabs(fc)<hf) { c0=cc; break; }
		pt=FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, ca, vy, qob, ete, cem, 0, snm, dpct, vt, ektp, y0, ktpr, kttk, kost, ktpvo, vte, dmkvo, ecktp, h, hvo, R, T, A, Rt, Tt, At, tx, qx, dtx, sctx, tocras, por, 2*N, stchsr, Tpct, b0, c0, mk, svfdv, 0);
		k=3; fa=pt[k]; delete[]pt;
		if (fabs(fa)<hf) { c0=ca; break; }
		pt=FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, cb, vy, qob, ete, cem, 0, snm, dpct, vt, ektp, y0, ktpr, kttk, kost, ktpvo, vte, dmkvo, ecktp, h, hvo, R, T, A, Rt, Tt, At, tx, qx, dtx, sctx, tocras, por, 2*N, stchsr, Tpct, b0, c0, mk, svfdv, 0);
		k=3; fb=pt[k]; delete[]pt;
		if (fabs(fb)<hf) { c0 = cb; break; }
		if ((fa*fc)<0.0) { cb = cc; q = 0; }
		if ((fb*fc)<0.0) { ca = cc; q = 0; }
		if (((fb*fc)<0.0) && ((fa*fc)<0.0)) { cb = cc; q = 0; }
		if (((fb*fc)>0.0) && ((fa*fc)>0.0)) q = 1;
		m++; raz=fabs(cb-ca);
	}
	c0=cc; c00=c0; co=c00; cout << "c0 = " << c0 << "\ttb = " << tb << "\ttc = " << tc; 
	m=0; ba=bn; bb=bk; l=vt; ta=thol[l]; tb=tgor[l]; tc=temras[sctx];
	while ((m<mma) && ((bb-ba)>1e-1)) {
		bc=(ba+bb)/2e0; b0=bc; c0=c00; cc=c0; 
		pt=FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, c0, vy, qob, ete, cem, 0, snm, dpct, vt, ektp, y0, ktpr, kttk, kost, ktpvo, vte, dmkvo, ecktp, h, hvo, R, T, A, Rt, Tt, At, tx, qx, dtx, sctx, tocras, por, 2*N, stchsr, Tpct, b0, c0, mk, svfdv, 0);
		k=9; s=pt[k]; k=7; mk=pt[k]; k=5; b0=pt[k]; k=4; tt=pt[k]; k=10; r=pt[k]; delete[]pt;
		if ((mk<mkm) && (mk>0.0)) { mkm=mk; bo=bc; b0=bc; lx[sctxv] = r; oox[sctxv] = s; tx[sctx] = tt; qx[sctxv] = qob; } //cout << "\tbc = " << bc << "\tmk = " << mk << "\tmkm = " << mkm << endl;
		m++; if (bc<hf) hb=hb0/kdb; else hb=hb0; bb=bb-hb;
	} cout << "\tbop = " << bo << "\tmop = " << mkm << "\tcop = " << co << endl;
	if (fabs(tnr-tt)>1e-1)
		tnr=tt; else {
			cn=co-hc0; if (cn<=0.0) cn = 1e-2;
			ck=co+2e0*hc0;
			bn=bo-hb0; if (bn<=0.0) bn=1e-2;
			bk=bo+2e0*hb0;
			tnr=tt;
		}
	delete[]ao; return rez;
}
double KorrZnachVozdProsSham(double hps, double ksf, double por)
{
	int j = 0, k = 1000; double pa = 1e-3, pb = 1e0, pc, ra = fabs(pa - pb); 
	dpcts=opredKTPTKTochSha(dpctsm, etesha, etesha[vtsh], cem);
	double fa, fb, fc, tca, tcb, tcc, ka = hps*pa / por, kb = pb*hps / por, kc; //cout << "hps = " << hps << "\tksuf = " << ksf << "\tpor = " << por << endl;
	while ((ra > tocrassha) && (j < k)) { //подтягиваем пористость к значению, которое задали изначально, во время подстройки ЭКТП
		pc = (pa + pb) / 2e0;
		kc = hps*pc / por;
		tcc = kc*(1e0 - pc) / pc;
		fc = (1e0 - dpcts)*kc / (kc + tcc) - por;
		ka = hps*pa / por;
		tca = ka*(1e0 - pa) / pa;
		fa = (1e0 - dpcts)*ka / (ka + tca) - por;
		kb = hps*pb / por;
		tcb = kb*(1e0 - pb) / pb;
		fb = (1e0 - dpcts)*kb / (kb + tcb) - por;
		if ((fc*fb > 0.0) && (fa*fc<0.0)) pb = pc;
		if ((fc*fa>0.0) && (fb*fc < 0.0)) pa = pc;
		j++; ra = fabs(pa - pb);
	} //скорректированное значение толщины воздушной прослойки (размер поры), когда ввели перемычку
	return kc;
}
double RaschAlphaTvKar(int vfv, int vsv, int isrp, int vpkf, double por)
{
	int l=3, k=0, fr=0, nfr=4, d=10, q=0, cfk=l*nfr, cei=0;
	long j=0, jk=100*100*100, pj=0, h=0;
	double x=0.0, y=0.0, p=0.0, xt=0.0, yp=0.0, hf=1e0, te0=273.15, tc=22.0+te0, altc=1e-6;
	double *srra=NULL, *porm=new double[nfr], srp=0.0, e=1e-6;
	double *pn=new double[d], *alx=new double[d], *rcf=new double[cfk];
	double *rapo=NULL, **unau=NULL, ce=0.0, cet=0.0, *pp=NULL, lamtem=0.0; 
	double *legr=NULL, *prgr=NULL, dvpn=1448.0*altc/tc/4e0; //меньше 1 %
	unau=rasPorpoRazVer(por, vfv, 1, vsv, isrp, vpkf);
	k=0; rapo=unau[k]; k++; srra=unau[k]; k++; prgr=unau[k]; k++; 
	legr=unau[k]; k++; pp=unau[k]; j=0; srp=pp[j]; if (pp) delete[]pp; k++; pp=unau[k];
	ce=pp[j]; cet=ce; if (pp) delete[]pp; if (unau) delete[]unau;
	j = 0; x=e; while (x<cet) { x=x+hf; j++; } cei = j; //cout << "cem_srp = " << cei << "\tsrp = " << srp << "\t"; cout << endl; //for (j=0; j<cei; j++) if (j<10) cout << "j = " << j << "\trpr = " << rapo[j] << "\tlegr = " << legr[j] << "\t"; cout << endl;
	if (porm) {
		x = 0.0; yp = 0.0; xt = 0.0;
		for (j = 0; j < cei; j++) {
			if (yp <= 3e1) x = x + srra[j] * rapo[j];
			if (yp <= 8e1) xt = xt + srra[j] * rapo[j]; yp = yp + hf; }
		k=0; porm[k]=x*por/srp; k++; porm[k]=xt*por/srp; k++; porm[k]=porver; k++; porv[k]=por; }
	else { cout << snmv << endl; j = getchar(); exit(1); } //for (j=0; j<nfr; j++) cout << "por = " << porv[j] << "\t";
	double pr = 0.0, rmf = 0.0, prf = 0.0, po = 0.0, *alsf = new double[cfk];
	for (j = 0; j < RAND_MAX; j++) rmf = rmf + hf; for (j = 0; j < jk; j++) po = po + hf;
	if ((pn) && (alx) && (rcf) && (alsf))
	{
		for (j = 0; j < d; j++) { pn[j] = 0.0; alx[j] = 0.0; }
		for (j = 0; j < cfk; j++) { rcf[j] = 0.0; alsf[j] = 0.0; }
	}
	else { cout << snmv << endl; j = getchar(); exit(1); }
	long lt; unsigned int st; lt = time(NULL); st = (unsigned int)(lt - (lt % 2)) / 2; srand(st);
	for (fr = 0; fr < nfr; fr++) {
		for (k = 0; k < l; k++) {
			if (!fr) x = RasFracXeffVer60(k);
			if (fr == 1) x = RasFracXeffVer60_100(k);
			if (fr == 2) x = RasFracXeffVer100_150(k);
			if (fr == 3) x = RasFracXeffVer150_200(k); //размер частицы
			rcf[k + fr*l] = x*(1e-6); y = x*porv[fr]; //размер поры //cout << "x = " << x << "\ty = " << y << "\tp = " << porv[fr] << endl;
			for (j = 0; j < d; j++) pn[j] = 0.0;
			for (j = 0; j < jk; j++) {
				pj = rand(); prf = 0.0; for (h = 0; h < pj; h++) prf = prf + hf; pr = prf / rmf;
				yp = y*pr; xt = yp*(hf - porv[fr]) / porv[fr];
				p = x / (xt + yp);
				pr = 0.0; for (h = 0; h < d; h++)	{ xt = pr + hf; if ((p >= pr) && (p < xt)) pn[h] = pn[h] + hf; pr = xt; }
			}
			pr = 0.0; for (j = 0; j < d; j++) { pn[j] = pn[j] / po; pr = pr + pn[j]; } //cout << "Summa = " << pr << endl; for (j=0; j<d; j++) cout << "pn ( " << j << " ) = " << pn[j] << "\t"; cout << endl;
			for (j = 0; j < d; j++) pn[j] = pn[j] / pr;
			for (j = 2; j<d; j++) {
				p = 0.0; for (h = 0; h<j; h++) p = p + hf;
				yp = porv[fr] * x*(1e-6) / (p - hf); xt = (hf - porv[fr])*x*(1e-6) / p;
				if (yp>dvpn) {
					lamtem = F0_lamT(yp*tc); lamtem = bbfn(yp*tc*1e6); 
					if (lamtem <= 0.0) lamtem = 0.0; if (lamtem>hf) lamtem = hf;
					alx[j] = BolTochRasAlpha(0, j, yp, xt, tc, ktpvover, vtever, etev, cemv, dmkvover, snmv, Tver, Rver, 2 * N, kttkv, Aver, 0, Raver, Taver, Aaver, Rtver, Ttver, Atver)*lamtem;
				}
				else alx[j] = 0.0;
			}
			unau = RaschRTAVer(d, xt, 0.0, 0.0, 1, 0, yp, 1, tc, 0, 0, 0); pp = unau[0]; altc = pp[0]; delete[]unau; delete[]pp;
			alx[0] = 0.0; alx[1] = altc / (1e0 - porv[fr]); //for (j=0; j<d; j++) cout << "j = " << j << "\talx = " << alx[j] << "\t"; cout << endl;
			q = k + fr*l; for (j = 1; j<d; j++) {
				p = 0.0; for (h = 0; h <= j; h++) p = p + hf;
				yp = porv[fr] * x*(1e-6) / (p - hf); if (j>1) xt = oprProcSoderpoPoris(rapo, legr, yp, cei); else xt = 1e0;
				alsf[q] = alsf[q] + pn[j] * alx[j] * xt;
			}
		}
	}
	x = 0.0; yp = 0.0; for (j = 0; j < cfk; j++) { x = x + alsf[j]; yp = yp + hf; } x = x / yp; //for (j=0; j<cfk; j++) cout << "j = " << j << "\tal_sr = " << alsf[j] << "\t"; cout << endl;
	delete[]rapo; delete[]srra; delete[]rcf; delete[]pn; delete[]alx; 
	delete[]alsf; delete[]legr; delete[]porv; delete[]prgr;
	x = x / altc; printf("sr_okp = %0.10lf\n", x); //ослабление КП за счет пористой структуры вермикулита
	return x;
}