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
const int N = 13, ks = 10;
const double te0=273.15, tocrassha=1e-8;
struct derevo {
	int otre; //1 - отражение или 0 - пропускание
	int ste; //номер стенки
	int vis; //видимость: 1 - виден, 0 - нет
	int lev; //номер уровня
	struct derevo *back; //указатель назад
	struct derevo *next; //указатель вперед
};
class KoePog { public: double alp, tem; KoePog *nex; };
double *dkosct = NULL, *kttk = NULL, *stchsrsha = NULL, *ktr = NULL, *ecktpsha = NULL, *Tpct = NULL, *dkoscm = NULL; //dkoscm - отклонение СЧ (или ПС) от литературных данных
double *mko, *Th = NULL, *Tg = NULL, *kektp = NULL, *ektpsha = NULL, *ktpvosha = NULL, *vtesha = NULL, tmasha, tmisha, dpcts = 0.314148; //dpcts - доля площади перемычки, через которую переносится чистой теплопроводностью ТК
double *alphaSha = NULL, *txsha = NULL, *kxsha = NULL, cop = 0.0, bop = 0.0, tnr = 0.0, bn = 1e-2, bk = 2e0 / pow(1e0 - dpcts, 2e0);
double cn = 1e-2, ck = 4e0, *temrass=NULL, *ooxsha=NULL, sfeosha = 1.24e-2, smgosha = 0.29e-2, salosha = 41.9e-2, ssiosha = 54e-2, *dpctsm;
int vtsh = 0, sctx = 0, dkoscl = 6, vystsha = 1, vybsha = 3;
char *sns = NULL, *sss = NULL, *skpts = NULL, *skptsk = NULL, *svsk = NULL, *svs = NULL, *snms = NULL, *sfnos = NULL;
//----------------------------------------
double *podschchieleSha(int, int, int, double *, double *);
double rasotprpovsSha(int *, int *, int, int, int, int, double *, double *);
double otrprovsSha(int, int, double *, double *);
double *koefPribSha(double *, double *, int, double *);
double *opredKoefOtr(double *, double *, int, double, int, int, double *, double *, int);
void zadrktShaNac();
void NapMasVozdSha(double *, double *, int);
double *opredKTPTverKarkSha(double *, double *, double, double, double, double, double, double, int, double *, double *, double *, int);
void opredtemphc(double *, double *, int, int, double);
double opredKTPTKTochSha(double *, double *, double, int);
double *usrednen(double *, double *, int, int);
void vyvodfile(double *, int, int, double, char *);
void napstrsha();
double *reshnewtrafs(double *, double *, double *, double *, int, double *, double *, int, double *, double *, int, double *, double *, int, double, double *, double *, int, double, double, double, double *, int, double, double, double);
double *reshnewrafsokon(double *, double **, double **, double *, int, int, double, double *, double *, int, double, double *, int, double *, int, double, double, double, double, double, double *, double *, double *, int);
void PrintMatr(double **, int);
void osvopam();
void iniarrele(int, double, double, double);
double **opredTempLPStenSha(double *, double *, double *, double *, double, int, double *, double *, double *, double *, int, int, double, double, double, int, char *);
double **RaschVneshIzluchSha(double *, double *, double *, double *, double *, int, double, double, char *);
double **RasLuchPloTepPot(int, double *, double *, double *, double *, double *, double *, int, double *, double *, double *, char *);
double **RaschSobLuchPlotTepPot(int, double *, double *, double *, double *, double *, double *, double, double, double *, int, char *);
double opredUrovPodderM03(double);
double opredUrovOtsech(double);
double *opredTempStenShaFragm(double *, int, double *, double *, double *, double *, int, int, double, double, double, double, double);
double *KoefPoglRosselNac(double *, int, int, double, double, double, double *, double *, int, double, double, double *, double *, int, int, int);
double LuchKTPChudnovsky(double *, double, int, double);
double epsisredsha(double, double *, double *, int, double *, double *, int);
double **RaschRTASha(int, double, double, double, double, double, int, int, int, int);
double *oprRasTemNach(int, int, double *, double *, int, double *, double, double, double, double, double, int);
double *oprEffDoliTepPeren(double, double, double);
double *KorrZnachVozdPros(double, double, double, int);
double *koefoslab(double, double, double, double *, int, double *);
double vychNevyaz(double **, double **, double *, double *, int);
double *opredPolTempTvKarShaFragm(double *, int, double *, double *, int, double, double, double, double, double);
double *EffectTols(double *, double *, double *, double, double, int);
double SeryeStenkiRasIzl(double, double, double, double *, double *, double *, double *, double *, int *, int);
double oprProcSoderpoPoris(double *, double *, double, int);
double BolTochRasAlpha(int, int, double, double, double, double *, double *, double *, int, int, char *, double *, double *, int, double *, double *, int, double *, double *, double *, double *, double *, double *);
double RasFracXeffSha60(int);
double RaschAlphaTvKarSha();
double F0_lamT(double);
double RasIzlSerStenNac(double *, double *, double *, double *, double *, double *, double *, double, double, double, double, int, double, double, double, double, int, double *, double *, int, int, double *, double *, double *, char *, int, double *, double *, int);
double KorrZnachVozdProsSham(double, double, double);
double *rasMetNewtRafs(double *, double **, double **, double *, double, double, int);
double *FuncRaschIzl(double, double, double, double, int, double *, double, double *, double, double, double, double, int, double *, double *, int, int, char *, double, int, double *, double, double *, double *, int, double *, double *, int, double *, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, int, double, double, int, double *, double *, double, double, double, char *, int, double);
double *PoiskZavVelTemVer(int, double *);
void zapisvfile(double *, int, char *);
//----------------------------------------
double *podschchieleSha(int no, int kste, int ocs, double *Ref, double *Tra)
{
	int kol = 0, *ll = new int[N], *kl = new int[N], p, pe = 0, ee, qs = 0, kk = 0, k = 0, j, jj = 1, w, b = 0, bb = 0, gg = 1, *st = new int[N], *ot = new int[N];
	double pp, *na = new double[4], nat = 0.0; struct derevo *prev, *prevv = NULL, *roo, *poir, *poil, *poi, **pruk = new struct derevo*[1], **ceuk = new struct derevo*[2], **unpe = new struct derevo*[N];
	prev = new struct derevo;
	if ((!na) || (!kl) || (!ll) || (!prev) || (!pruk) || (!ceuk) || (!st) || (!ot) || (!unpe)) { cout << snms << endl; p = getchar(); exit(1); }
	for (j = 0; j < 4; j++) na[j] = 0.0;
	kol = 0; kk = kol; prev->lev = gg; gg++; prev->ste = no; prev->vis = 1; prev->next = NULL; roo = prev; kol++; prev->otre = -1; prev->back = NULL; pruk[0] = prev; unpe[b] = prev; b++; kl[k] = kol - kk; k++;
	//-----------
	bb = 2; for (j = 0; j < bb; j++) ceuk[j] = NULL; kk = kol; pe = 0;
	poir = new struct derevo; poil = new struct derevo; if ((!poir) || (!poil)) { cout << snms << endl; j = getchar(); exit(1); }
	poil->back = prev; poil->otre = 0; poir->lev = gg; poil->lev = gg; poir->back = prev; poir->otre = 1; poir->ste = prev->ste + 1; poil->ste = prev->ste - 1;
	poir->vis = 1; if ((poir->ste<1) || (poir->ste>ocs)) poir->vis = 0; if (poir->vis == 1) { ceuk[pe] = poir; kol++; pe++; }
	else delete poir;
	poil->vis = 1; if ((poil->ste<1) || (poil->ste>ocs)) poil->vis = 0; if (poil->vis == 1) { ceuk[pe] = poil; kol++; pe++; }
	else delete poil;
	j = 0; poi = ceuk[j]; prev = unpe[b - 1]; gg++;
	while ((poi) && (j < pe))
	{
		if (poi->vis == 1)
		{
			prev->next = poi;
			prev = poi;
		}
		j++; poi = ceuk[j];
	}
	unpe[b] = prev; b++;
	kl[k] = kol - kk; k++;
	delete[]pruk; pruk = new struct derevo*[pe]; if (!pruk) { cout << snms << endl; j = getchar(); exit(1); }
	qs = 0; for (j = 0; j<pe; j++) { pruk[j] = NULL; poi = ceuk[j]; if ((poi) && (poi->vis == 1)) { pruk[qs] = poi; qs++; } } jj = qs; delete[]ceuk;
	//------------
	//0 - пропускание (left), 1 - отражение (right)
	pp = 3.0; for (p = 3; p <= N; p++) {
		w = floor(pow(2.0, pp - 1.0)); bb = 2 * pe; pp = pp + 1.0;
		if (bb>0) { ceuk = new struct derevo*[bb]; if (!ceuk) { cout << snms << endl; j = getchar(); exit(1); } else for (j = 0; j < bb; j++) ceuk[j] = NULL; }
		pe = 0; qs = 0; ee = 0; kk = kol;
		while (ee<w) {
			if ((qs<jj) && (pruk) && (ceuk)) {
				prev = pruk[qs];
				if (prev) {
					prevv = prev->back;
					if (prevv)
					{
						poir = new struct derevo; poil = new struct derevo; if ((!poir) || (!poil)) { cout << snms << endl; j = getchar(); exit(1); }
						poir->back = prev; poir->otre = 1; poil->back = prev; poil->otre = 0; poir->next = NULL; poil->next = NULL; poir->lev = gg; poil->lev = gg;
						if (prev->otre == 1) { if (prevv->ste>prev->ste) { poil->ste = prev->ste - 1; poir->ste = prev->ste + 1; } else { poil->ste = prev->ste + 1; poir->ste = prev->ste - 1; } }
						if (prev->otre == 0) if (prevv->ste>prev->ste) { poil->ste = prev->ste - 1; poir->ste = prev->ste + 1; }
						else { poil->ste = prev->ste + 1; poir->ste = prev->ste - 1; }
						if (prev->otre == -1) if (prevv->ste > prev->ste) { poil->ste = prev->ste - 1; poir->ste = prev->ste + 1; }
						else { poil->ste = prev->ste + 1; poir->ste = prev->ste - 1; }
						poir->vis = 1; if ((poir->ste<1) || (poir->ste>ocs) || (prev->vis == 0)) poir->vis = 0; if (poir->vis == 1) { ceuk[pe] = poir; kol++; pe++; }
						else delete poir;
						poil->vis = 1; if ((poil->ste<1) || (poil->ste>ocs) || (prev->vis == 0)) poil->vis = 0; if (poil->vis == 1) { ceuk[pe] = poil; kol++; pe++; }
						else delete poil;
						ee = ee + 2; qs++;
					}
				}
			} ee++;
		}
		pruk = new struct derevo*[pe]; if (!pruk) { cout << snms << endl; j = getchar(); exit(1); }
		qs = 0; for (j = 0; j < pe; j++) { poi = ceuk[j]; if ((poi) && (poi->vis == 1)) { pruk[qs] = poi; qs++; } } jj = qs;
		kl[k] = kol - kk; k++;
		j = 0; poi = pruk[j]; prev = unpe[b - 1];
		while ((poi) && (j < pe))
		{
			if (poi->vis == 1)
			{
				prev->next = poi; prev = poi;
			}
			j++; poi = pruk[j];
		} unpe[b] = prev; b++; gg++; delete[]ceuk;
	}
	j = 0; for (p = 0; p < N; p++) { j = j + kl[p]; ll[p] = j; }
	//------------------------------
	for (j = 0; j < 4; j++) na[j] = 0.0; bb = 1; poi = roo;
	for (k = 1; k <= N; k++) {
		kk = bb - 1; b = 0;
		while ((b < kl[kk]) && (poi)) {
			if (poi->lev == bb) {
				for (w = 0; w < N; w++) { ot[w] = 0; st[w] = 0; }
				prev = poi;
				for (w = 0; w<bb; w++)
				{
					if (poi) {
						st[kk - w] = poi->ste;
						ot[kk - w] = poi->otre;
						poi = poi->back;
					}
				}
				poi = prev;
				nat = 0.0; nat = rasotprpovsSha(st, ot, kk, kste, no, bb, Ref, Tra);
				if (nat>0) {
					w = 0; //cout << "nat = " << nat << endl; 
					if (st[kk] == kste)
					{
						if (st[kk - 1]<kste) w = 0; else w = 1; if (kste == 1) w = 1; if (kste == ocs) w = 0; na[w] = na[w] + nat;
						w = 2; if (st[1]<no) w = 2; else w = 3; if (no == 1) w = 3; if (no == ocs) w = 2; na[w] = na[w] + nat; //3 - излучение идет вправо от излучающей стенки, 2 - влево
					}
				} //0 - lr (излучение падает на рассматриваемую стенку слева, т.е. идет слева - направо), 1 - rl (излучение от других стенок идет справа, т.е. распространяется справа - налево)
				poi = poi->next; b++;
			}
		} bb++;
	}
	poi = roo; while (poi)  { prev = poi->next; delete poi; poi = prev; }
	delete[]ot; delete[]st; delete[]pruk; delete[]kl; delete[]ll; delete[]unpe; //getchar();
	return na;
}
double rasotprpovsSha(int *stt, int *ott, int gg, int kst, int izlst, int kost, double *R, double *T)
{
	int e = 1, ee = 0, f = 1, w; double po = 1.0, pot;
	if ((stt[gg] == kst) && (kost > 1)) {
		while (e <= gg) {
			w = stt[ee]; if (izlst == stt[0]) { if (f > 0) { e++; ee++; f = 0; continue; } }
			pot = otrprovsSha(w - 1, ott[e], R, T);
			po = po*pot; ee++; e++;
		}
	}
	else po = 0.0; return po;
}
double otrprovsSha(int ste, int otr, double *R, double *T)
{
	double ko = 1e0;
	if (otr == 1) { ko = ko*R[ste]; }
	else if (!otr) { ko = ko*T[ste]; }
	return ko;
}
double *koefPribSha(double *ktp, double *te, int le, double *ko)
{
	int k, kem = 3; double **A = new double*[kem], *AA, *b = new double[kem], yx2 = 0.0, yx = 0.0, p, hf = 1e0;
	double x4 = 0.0, x3 = 0.0, x2 = 0.0, x = 0.0, y = 0.0, de = 0.0, de1 = 0.0, de2 = 0.0, de3 = 0.0;
	if (A) { for (k = 0; k < kem; k++) { AA = new double[kem]; if (AA) A[k] = AA; else { cout << snms << endl; k = getchar(); exit(1); } } }
	else { cout << snms << endl; k = getchar(); exit(1); }
	if ((!ko) || (!b)) { cout << snms << endl; k = getchar(); exit(1); }
	for (k = 0; k < le; k++) {
		yx2 = yx2 + ktp[k] * pow(te[k], 2e0); yx = yx + ktp[k] * te[k]; y = y + ktp[k];
		x4 = x4 + pow(te[k], 4e0); x3 = x3 + pow(te[k], 3e0); x2 = x2 + pow(te[k], 2e0); x = x + te[k];
	}
	b[0] = yx2; b[1] = yx; b[2] = y; p = 0.0; for (k = 0; k < le; k++) p = p + hf;
	A[0][0] = x4; A[0][1] = x3; A[0][2] = x2; A[1][0] = x3; A[1][1] = x2; A[1][2] = x; A[2][0] = x2; A[2][1] = x; A[2][2] = p;
	de = A[0][0] * (A[2][2] * A[1][1] - A[2][1] * A[1][2]) - A[0][1] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]) + A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
	de1 = b[0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]) - A[0][1] * (b[1] * A[2][2] - b[2] * A[1][2]) + A[0][2] * (b[1] * A[2][1] - b[2] * A[1][1]);
	de2 = A[0][0] * (b[1] * A[2][2] - b[2] * A[1][2]) - b[0] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]) + A[0][2] * (A[1][0] * b[2] - A[2][0] * b[1]);
	de3 = A[0][0] * (A[1][1] * b[2] - A[2][1] * b[1]) - A[0][1] * (A[1][0] * b[2] - A[2][0] * b[1]) + b[0] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
	ko[2] = de1 / de; ko[1] = de2 / de; ko[0] = de3 / de;
	delete[]b; for (k = 0; k < kem; k++) { AA = A[k]; delete[]AA; } return ko;
}
void NapMasVozdSha(double *ktpvoz, double *tvoz, int dmvo)
{
	int k = 0, dlma = dmvo;
	ktpvoz[k] = 2.44; k++; ktpvoz[k] = 2.51; k++; ktpvoz[k] = 2.59; k++; ktpvoz[k] = 2.67; k++; ktpvoz[k] = 2.76; k++; ktpvoz[k] = 2.83; k++;
	ktpvoz[k] = 2.90; k++; ktpvoz[k] = 2.96; k++; ktpvoz[k] = 3.05; k++; ktpvoz[k] = 3.13; k++; ktpvoz[k] = 3.21; k++; ktpvoz[k] = 3.34; k++;
	ktpvoz[k] = 3.49; k++; ktpvoz[k] = 3.64; k++; ktpvoz[k] = 3.78; k++; ktpvoz[k] = 3.93; k++; ktpvoz[k] = 4.27; k++; ktpvoz[k] = 4.60; k++;
	ktpvoz[k] = 4.91; k++; ktpvoz[k] = 5.21; k++; ktpvoz[k] = 5.74; k++; ktpvoz[k] = 6.22; k++; ktpvoz[k] = 6.71; k++; ktpvoz[k] = 7.18; k++;
	ktpvoz[k] = 7.63; k++; ktpvoz[k] = 8.07; k++; ktpvoz[k] = 8.50; k++; ktpvoz[k] = 9.15;
	for (k = 0; k < dlma; k++) ktpvoz[k] = ktpvoz[k] / 1e2; //КТП воздуха
	k = 0;
	tvoz[k] = 0.0; k++; tvoz[k] = 1e1; k++; tvoz[k] = 2e1; k++; tvoz[k] = 3e1; k++; tvoz[k] = 4e1; k++; tvoz[k] = 5e1; k++;
	tvoz[k] = 6e1; k++; tvoz[k] = 7e1; k++; tvoz[k] = 8e1; k++; tvoz[k] = 9e1; k++; tvoz[k] = 1e2; k++; tvoz[k] = 12e1; k++;
	tvoz[k] = 14e1; k++; tvoz[k] = 16e1; k++; tvoz[k] = 18e1; k++; tvoz[k] = 2e2; k++; tvoz[k] = 25e1; k++; tvoz[k] = 3e2; k++;
	tvoz[k] = 35e1; k++; tvoz[k] = 4e2; k++; tvoz[k] = 5e2; k++; tvoz[k] = 6e2; k++; tvoz[k] = 7e2; k++; tvoz[k] = 8e2; k++;
	tvoz[k] = 9e2; k++; tvoz[k] = 1e3; k++; tvoz[k] = 11e2; k++; tvoz[k] = 12e2;
	for (k = 0; k < dlma; k++) tvoz[k] = tvoz[k] + te0; //температуры, при которых определяется КТП воздуха
}
double *FuncRaschIzl(double ta, double tb, double d, double tc, int m, double *prot, double de, double *ao, double rez, double rap, double b, double c, int vybo, double *qobve, double *eteve, int cemve, int vyve, char *snm, double dpct, int vyte, double *ektpve, double y0ve, double *ktrve, double *kttkve, int kost, double *ktpvove, double *vteve, int dmkvove, double *ecktpve, double htch, double hvoz, double *Refe, double *Trae, double *Abse, double *Reft, double *Trat, double *Abst, double *txve, double *qxve, double *dtxve, int sctve, double tocrasve, double poristost, int dmao, double *stchsrve, double *Tpctve, double b0, double c0, double mk, char *navyfa, int fla, double hksha)
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
	tevy = reshnewtrafs(prot, Tna, hlr1, hrl1, zf, Trae, Abse, kost, kttkve, eteve, cemve, ktpvove, vteve, dmkvove, qs, txve, qxve, sctve, htch, hvoz, tocrasve, ektpve, cels + 2 * kost, tb, ta, y0ve);
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
double *usrednen(double *ao, double *usr, int sao, int kost)
{
	double s = 0.0, k = 0.0, minz = tocrassha, hf = 1e0; int j;
	for (j = 0; j<kost; j++) if (fabs(usr[j])>minz) { s = s + usr[j]; k = k + hf; }
	s = s / k; ao[sao] = s; //средняя величина
	return ao;
}
double *reshnewtrafs(double *w, double *te, double *hlr, double *hrl, int zf, double *Tm, double *Am, int kost, double *ktptc, double *ete, int chel, double *ktpvo, double *tvo, int dmkvo, double qc, double *tx, double *qx, int scx, double htk, double hvo, double tora, double *laefm, int koelvyma, double tmax, double tmin, double tolob)
{
	int k, j = 0, q, p = 0, dm = 2 * kost, m = 0; double **m1 = new double*[dm], **m4 = new double*[dm], *te2 = new double[dm], *vykh;
	double sig = 5.67e-8, lav, latk, laef, t, wo, wr, wt, w0, w1, w2, w3, *t1, *t4, s, *ssc = new double[dm];
	if ((!m1) || (!m4) || (!ssc) || (!te2)) { cout << snms << endl; k = getchar(); exit(1); }
	for (k = 0; k < dm; k++) te2[k] = te[j];
	p = 0; for (k = 0; k<dm; k++) {
		q = 0;  //k - номера стенок для расчет температур, штрих - слева, два штриха - справа
		t1 = new double[dm]; t4 = new double[dm]; if ((!t1) || (!t4)) { cout << snms << endl; q = getchar(); exit(1); }
		if (k>0) p = p + ((k + 1) % 2); //cout << "p = " << p << endl;
		for (j = 0; j<dm; j++) { t1[j] = 0.0; t4[j] = 0.0; }
		s = 0.0; j = 0; m = 0;
		while (j<dm) {
			if (j>0) m = m + ((j + 1) % 2); //cout << "m = " << m << "\t"; //2 - излучение пошло влево от излучающей стенки, знак минус, 3 - вправо, знак плюс
			w0 = w[q]; w1 = w[q + 1]; w2 = w[q + 2]; w3 = w[q + 3]; wo = w2 + w3; //cout << "w0 = " << w0 << "\tw1 = " << w1 << "\tw2 = " << w2 << "\tw3 = " << w3 << endl;
			if (!(k % 2)) wr = w0 - w1*Tm[p]; else wr = w0*Tm[p] - w1;
			if (fabs(wo)>0.0) wt = Am[m] * wr*sig / wo; else wt = 0.0;
			if (!(j % 2)) { wt = wt*w2; }
			else { wt = wt*w3; } t4[j] = wt;
			if (fabs(wo) > 0.0) {
				s = s + wr*(w2 + w3)*(hlr[m] - hrl[m]) / wo; //cout << "wr = " << wr << "\tw2 = " << w2 << "\tw3 = " << w3 << "\thlr = " << hlr[m] << "\t"; 
			} if (j == k) {
				t = te[m]; latk = opredKTPTKTochSha(ktptc, ete, t, chel);
				lav = opredKTPTKTochSha(ktpvo, tvo, t, dmkvo); //1 - излучение приходит слева на конечную стенку, знак плюс, 0 - справа, знак минус
				if (!j) { t1[j] = latk / htk; t1[j + 1] = -t1[j]; }
				else if (j == (dm - 1)) { t1[dm - 2] = latk / htk; t1[dm - 1] = -t1[dm - 2]; }
				else if (!(j % 2)) { t1[j - 1] = lav / hvo; t1[j] = latk / htk - lav / hvo; t1[j + 1] = -latk / htk; }
				else { t1[j - 1] = latk / htk; t1[j] = lav / hvo - latk / htk; t1[j + 1] = -lav / hvo; }
			}
			if (j > 0) q = q + 4 * ((j - 1) % 2); j++;
		}
		ssc[k] = -s; m1[k] = t1; m4[k] = t4;
	} //cout << "m1 = " << endl; PrintMatr(m1,dm); //cout << "m4 = " << endl; PrintMatr(m4,dm); for (j=0; j<dm; j++) cout << "ssc (" << j << ") = " << ssc[j] << "\t" << "te (" << j << ") = " << te2[j] << "\t"; cout << endl; 
	vykh = reshnewrafsokon(te2, m1, m4, ssc, zf, dm, qc, tx, qx, scx, tora, laefm, chel, ete, koelvyma, tmax, tmin, tolob, htk, hvo, ktptc, ktpvo, tvo, dmkvo);
	for (k = 0; k < dm; k++) { t1 = m1[k]; delete[]t1; } delete[]m1;
	for (k = 0; k < dm; k++) { t4 = m4[k]; delete[]t4; } delete[]m4;
	delete[]ssc; delete[]te2; return vykh;
}
double *reshnewrafsokon(double *te, double **ms1, double **ms4, double *ssc, int zf, int ksu, double qc, double *tx, double *qx, int scx, double tora, double *laefm, int chis, double *ete, int koelvyma, double tmax, double tmin, double tolob, double htch, double hvozd, double *ktptk, double *ktpvo, double *vte, int dmv)
{ //cout << "Newton-Rafsons is on" << endl;
	int qk = 1, k, q, j=10, p, f=4; for (p=0; p<f; p++) qk=qk*j;
	double *mvyz = new double[koelvyma], ta = tmin, tb = tmax, tc, ep=1e-6;
	double *tet = new double[ksu], ksuf = 0.0, ht=1e0, tempplav=18e2;
	double fa, fb, fc, tmaa, tmab, e = 1e1, ra = 1e3, t, r, eps=1e-10;
	if ((!mvyz) || (!tet)) { cout << snms << endl; k = getchar(); exit(1); }
	for (k = 0; k < ksu; k++) { tet[k] = te[k]; ksuf = ksuf + ht; }
	tet = rasMetNewtRafs(tet, ms1, ms4, ssc, ksuf, tora, ksu);
	for (k=0; k<ksu; k++) { if ((tet[k]<ta) || (tet[k]>tb)) tet[k]=0.0; } //zapisvfile(tet, ksu, "C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\2022\\Temp_Stenok.txt"); zapisvfile(ssc, ksu, "C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\2022\\Svob_ch.txt"); for (k=0; k<ksu; k++) zapisvfile(ms1[k], ksu, "C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\2022\\Matr_1.txt"); for (k=0; k<ksu; k++) zapisvfile(ms4[k], ksu, "C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\2022\\Matr_4.txt");
	r = 0.0; fa=r; fb=r; for (j = 0; j < ksu; j++) { t = tet[j]; if (r < t) r = t; if (t>ep) { fa=fa+t; fb=fb+ht; } } j=5; mvyz[j] = r; 
	if (fabs(fb)>eps) fa=fa/fb; else fa=0.0;
	cout << "\tSred Temp NR = " << fa; //FILE *fo = fopen("C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\2022\\Temp_NR.txt", "at"); fprintf(fo, "%0.15lf\n", fa); fclose(fo);
	for (k = 0; k < ksu; k++) tet[k] = te[k];
	q = 0; while ((q<qk) && (ra>e)) {
		tc = (ta + tb) / 2e0;
		tet = opredTempStenShaFragm(tet, ksu, ktpvo, vte, ete, ktptk, chis, dmv, htch, hvozd, qc, ta, ht);
		fa = vychNevyaz(ms1, ms4, tet, ssc, ksu); tmaa = 0.0; for (k = 0; k < ksu; k++) if (tmaa < tet[k]) tmaa = tet[k];
		tet = opredTempStenShaFragm(tet, ksu, ktpvo, vte, ete, ktptk, chis, dmv, htch, hvozd, qc, tb, -ht);
		fb = vychNevyaz(ms1, ms4, tet, ssc, ksu); tmab = 0.0; for (k = 0; k < ksu; k++) if (tmab<tet[k]) tmab = tet[k];
		tet = opredTempStenShaFragm(tet, ksu, ktpvo, vte, ete, ktptk, chis, dmv, htch, hvozd, qc, tc, ht);
		fc = vychNevyaz(ms1, ms4, tet, ssc, ksu);
		if ((fc*fb>0.0) && (fa*fc<0.0)) tb = tc; if ((fc*fa>0.0) && (fb*fc < 0.0)) ta = tc;
		ra = fabs(tmab - tmaa);
		q++;
	}
	tet = opredTempStenShaFragm(tet, ksu, ktpvo, vte, ete, ktptk, chis, dmv, htch, hvozd, qc, tc, 1.0);
	t = 0.0; for (k = 0; k < ksu; k++) for (j = 0; j<ksu; j++) { tc = tet[k] - tet[j]; if (tc>t) t = fabs(tc); } //максимальная разность температур
	j=0; mvyz[j] = t;
	j=4; mvyz[j] = t / (htch*ksuf + hvozd*(ksuf - 1e0));
	j=1; mvyz[j] = qc; //плотность теплового потока
	t = 0.0; for (k = 0; k < ksu; k++) t = t + tet[k]; t = t / ksuf; //средняя температура конечная //cout << "T_sred_tet = " << t << "\tT_sred = " << tc << "\ttmin = " << tmaa << "\ttmax = " << tmab << "\tra_tem = " << ra << endl;
	j=2; mvyz[j] = t;
	t = 0.0; for (k = 0; k < ksu; k++) t = t + te[k]; t = t / ksuf; //средняя температура начальная
	j=3; mvyz[j] = t;
	f = koelvyma - ksu; for (j = f; j < koelvyma; j++) mvyz[j] = tet[j - f]; //распределение температур
	delete[]tet;
	return mvyz;
}
void vyvodfile(double *ao, int sao, int zf, double hvv, char *nafa)
{
	int k; FILE *fo = fopen(nafa, "at");
	if (!fo) { cout << sfnos << endl; k = getchar(); exit(1); }
	if (zf < 1) {
		fprintf(fo, "\nhv = %0.10lf\n", hvv); 
	} for (k = 0; k < sao; k++) fprintf(fo, "%0.15lf\n", ao[k]); fprintf(fo, "\n"); fclose(fo);
}
double opredKTPTKTochSha(double *ktptks, double *te, double temp, int ce)
{
	int n = ce, f = 1, p = 0, k; 
	double e=1e-4, ktp = 0.0, ko = 0.0, x1=0.0, x2=0.0, y1=0.0, y2=0.0, b=0.0, dt=0.0; 
	if ((temp>=te[0]) && (temp<=te[n-1])) {
	for (k = 0; k<n; k++)
		if ((te[k] >= temp) && (f>0)) { p = k; f = 0; break; } 
	if (!p) { if (!f) p = 1; else { p = n - 1; f = 0; }	} }
	else if (temp<te[0]) { p=1; f=0; }
	else if (temp>te[n-1]) { p=n-1; f=0; }
	if ((!f) && (p>0)) {
		x2=te[p];
		x1=te[p - 1];
		dt = x2 - x1;
		if (fabs(dt) > e) {
			y2=ktptks[p];
			y1=ktptks[p - 1];
			b=y1;
			if ((n-1)==p) b=y2;
			ko = (y2 - y1) / dt;
			ktp = b + ko*(temp - x1);
		}
	}
	return ktp;
}
double **opredTempLPStenSha(double *Ts, double *Tss, double *Tsr, double *Tna, double Tnach, int n, double *ktpvo, double *vte, double *ete, double *ktptk, int ce, int dmv, double htk, double hnf, double ptpob, int kost, char *snm)
{
	int q = 4, m, k; double lam, tt, lamvo, lamt, ht, **mu = new double*[q];;
	Ts = new double[kost]; Tss = new double[kost]; Tsr = new double[kost]; Tna = new double[n];
	if ((!Tna) || (!Ts) || (!Tss) || (!mu)) { cout << snm << endl; k = getchar(); exit(1); }
	for (k = 0; k < kost; k++) { Ts[k] = 0.0; Tss[k] = 0.0; Tsr[k] = 0.0; } for (k = 0; k < n; k++) Tna[k] = 0.0;
	Tna[0] = Tnach; Ts[0] = Tnach; q = 1; m = 0;
	for (k = 1; k < n; k++)  {
		tt = Tna[k - 1]; lamt = opredKTPTKTochSha(ktptk, ete, tt, ce); if (lamt < 0.0) lamt = 0.0; //ищем КТП твердого скелета при данной температуре
		lamvo = opredKTPTKTochSha(ktpvo, vte, tt, dmv); if (lamvo<0.0) lamvo = 0.0; //ищем КТП воздуха при данной температуре
		if (k % 2) { ht = htk; lam = lamt; }
		else { ht = hnf; lam = lamvo; }
		if (fabs(lam)>0.0) Tna[k] = tt - ptpob*ht / lam; else Tna[k] = Tna[k - 1]; //рассчитываем температуру на границах стенок 
		if ((k % 2)) { Tss[m] = Tna[k]; m++; }
		else { Ts[q] = Tna[k]; q++; }
	}
	for (k = 0; k < kost; k++) Tsr[k] = (Ts[k] + Tss[k]) / 2e0; //рассчитываем температуру внутри стенок 
	k = 0; mu[k] = Ts; k++; mu[k] = Tss; k++; mu[k] = Tsr; k++; mu[k] = Tna; return mu;
}
double **RaschVneshIzluchSha(double *Tra, double *Ref, double *prot, double *hlr1, double *hrl1, int kost, double hlr0, double hrl0, char *snm)
{
	int k = 3, j, m, q = 2; double *khrl = new double[kost], *khlr = new double[kost], *ri = new double[q], **mu = new double*[k];
	hrl1 = new double[kost]; hlr1 = new double[kost];
	if ((!khrl) || (!khlr) || (!ri) || (!hrl1) || (!hlr1) || (!mu)) { cout << snm << endl; k = getchar(); exit(1); }
	for (k = 0; k < kost; k++) { khrl[k] = 0.0; khlr[k] = 0.0; hrl1[k] = 0.0; hlr1[k] = 0.0; } for (k = 0; k < q; k++) ri[k] = 0.0;
	double w0, w1, w2, w3, wk, slr, srl, tta, ttb, tt;
	khlr[0] = 1e0; j = kost - 1; khrl[j] = 1e0;
	for (k = 1; k < kost; k++) {
		khlr[k] = khlr[k - 1] * Tra[k - 1];
		khrl[j - 1] = khrl[j] * Tra[j];
		j--;
	} //for (k=0; k<kost; k++) cout << "khrl " << k << " = " << khrl[k] << "\tkhlr " << k << " = " << khlr[k] << endl;
	q = 0; slr = 0.0; srl = 0.0; for (j = 0; j < kost; j++) //узнаем, какая доля внешнего излучения дошла до j-ой стенки - учитываем вклад от всех стенок
	{
		hlr1[j] = 0.0; hrl1[j] = 0.0; //внешнее излучение падает с обеих сторон
		tta = 0.0; ttb = 0.0; m = q;
		for (k = 0; k<kost; k++) {
			w0 = prot[q]; w1 = prot[q + 1]; w2 = prot[q + 2]; w3 = prot[q + 3]; wk = w2 + w3; //для внешнего излучения, распространяющегося слева направо
			tt = Ref[k] * w2 + Tra[k] * w3; //доля энергии, отраженная от k-ой стенки, и прошедшая через нее, пришла справа к j-ой стенке
			if (fabs(wk)>0.0) {
				tta = tta + tt*khlr[k] * w0 / wk; //доля энергии, отраженная от k-ой стенки, и прошедшая через нее, пришла слева к j-ой стенке
				ttb = ttb + tt*khlr[k] * w1 / wk;
			} //для учета при подсчете лучистой энергии, распространяющейся справа налево
			q = q + 4;
		} q = m + 4 * kost - 1; m = q; //считаем, что стенки, до которых дошло излучение, как бы "сами" излучают, из-за того, что на них попало внешнее излучение
		for (k = (kost - 1); k >= 0; k--) {
			w0 = prot[q - 3]; w1 = prot[q - 2]; w2 = prot[q - 1]; w3 = prot[q]; wk = w2 + w3; q = q - 4; //для внешнего излучения, распространяющегося справа налево
			tt = Ref[k] * w3 + Tra[k] * w2; //деление на два потока
			if (fabs(wk) > 0.0) {
				tta = tta + tt*khrl[k] * w0 / wk; //для учета лучистой энергии, распространяющейся слева направо, падает слева от j-ой стенки
				ttb = ttb + tt*khrl[k] * w1 / wk;
			} //то, что попадает справа налево на стенку, справа от j-ой стенки
		} q = m + 1;
		hlr1[j] = tta*hlr0; //узнаем величину в абсолютных единицах этой лучистой энергии - слева от пластинки, падает слева направо
		hrl1[j] = ttb*hrl0; //справа от пластинки, падает справа налево
	} slr = slr + Tra[kost - 1] * hlr1[kost - 1] + hrl0*Ref[kost - 1]; ri[0] = slr; //то, что вышло из последней стенки вправо
	srl = srl + Tra[0] * hrl1[0] + Ref[0] * hlr1[0]; ri[1] = srl; //то, что вышло из первой стенки влево
	hlr1[0] = hlr0; hrl1[kost - 1] = hrl0; //for (k=0; k<ks; k++) cout << "hrl " << k << " = " << hrl1[k] << "\thlr " << k << " = " << hlr1[k] << endl;
	delete[]khlr; delete[]khrl; //cout << "Ras Vn Iz Sha k" << endl;
	k = 0; mu[k] = hlr1; k++; mu[k] = hrl1; k++; mu[k] = ri; return mu;
}
double **RaschSobLuchPlotTepPot(int kost, double *prot, double *Ts, double *Tss, double *Tra, double *Ref, double *Ab, double slr, double srl, double *ao, int sao, char *snm)
{
	int j, q = 3, k, p = sao;
	double tta, ttb, ttc, tt, w0, w1, w2, w3, wk, sig = 5.67e-8, **mu = new double*[q];
	double *sislr = new double[kost], *sisrl = new double[kost];
	if ((!sisrl) || (!sislr) || (!mu)) { cout << snm << endl; k = getchar(); exit(1); }
	q = 0; for (j = 0; j < kost; j++) { //узнаем, какая энергия собственного излучения дошла до j-ой стенки
		sislr[j] = 0.0; sisrl[j] = 0.0; tta = 0.0; ttb = 0.0;
		for (k = 0; k<kost; k++) {
			w0 = prot[q]; w1 = prot[q + 1]; w2 = prot[q + 2]; w3 = prot[q + 3]; wk = w2 + w3;
			tt = pow(Ts[k], 4e0)*w2 + pow(Tss[k], 4e0)*w3;
			if (fabs(wk)>0.0) {
				ttc = Ab[k] * w0*tt / wk; tta = tta + ttc; // то, что попадает слева направо
				tt = Ab[k] * w1*tt / wk; ttb = ttb + tt;
			} // то, что попадает справа налево
			q = q + 4;
		}
		sislr[j] = sig*tta; //собственное излучение стенок слева направо и справа налево
		sisrl[j] = sig*ttb;
	}
	slr = slr + Tra[kost - 1] * sislr[kost - 1];
	slr = slr + Ab[kost - 1] * sig*pow(Tss[kost - 1], 4e0);
	srl = srl + Tra[0] * sisrl[0];
	srl = srl + Ab[0] * sig*pow(Ts[0], 4e0);
	sislr[0] = 0.0; sisrl[kost - 1] = 0.0; //стенок вне многослойной стенки нет, излучение ниоткуда не падает, есть только внешнее излучение //ra=slr-srl; if (!vyb) cout << "slr = " << slr << "\tsrl = " << srl << "\tra = " << ra << "\tqc = " << qobsha << endl; //излучение за крайней правой стенкой и за крайней левой
	ao[p] = slr; p++; ao[p] = srl;
	k = 0; mu[k] = ao; k++; mu[k] = sislr; k++; mu[k] = sisrl; return mu;
}
double **RasLuchPloTepPot(int kost, double *hrl1, double *hlr1, double *ao, double *Tra, double *Ref, double *Ab, int sao, double *sislr, double *sisrl, double *Tsr, char *snm)
{
	int j, k = 3, saot = sao;
	double *Erlr = new double[kost], *Elrr = new double[kost], *Erl = new double[kost], *Elr = new double[kost], *Etem = new double[kost];
	double *Eres = new double[kost], *Eresl = new double[kost], *Eresr = new double[kost], **mu = new double*[k];
	double srl, slr, tth, tat, ttb, tta;
	if ((!Erlr) || (!Elrr) || (!Erl) || (!Elr) || (!Eresr) || (!Eresl) || (!Eres) || (!Etem)) { cout << snm << endl; k = getchar(); exit(1); }
	for (j = 0; j < kost; j++) { Erlr[j] = 0.0; Elrr[j] = 0.0; Erl[j] = 0.0; Elr[j] = 0.0; Eres[j] = 0.0; Eresl[j] = 0.0; Eresr[j] = 0.0; }
	srl = 0.0; srl = 0.0; tth = 0.0; tat = 0.0;
	for (k = 1; k < kost; k++) {
		srl = hrl1[k] * Tra[k]; //слева от k-ой стенки
		slr = hlr1[k]; tta = slr + srl*Ref[k - 1]; Elr[k] = Elr[k] + tta; //в воздушном промежутке между (k-1)-ой и k-ой стенками - излучение, идущее слева направо
		ttb = srl + slr*Ref[k]; Erl[k] = Erl[k] + ttb; //излучение, идущее справа налево
		Eresl[k] = Eresl[k] + tta - ttb; //вправо - влево
		slr = sislr[k]; srl = sisrl[k] * Tra[k]; tta = slr + srl*Ref[k - 1]; //слева - направо
		Elr[k] = Elr[k] + tta; ttb = srl + slr*Ref[k]; //справа - налево
		Erl[k] = Erl[k] + ttb; Eresl[k] = Eresl[k] + tta - ttb; //результирующее излучение слева от k-ой стенки
	} for (k = 0; k < (kost - 1); k++) {
		slr = hlr1[k] * Tra[k]; // справа от k-ой стенки
		srl = hrl1[k + 1]; tta = srl + slr*Ref[k + 1]; //влево
		Erlr[k] = Erlr[k] + tta; //справа от стенки, излучение идет справа налево
		ttb = srl*Ref[k] + slr; Elrr[k] = Elrr[k] + ttb; //вправо
		Eresr[k] = Eresr[k] + ttb - tta; //результирующее излучение справа от k-ой стенки
		slr = sislr[k] * Tra[k]; srl = sisrl[k + 1]; tta = slr + srl*Ref[k]; //вправо
		Elrr[k] = Elrr[k] + tta; ttb = srl + slr*Ref[k + 1]; //влево
		Erlr[k] = Erlr[k] + ttb; Eresr[k] = Eresr[k] + tta - ttb;
	}
	for (k = 0; k < ks; k++) Eres[k] = Eresl[k] - Eresr[k];
	ao = usrednen(ao, Tsr, saot, kost); saot++; //2
	ao = usrednen(ao, Eresl, saot, kost - 1); saot++; //3
	ao = usrednen(ao, Eresr, saot, kost - 1); saot++; //4
	ao = usrednen(ao, Eres, saot, kost - 1); saot++; //5
	ao = usrednen(ao, Erl, saot, kost - 1); saot++; //6
	ao = usrednen(ao, Elr, saot, kost - 1); ao[saot] = ao[saot] - ao[saot - 1]; saot++; //7
	ao = usrednen(ao, Elrr, saot, kost - 1); saot++; //8
	ao = usrednen(ao, Erlr, saot, kost - 1); ao[saot] = ao[saot - 1] - ao[saot]; saot++; //9
	ao[saot] = (ao[saot - 1] + ao[saot - 3]) / 2.0; //10 //vyvodfile(ao,saot,zf,hvsha,sfo);
	for (k=0; k<kost; k++) { if (!k) tta=Eresl[k]; else if (!(k+1-kost)) tta=Eresr[k]; else tta=(Eresr[k-1]+Eresl[k])/2e0; Etem[k]=tta; } 
	delete[]Eresl; delete[]Eresr; delete[]Erl; delete[]Elr; delete[]Elrr; delete[]Erlr;
	k = 0; mu[k] = ao; k++; mu[k] = Eres; k++; mu[k]=Etem; return mu;
}
double *opredTempStenShaFragm(double *Tna, int n, double *ktpvo, double *vte, double *ete, double *ktptk, int ce, int dmv, double htk, double hnf, double ptpob, double tenac, double zn)
{
	int k; double lamt, lamvo, lam, ht, tt; Tna[0] = tenac;
	for (k = 1; k < n; k++)  {
		tt = Tna[k - 1];
		lamt = opredKTPTKTochSha(ktptk, ete, tt, ce); //ищем КТП твердого скелета при данной температуре
		lamvo = opredKTPTKTochSha(ktpvo, vte, tt, dmv); //ищем КТП воздуха при данной температуре
		if (k % 2) { ht = htk; lam = lamt; }
		else { ht = hnf; lam = lamvo; }
		Tna[k] = tt + zn*ptpob*ht / lam;  //рассчитываем температуру на границах стенок 
	} return Tna;
}
double *opredPolTempTvKarShaFragm(double *Tna, int n, double *ete, double *ktptk, int ce, double htk, double hvp, double ptpob, double tenac, double zn)
{
	int k; double lamt, hh = hvp + htk, ht = htk + hvp / 2e0, tt; Tna[0] = tenac;
	for (k = 1; k < n; k++)  {
		tt = Tna[k - 1];
		lamt = opredKTPTKTochSha(ktptk, ete, tt, ce); //ищем КТП твердого скелета при данной температуре
		Tna[k] = tt + zn*ptpob*ht / lamt;  ht = ht + hh;
	} return Tna;
}  //рассчитываем температуру на границах стенок
double RasIzlSerStenNac(double *Refl, double *Tran, double *Abso, double *mpo, double *Reflt, double *Trant, double *Absot, double temc, double hlr0, double hrl0, double b, int vyve, double poristo, double dpct, double htch, double hvoz, int kost, double *ete, double *kttkve, int cemve, int vtve, double *qobve, double *ktpvove, double *tevove, char *snm, int dmao, double *stchsrve, double *Tpctve, int dmkvove)
{
	int j = 2, k, f = 6, q, *nost = new int[j];
	double p, r, hf = 1e0, ra, qs, d = 0.0, **mu, qobvee, *Elr, p1=hf, p2=hf, *Evsv; for (j = 0; j < kost; j++) d = d + hf; //cout << "d = " << d << endl;
	q = 2 * kost; double *tvs = new double[f], *tbs = new double[q], *epst = new double[f], *hvi = new double[f];
	double *Rs = new double[f], *Ap = new double[f], rsf = 0.0;
	if ((!tvs) || (!tbs) || (!epst) || (!hvi) || (!Rs) || (!Ap)) { cout << snm << endl; k = getchar(); exit(1); }
	else for (j = 0; j < f; j++) { tvs[j] = 0.0; epst[j] = 0.0; hvi[j] = 0.0; Rs[j] = 0.0; Ap[j] = 0.0; }
	qobvee = opredKTPTKTochSha(qobve, ete, temc, cemve);
	p = dpct; Tpctve = opredPolTempTvKarShaFragm(Tpctve, kost, ete, kttkve, cemve, htch, hvoz, p*qobvee, temc, -hf);
	tbs = opredTempStenShaFragm(tbs, q, ktpvove, tevove, ete, kttkve, cemve, dmkvove, htch, hvoz, (hf - p)*qobvee, temc, -hf);
	p1 = b; if ((1e0 - dpct)*b <= 1e0) p2 = 1e0; else if (pow(1e0 - dpct,2e0)*b <= 1e0) p2 = p1;
	double *hlr1 = NULL, *hrl1 = NULL, *reiz = NULL, slr, srl, h = hvoz, w = p2*h, l = p1*h;
	Ap[0] = w*l; Ap[1] = h*l; Ap[2] = Ap[0]; Ap[3] = Ap[1]; Ap[4] = w*h; Ap[5] = Ap[4]; nost[0] = 0; nost[1] = 2;
	double rhlr0 = hlr0 - hrl0, *Ts = NULL, *Tss = NULL, *Tsr = NULL, *Tna = NULL;
	qs = (1e0 - dpct)*qobvee; if (qs < rhlr0) { rhlr0 = qs; hlr0 = rhlr0 + hrl0; }
	mu = opredTempLPStenSha(Ts, Tss, Tsr, Tna, temc, q, ktpvove, tevove, ete, kttkve, cemve, dmkvove, htch, hvoz, qs - rhlr0, kost, snm);
	k = 0; Ts = mu[k]; k++; Tss = mu[k]; k++; Tsr = mu[k]; k++; Tna = mu[k]; delete[]mu;
	mu = RaschVneshIzluchSha(Tran, Refl, mpo, hlr1, hrl1, kost, hlr0, hrl0, snm);
	k = 0; hlr1 = mu[k]; k++; hrl1 = mu[k]; k++; reiz = mu[k]; delete[]mu;
	double *ao = new double[dmao], *sislr = NULL, *sisrl = NULL; k = 0; slr = reiz[k]; k++; srl = reiz[k];
	if (!ao) { cout << snm << endl; k = getchar(); exit(1); } for (k = 0; k < dmao; k++) ao[k] = 0.0;
	mu = RaschSobLuchPlotTepPot(ks, mpo, Ts, Tss, Tran, Refl, Abso, slr, srl, ao, 0, snm);
	k = 0; ao = mu[k]; k++; sislr = mu[k]; k++; sisrl = mu[k]; delete[]mu;
	mu = RasLuchPloTepPot(kost, hrl1, hlr1, ao, Tran, Refl, Abso, 2, sislr, sisrl, Tsr, snm);
	k = 0; ao = mu[k]; k++; Elr = mu[k]; k++; Evsv=mu[k]; delete[]mu;
	q = 0; rsf = 0.0; r = 0.0; for (j = 0; j < (kost - 1); j++) {
		for (k = 0; k < f; k++) tvs[k] = Tpctve[j]; tvs[0] = tbs[j]; tvs[2] = tbs[j + 1]; //for (k = 0; k < f; k++) Rs[k] = Reflt[j]; Rs[0] = Refl[j]; Rs[2] = Refl[j + 1]; 
		for (k = 0; k<f; k++) epst[k]=Absot[k]; epst[0]=Abso[j]+Tran[j]; epst[2]=Abso[j+1]+Tran[j+1]; //for (k = 0; k<f; k++) epst[k] = opredKTPTKTochSha(stchsrve, ete, tvs[k], cemve); epst[0] = Abso[j]; epst[2] = Abso[j + 1];
		for (k = 0; k<f; k++) Rs[k] = 1e0-epst[k]; //for (k = 0; k<f; k++) epst[k] = 1e0-Rs[k];
		for (k = 0; k<f; k++) hvi[k] = 0.0; hvi[0] = Elr[j]; //hvi[2]=-Elr[j];
		ra = SeryeStenkiRasIzl(w, h, l, tvs, epst, hvi, Rs, Ap, nost, f); //cout << "ra ( " << j << " ) = " << ra << endl; 
		if (fabs(ra)>0.0) { rsf = ra + rsf; r = r + hf; } q = q + 2; 
	}
	delete[]Ts; delete[]Tss; delete[]Tna; delete[]Tsr; delete[]hlr1;
	delete[]hrl1; delete[]reiz; delete[]sislr; delete[]sisrl;
	delete[]Elr; delete[]Evsv; delete[]ao; delete[]tvs; delete[]tbs; delete[]epst;
	delete[]hvi; delete[]Rs; delete[]nost; delete[]Ap; if (fabs(r)>0.0) rsf = rsf / r; else rsf = 0.0; return rsf;
}