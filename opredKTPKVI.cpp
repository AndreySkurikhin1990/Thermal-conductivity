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
const double tek0 = 273.15, tnosck = 3e2, dtosck = 1e2, tnack = 2e2, detek = 1e2;
const double y0kvi = 3e1*1e-3, tocraskvi = 1e-4, templak = 175e1 + tek0, dkoalk = 0.639201;
const double porkvi400=51.74*1e-2, porkvi500=52.07*1e-2, porkvi600=51.51*1e-2, porkvi700=39.75*1e-2;
const double porkvi800=40.85*1e-2, porkvi900=39.37*1e-2, porkvi1000=36.07*1e-2;
const int dsk = 60, ks = 10, dmkoosck = 14, dmkvokvi = 28, dmkok = 3, dkosckl = 6;
const int N = 13, vmik = 0, vybkvi = 6, nnxtk = 0, nxtk = 30; //vmik: 0 - керамовермикулит, 1 - ИТОМ, vybkvi: 4 - 400, 5 - 500, 6 - 600, 7 - 700, 8 - 800, 9 - 900, 10 - 1000
const char nfk = 'C';
struct derevo {
	int otre; //1 - отражение или 0 - пропускание
	int ste; //номер стенки
	int vis; //видимость: 1 - виден, 0 - нет
	int lev; //номер уровня
	struct derevo *back; //указатель назад
	struct derevo *next; //указатель вперед
};
class KoePog { public: double alp, tem; KoePog *nex; };
void RasKTPkvi();
void napstrkvi();
void initarrkvi(int);
double *koefoslab(double, double, double, double *, int, double *);
void NapMasVozdSha(double *, double *, int);
double epsisredkvi(double, double *, double *, int, double *, double *, int, int);
double *opredKTPTverKarkkvi(double *, double *, double, double, double, double, int, int, double, double, int, double *, double *, double *, int, int, int);
double **opredTempHolGor(double *, double *, int, int, double, int, double **, int, double *, double *);
double *NapMasKTPkvi(double *, int, int);
double **PoiskZavVelTemVer(int, double **, int, double *, double *, int);
double *koefPribSha(double *, double *, int, double *);
void zadrktkviNac();
void napMasEKTPkviNac(double, double, double);
double *vydelPolkvi(double *, double *, int, int);
double RasFracXeffkvi(int);
double *EffectTols(double *, double *, double *, double, double, int);
double RaschAlphaTvKarkvi();
double **rasPorpoRazkvi(int);
double BolTochRasAlpha(int, int, double, double, double, double *, double *, double *, int, int, char *, double *, double *, int, double *, double *, int, double *, double *, double *, double *, double *, double *);
double **RaschRTAkvi(int, double, double, double, int, int, double, int, double, int, int, int);
double *kopokvi(double, double, int, double, int, int, double, double, int, double *);
double **izmRTAkvi(double *, int, int, double *, double *, double *, double *, double *, double *, int);
double **chaRTAkvi(int, double *, double *, double *, double *, double *, double *, int);
double *oprRasTemNachkvi(int, int, double *, double *, int, double *, double, double, double, double, double, int);
double opredKTPTKTochSha(double *, double *, double, int);
double *zadrktkvi(int, int, double, int, double, double, int, int, double, int, int, double *);
double *izstNkvi(int, int, int, int);
double F0_lamT(double);
double *KoefPoglRosselNac(double *, int, int, double, double, double, double *, double *, int, double, double, double *, double *, int, int, int);
double **opredTempLPStenSha(double *, double *, double *, double *, double, int, double *, double *, double *, double *, int, int, double, double, double, int, char *);
double *FuncRaschIzl(double, double, double, double, int, double *, double, double *, double, double, double, double, int, double *, double *, int, int, char *, double, int, double *, double, double *, double *, int, double *, double *, int, double *, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, int, double, double, int, double *, double *, double, double, double, char *, int);
double *podschchieleSha(int, int, int, double *, double *);
double oprProcSoderpoPoris(double *, double *, double, int);
double opredLuchSostkvi(double *, double *, int, int, int);
void opredtemphckvi(double *, double *, double *, double *, double *, int, int, double);
double *opredKoefOtr(double *, double *, int, double, int, int, double *, double *, int);
double KorrZnachVozdProskviKon(double, double, double);
double *oprEffDoliTepPerenkvi(double, double, double);
double KorrZnachVozdProskvi(double, double, double, int);
double *opredTempStenShaFragm(double *, int, double *, double *, double *, double *, int, int, double, double, double, double, double);
void vyvodfile(double *, int, int, double, char *);
double bbfn(double);
double *poisMasKoefkvi(int);
void oprsodoxkvi();
void zapisvfile(double *, int, char *);
//-----------
char *sfatk = NULL, *sfnok = NULL, *snk = NULL, *ssk = NULL, *skptk = NULL, *svskv = NULL, *snmk = NULL;
char *sfok = NULL, *svfdku = NULL, *svfdk = NULL;
double *dkosckm = NULL, *dkosckt = NULL, *qobk = NULL, *etek = NULL, *kttkk = NULL, *mkok = NULL, *ektpk = NULL;
double *tgork = NULL, *tholk = NULL, *kektpk = NULL, *tsredkvi = NULL, *stchsrkvi = NULL, *Tkvi = NULL;
double *Akvi = NULL, *Rkvi = NULL, *Rakvi = NULL, *Takvi = NULL, tnak = tnack + tek0, *Aakvi = NULL;
double *Rtkvi = NULL, *Ttkvi = NULL, *Atkvi = NULL, *alphakvi = NULL, *Tpctk = NULL, porkvi, smgok, salok;
double ssiok, *txkvi = NULL, *qxkvi = NULL, *lxkvi = NULL, *kxkvi = NULL, *dtxkvi = NULL, *tkusck = NULL;
double *kusck = NULL, *ktpvokvi = NULL, *vtekvi = NULL, *ktrk = NULL, qobkvi, dkospk = 1e0, hkkvi, hkvi;
double hvokt, tmikvi, tmakvi, dpctk = 1e0, *ecktpkvi = NULL, tmak, tmik, cok, bok, *sodoxkvi;
double cnk = 1e-2, ckk = 1e2, bnk = 1e-2, bkk = 1e2, *temrask = NULL, tnrk, *ooxkvi = NULL, *mtsk = NULL;
int cemk = 11, vtk = 0, sctxk = 0, vtkk=3, vtkn = 0;
//------------
void zadrktkviNac()
{
	int j, jk = nxtk, jn = nnxtk, k, q, f = 6; initarrkvi(jk);
	double hf = 1e0, nf = 0.0, ka, kb, hnkvi = 0.0, dhk = 1e-3, hvh = 1e-6;
	double hvko = (13e1)*hvh, hvna, p=0.0, r, d = 0.0, *atr = NULL, t, ko, *po, **mu; 
	for (j = 0; j < ks; j++) d = d + hf; //cout << "d = " << d << "\t";
	int cek=0, u=0, w=0;
	double x = 0.0, *srra=NULL, *legr=NULL, *prgr=NULL;
	double *rapo=NULL, *rpr=NULL, srp=0.0, ce=0.0, cet=0.0, *pp=NULL;
	mu = rasPorpoRazkvi(vybkvi);
	k = 0; rpr = mu[k]; k++; rapo = mu[k]; k++; //0, 1
	srra=mu[k]; k++; prgr=mu[k]; k++; legr=mu[k]; k++; //2, 3, 4
	j=0; pp=mu[k]; srp = pp[j]; k++; if (pp) delete[]pp; //5
	pp=mu[k]; ce = pp[j]; k++; if (pp) delete[]pp; if (mu) delete[]mu; //6
	cet = ce; j = 0; while (cet > 0.0) { cet = cet - hf; j++; } cek = j; delete[]mu; //k=31; x = 0.0; for (j = 0; j < cek; j++) if (j <= k) x = x + srra[j] * rapo[j]; else break; //for (j=0; j<k; j++) cout << "sr_ra = " << srra[j] << "\tlegr = " << legr[j] << endl;  
	//--------------
	dkospk = RaschAlphaTvKarkvi(); vtkk=cemk;
	for (vtk = vtkn; vtk<vtkk; vtk++) { //пробегаем по температуре
		
		ko = srp;
		cout << "Sred razm por = " << ko << "\t";
		r = d; t = ko; ko = KorrZnachVozdProskvi(ko, r, porkvi, 0); cout << "Korr sred razm por = " << ko << "\t";
		r = d; ko = t; p = KorrZnachVozdProskvi(ko, r, porkvi, 1); cout << "Dol Plo = " << p << "\t";
		r = d; ko = t; ko = KorrZnachVozdProskviKon(ko, r, porkvi); hvna = ko; hvko = ko;
		if (vtk>vtkn) {
			q = jk - jn; ka = (tholk[vtk] - tgork[vtk]) / y0kvi; kb = tgork[vtk]; hkkvi = hnkvi;
			for (k = 0; k < q; k++) { temrask[k] = kb + ka*hkkvi; hkkvi = hkkvi + dhk; }
		}
		hvokt = hvna; hkvi = hvokt*(1e0 - porkvi) / porkvi;
		while (hvokt <= hvko) {
			j = jn; hkkvi = hnkvi*dhk; //пробегаем по размерам пор
			while ((j < jk) && (hkkvi < y0kvi)) { //cout << "hk = " << hkkvi << "\tte = " << etek[vtk] << endl; 
				sctxk = j; u=0; //пробегаем по координате
				mu = RaschRTAkvi(ks, hkvi, 0.0, 0.0, 1, vtk, hvokt, u, etek[vtk], u, u, 1);
				Rakvi = mu[u]; u++; Takvi = mu[u]; u++; Aakvi = mu[u]; u++; 
				Rkvi = mu[u];  u++; Tkvi = mu[u];  u++; Akvi = mu[u];  u++; 
				Rtkvi = mu[u]; u++; Ttkvi = mu[u]; u++; Atkvi = mu[u]; delete[]mu;
				u=0; w=1; atr = zadrktkvi(j, ks, d, vtk, hkvi, hvokt, w, u, etek[vtk], u, u, atr);
				kxkvi[j - jn] = hkkvi;
				hkkvi = hkkvi + dhk;
				j++; delete[]atr;
			}
			vyvodfile(lxkvi, jk - jn, 0, hvokt, sfok);
			for (j = 0; j<jk - jn - 1; j++) {
				p = kxkvi[j + 1] - kxkvi[j];
				r = txkvi[j] - txkvi[j + 1]; if (fabs(p)>0.0) r = fabs(r / p); else r = 0.0;
				p = fabs(qxkvi[j + 1] + qxkvi[j]) / 2e0; lxkvi[j] = p / r;
			}
			vyvodfile(ooxkvi, jk - jn, 2, hvokt, sfok);
			vyvodfile(txkvi, jk - jn, 2, hvokt, sfok);
			vyvodfile(lxkvi, jk - jn - 1, 2, hvokt, sfok);
			hvokt = hvokt + hvh;
		}
	}
	if (legr) delete[]legr; if (prgr) delete[]prgr; if (srra) delete[]srra;
	if (rpr) delete[]rpr; if (rapo) delete[]rapo;
}
void napstrkvi()
{
	if ((!sfatk) || (!sfok) || (!snk) || (!ssk) || (!skptk) || (!svskv) || (!snmk) || (!sfnok)) 
	{ cout << "No_memory!" << endl; getchar(); exit(1); }
	int k = 0;
	snk[k] = 'D'; k++; snk[k] = ':'; k++; snk[k] = '\\'; k++; snk[k] = '\\'; k++; snk[k] = '_';  k++; snk[k] = 'А';  k++;
	snk[k] = 'с'; k++; snk[k] = 'п'; k++; snk[k] = 'и';  k++; snk[k] = 'р';  k++; snk[k] = 'а';  k++; snk[k] = 'н';  k++;
	snk[k] = 'т'; k++; snk[k] = 'у'; k++; snk[k] = 'р';  k++; snk[k] = 'а';  k++; snk[k] = '\\'; k++; snk[k] = '\\'; k++;
	snk[k] = 't'; k++; snk[k] = 'm'; k++; snk[k] = 'p';  k++; snk[k] = '\\'; k++; snk[k] = '\\'; k++; snk[k] = '\0';
	k = 0;
	ssk[k] = 'C';  k++; ssk[k] = ':'; k++;  ssk[k] = '\\'; k++; ssk[k] = '\\'; k++; ssk[k] = 'U';  k++; ssk[k] = 's';  k++;
	ssk[k] = 'e';  k++; ssk[k] = 'r'; k++;  ssk[k] = 's';  k++; ssk[k] = '\\'; k++; ssk[k] = '\\'; k++; ssk[k] = 'А';  k++;
	ssk[k] = 'н';  k++; ssk[k] = 'д'; k++;  ssk[k] = 'р';  k++; ssk[k] = 'е';  k++; ssk[k] = 'й';  k++; ssk[k] = '\\'; k++;
	ssk[k] = '\\'; k++; ssk[k] = 'D'; k++;  ssk[k] = 'o';  k++; ssk[k] = 'c';  k++; ssk[k] = 'u';  k++; ssk[k] = 'm';  k++;
	ssk[k] = 'e';  k++; ssk[k] = 'n'; k++;  ssk[k] = 't';  k++; ssk[k] = 's';  k++; ssk[k] = '\\'; k++; ssk[k] = '\\'; k++;
	ssk[k] = '_';  k++; ssk[k] = 'А'; k++;  ssk[k] = 'с';  k++; ssk[k] = 'п';  k++; ssk[k] = 'и';  k++; ssk[k] = 'р';  k++;
	ssk[k] = 'а';  k++; ssk[k] = 'н'; k++;  ssk[k] = 'т';  k++; ssk[k] = 'у';  k++; ssk[k] = 'р';  k++; ssk[k] = 'а';  k++;
	ssk[k] = '\\'; k++; ssk[k] = '\\'; k++; ssk[k] = 't';  k++; ssk[k] = 'm';  k++; ssk[k] = 'p';  k++; ssk[k] = '\\'; k++;
	ssk[k] = '\\'; k++; ssk[k] = '\0';
	k = 0;
	skptk[k] = 'K';  k++; skptk[k] = 'o'; k++; skptk[k] = 'e'; k++; skptk[k] = 'f'; k++; skptk[k] = 'f'; k++; skptk[k] = 'i'; k++;
	skptk[k] = 'c';  k++; skptk[k] = 'i'; k++; skptk[k] = 'e'; k++; skptk[k] = 'n'; k++; skptk[k] = 't'; k++; skptk[k] = '_'; k++;
	skptk[k] = 'p';  k++; skptk[k] = 'o'; k++; skptk[k] = 'g'; k++; skptk[k] = 'l'; k++; skptk[k] = 'o'; k++; skptk[k] = 's'; k++;
	skptk[k] = 'c';  k++; skptk[k] = 'h'; k++; skptk[k] = 'e'; k++; skptk[k] = 'n'; k++; skptk[k] = 'i'; k++; skptk[k] = 'y'; k++;
	skptk[k] = 'a';  k++; skptk[k] = '_'; k++; skptk[k] = 'k'; k++; skptk[k] = 'v'; k++; skptk[k] = 'i'; k++; 
	if (vybkvi == 4) { skptk[k] = '4'; k++; skptk[k] = '0'; k++; skptk[k] = '0'; k++; }
	else if (vybkvi == 5) { skptk[k] = '5'; k++; skptk[k] = '0'; k++; skptk[k] = '0'; k++; }
	else if (vybkvi == 6) { skptk[k] = '6'; k++; skptk[k] = '0'; k++; skptk[k] = '0'; k++; }
	else if (vybkvi == 7) { skptk[k] = '7'; k++; skptk[k] = '0'; k++; skptk[k] = '0'; k++; }
	else if (vybkvi == 8) { skptk[k] = '8'; k++; skptk[k] = '0'; k++; skptk[k] = '0'; k++; }
	else if (vybkvi == 9) { skptk[k] = '9'; k++; skptk[k] = '0'; k++; skptk[k] = '0'; k++; }
	else if (vybkvi == 10) { skptk[k] = '1'; k++; skptk[k] = '0'; k++; skptk[k] = '0'; k++; skptk[k] = '0'; k++; }
	skptk[k] = '_'; k++; skptk[k] = 'T'; k++; skptk[k] = '.';  k++; skptk[k] = 't'; k++; 
	skptk[k] = 'x'; k++; skptk[k] = 't'; k++; skptk[k] = '\0'; k++;
	k = 0;
	svskv[k] = 'D'; k++; svskv[k] = 'o'; k++; svskv[k] = 'l'; k++; svskv[k] = 'i'; k++; 
	svskv[k] = '_'; k++; svskv[k] = 'p'; k++; svskv[k] = 'r'; k++; svskv[k] = 'o'; k++; 
	svskv[k] = 'p'; k++; svskv[k] = '_'; k++; svskv[k] = 'k'; k++; svskv[k] = 'v'; k++; 
	svskv[k] = 'i'; k++; svskv[k] = '-'; k++; svskv[k] = nfk; k++; svskv[k] = '.'; k++; 
	svskv[k] = 't'; k++; svskv[k] = 'x'; k++; svskv[k] = 't'; k++; svskv[k] = '\0';
	k = 0;
	snmk[k] = 'N'; k++; snmk[k] = 'o'; k++; snmk[k] = '_'; k++; snmk[k] = 'm'; k++; 
	snmk[k] = 'e'; k++; snmk[k] = 'm'; k++; snmk[k] = 'o'; k++; snmk[k] = 'r'; k++; 
	snmk[k] = 'y'; k++; snmk[k] = '!'; k++; snmk[k] = '\0';
	k = 0;
	sfnok[k] = 'F'; k++; sfnok[k] = 'i'; k++; sfnok[k] = 'l'; k++; sfnok[k] = 'e'; k++; 
	sfnok[k] = '_'; k++; sfnok[k] = 'i'; k++; sfnok[k] = 's'; k++; sfnok[k] = '_'; k++; 
	sfnok[k] = 'n'; k++; sfnok[k] = 'o'; k++; sfnok[k] = 't'; k++; sfnok[k] = '_'; k++;
	sfnok[k] = 'o'; k++; sfnok[k] = 'p'; k++; sfnok[k] = 'e'; k++; sfnok[k] = 'n'; k++; 
	sfnok[k] = '!'; k++; sfnok[k] = '\0';
	k = 0;
	svfdku[k] = 'V'; k++; svfdku[k] = 'y'; k++; svfdku[k] = 'v'; k++; svfdku[k] = 'o'; k++; 
	svfdku[k] = 'd'; k++; svfdku[k] = 'v'; k++; svfdku[k] = 'F'; k++; svfdku[k] = 'i'; k++; 
	svfdku[k] = 'l'; k++; svfdku[k] = 'e'; k++; svfdku[k] = '.'; k++; svfdku[k] = 't'; k++;
	svfdku[k] = 'x'; k++; svfdku[k] = 't'; k++; svfdku[k] = '\0';
	for (k = 0; k < (2 * dsk); k++) { sfatk[k] = '\0'; sfok[k] = '\0'; }
	strcpy(sfatk, ssk); strcat(sfatk, skptk); sfatk[strlen(sfatk) + 1] = '\0';
	strcpy(sfok, ssk); strcat(sfok, svskv); sfok[strlen(sfok) + 1] = '\0';
	strcpy(svfdk, ssk); strcat(svfdk, svfdku); svfdk[strlen(svfdk) + 1] = '\0';
}
void initarrkvi(int koel)
{
	double wmg, wsi, wal; int k, j;
	dkosckm = new double[dkosckl]; dkosckt = new double[dkosckl]; sodoxkvi=new double[6];
	if ((!dkosckm) || (!dkosckt) || (!sodoxkvi)) { cout << snmk << endl; k = getchar(); exit(1); }
	oprsodoxkvi(); k=0; wal=sodoxkvi[k]; k++; wsi=sodoxkvi[k]; k++; wmg=sodoxkvi[k]; delete []sodoxkvi; 
	snk = new char[dsk]; ssk = new char[dsk]; skptk = new char[dsk]; svskv = new char[dsk];
	snmk = new char[dsk]; sfnok = new char[dsk]; svfdku = new char[dsk];
	if ((!snk) || (!ssk) || (!skptk) || (!svskv) || (!snmk) || (!sfnok) || (!svfdku)) 
	{ cout << "No memory!" << endl; j = getchar(); exit(1); }
	for (j = 0; j < dsk; j++) { snk[j] = '\0'; ssk[j] = '\0'; skptk[j] = '\0'; svskv[j] = '\0'; snmk[j] = '\0'; 
	sfnok[j] = '\0'; svfdku[j] = '\0'; }
	j = 2 * dsk; sfatk = new char[j]; sfok = new char[j]; svfdk = new char[j];
	if ((!sfatk) || (!sfok) || (!svfdk)) { cout << "No memory!" << endl; j = getchar(); exit(1); }
	for (j = 0; j < (2 * dsk); j++) { sfatk[j] = '\0'; sfok[j] = '\0'; svfdk[j] = '\0'; }
	napstrkvi(); 
	txkvi = new double[koel]; qxkvi = new double[koel]; lxkvi = new double[koel]; kxkvi = new double[koel]; 
	dtxkvi = new double[koel]; ooxkvi = new double[koel];
	if ((!txkvi) || (!qxkvi) || (!lxkvi) || (!kxkvi) || (!dtxkvi) || (!ooxkvi)) { cout << snmk << endl; j = getchar(); exit(1); }
	for (j = 0; j < koel; j++) { txkvi[j] = 0.0; qxkvi[j] = 0.0; lxkvi[j] = 0.0; kxkvi[j] = 0.0; dtxkvi[j] = 0.0; ooxkvi[j] = 0.0; }
	qobk = new double[cemk]; etek = new double[cemk]; kttkk = new double[cemk]; 
	mkok = new double[cemk]; ektpk = new double[cemk]; tgork = new double[cemk]; 
	tholk = new double[cemk]; tsredkvi = new double[cemk]; stchsrkvi = new double[cemk];
	Tkvi = new double[ks]; Akvi = new double[ks]; Rkvi = new double[ks]; Rakvi = new double[ks];
	Takvi = new double[ks]; Aakvi = new double[ks]; Rtkvi = new double[ks]; Ttkvi = new double[ks];
	Atkvi = new double[ks]; alphakvi = new double[ks]; Tpctk = new double[ks]; 
	if ((!Tkvi) || (!Akvi) || (!Rkvi) || (!Rakvi) || (!Takvi) || (!Aakvi) || (!Ttkvi) || (!Rtkvi) || (!Atkvi) || (!alphakvi) || (!Tpctk)) 
	{ cout << snmk << endl; k = getchar(); exit(1); }
	if ((!qobk) || (!etek) || (!kttkk) || (!mkok) || (!ektpk) || (!tgork) || (!tholk) || (!stchsrkvi) || (!tsredkvi)) 
	{ cout << snmk << endl; k = getchar(); exit(1); }
	for (j = 0; j < cemk; j++) {
		qobk[j] = 0.0; etek[j] = 0.0; kttkk[j] = 0.0; mkok[j] = 0.0;
		ektpk[j] = 0.0; stchsrkvi[j] = 0.0; tgork[j] = 0.0; tholk[j] = 0.0; tsredkvi[j] = 0.0;
	} 
	for (j = 0; j < ks; j++) {
		Akvi[j] = 0.0; Rkvi[j] = 0.0; Rakvi[j] = 0.0; Takvi[j] = 0.0;
		Aakvi[j] = 0.0; Rtkvi[j] = 0.0; Ttkvi[j] = 0.0; Atkvi[j] = 0.0; Tkvi[j] = 0.0;
	}
	txkvi = new double[koel]; qxkvi = new double[koel]; lxkvi = new double[koel]; 
	kxkvi = new double[koel]; dtxkvi = new double[koel];
	if ((!txkvi) || (!qxkvi) || (!lxkvi) || (!kxkvi) || (!dtxkvi)) { cout << snmk << endl; j = getchar(); exit(1); }
	for (j = 0; j < koel; j++) { txkvi[j] = 0.0; qxkvi[j] = 0.0; lxkvi[j] = 0.0; kxkvi[j] = 0.0; dtxkvi[j] = 0.0; }
	j=0; etek[j] = tnak; for (j = 1; j < cemk; j++) etek[j] = etek[j - 1] + detek;
	tkusck = new double[dmkoosck]; kusck = new double[dmkoosck];
	if ((!tkusck) || (!kusck)) { cout << snmk << endl; j = getchar(); exit(1); }
	k=0; tkusck[k] = tnosck; for (k = 1; k < dmkoosck; k++) tkusck[k] = tkusck[k - 1] + dtosck; //cout << "wsi = " << wsi << "\twal = " << wal << "\twmg = " << wmg << endl;
	kusck = koefoslab(wmg, wsi, wal, tkusck, dmkoosck, kusck); 
	ktpvokvi = new double[dmkvokvi]; vtekvi = new double[dmkvokvi];
	for (j = 0; j < dmkvokvi; j++) { ktpvokvi[j] = 0.0; vtekvi[j] = 0.0; } //for (k = 0; k < dmkoosck; k++) cout << "tkusck = " << tkusck[k] << "\tkusck = " << kusck[k] << "\t"; cout << endl; //for (k = 0; k < dkosckl; k++) cout << "dkosckt = " << dkosckt[k] << "\tdkosckm = " << dkosckm[k] << "\t"; cout << endl;
	NapMasVozdSha(ktpvokvi, vtekvi, dmkvokvi); //for (k=0; k<dmkvokvi; k++) cout << "te = " << vtekvi[k] << "\tktp vozd ( " << k << " ) = " << ktpvokvi[k] << "\t";
	double s; for (k = 0; k < cemk; k++) { s = epsisredkvi(etek[k], tkusck, kusck, dmkoosck, dkosckt, dkosckm, dkosckl, vybkvi); stchsrkvi[k] = s; } //for (k=0; k<cemk; k++) cout << "Step cher ( " << k << " ) = " << stchsrkvi[k] << "\t"; cout << endl; 
	napMasEKTPkviNac(wmg, wsi, wal); //cout << "j = " << j << endl; 
}
double *poisMasKoefkvi(int vyb)
{
	int k = 3, j; double *kktp = new double[k], t1=25e0, t2=5e2, dt=t2-t1, kn;
	if (!kktp) { cout << snmk << endl; j = getchar(); exit(1); }
	kktp[2] = 0.0;
	if (vyb == 3) { kktp[1] = 0.00015; kktp[0] = 0.068; } //350
	else if (vyb == 4) { kktp[1] = 0.000125; kktp[0] = 0.082; //1
	kn=(0.156-0.087)/dt; kktp[1]=kn; kn=kn*t2; kktp[0] = 0.156-kn; //2
	kn=(0.155-0.087)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.155-kn; //3 
	kktp[1] = 0.00016; kktp[0] = 0.14; //4  //Spirina //выбор
	//kn=(0.146-0.087)/dt; kktp[1]=kn; kn=kn*t2; kktp[0] = 0.146-kn; //к вопросу о стандартизации КВИ
	} //400
	else if (vyb == 5) { kktp[1] = 0.0001; kktp[0] = 0.103; //из Ахтямова, 1991
	kn=(0.165-0.105)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.165-kn; //выбор
	//kn=(0.178-0.105)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.178-kn; //Двухслойные
	//kn=(0.152-0.105)/dt; kktp[1]=kn; kn=kn*t2; kktp[0] = 0.152-kn; //к вопросу о стандартизации КВИ
	} //500
	else if (vyb == 6) { kktp[1] = 0.00015; kktp[0] = 0.116; //выбор
	kn=(0.201-0.12)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.201-kn;
	kn=(0.195-0.12)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.195-kn;
	kktp[1] = 0.00015; kktp[0] = 0.17; //Spirina
	kn=(0.196-0.12)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.196-kn; //Двухслойные //выбор
	//kn=(0.178-0.12)/dt; kktp[1]=kn; kn=kn*t2; kktp[0] = 0.178-kn; //к вопросу о стандартизации КВИ
	} //600
	else if (vyb == 7) { kktp[1] = 0.00017; kktp[0] = 0.146; 
	kn=(0.216-0.15)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.216-kn; 
	kn=(0.235-0.15)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.235-kn; //к вопросу о стандартизации КВИ
	kn=(0.251-0.16)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.251-kn; //Двухслойные //выбор
	} //700
	else if (vyb == 8) { kktp[1] = 0.00018;  kktp[0] = 0.156; 
	kn=(0.226-0.16)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.226-kn; //выбор
	//kn=(0.25-0.16)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.25-kn;
	//kktp[1] = 0.00014; kktp[0] = 0.21; //Spirina
	//kn=(0.23-0.16)/dt; kktp[1]=kn; kn=kn*t2; kktp[0] = 0.23-kn; //к вопросу о стандартизации КВИ
	} //800
	else if (vyb == 9) { kktp[1] = 0.00019;  kktp[0] = 0.185; 
	kn=(0.23-0.195)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.23-kn; //выбор
	//kn=(0.29-0.195)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.29-kn;
	} //900
	else if (vyb == 10) { kktp[1] = 0.00025;  kktp[0] = 0.246; 
	kn=(0.287-0.25)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.287-kn; //выбор
	//kn=(0.36-0.25)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.36-kn;
	//kn=(0.35-0.25)/dt; kktp[1]=kn; kn=kn*t2; kktp[0] = 0.35-kn; //к вопросу о стандартизации КВИ
	} //1000
	return kktp;
}
void napMasEKTPkviNac(double wmg, double wsi, double wal)
{
	int vpv=vybkvi; 
	kektpk = poisMasKoefkvi(vybkvi);
	int k, nk, j; double hf = 1e0, g, s, t;
	double *ektpks = new double[cemk], *eteks = new double[cemk], **mu;
	if ((!ektpks) || (!eteks)) { cout << snmk << endl; k = getchar(); exit(1); }
	for (k = 0; k < cemk; k++) { t = etek[k]; eteks[k] = t; t=t-tek0; s = 0.0; g = 0.0; 
	for (j = 0; j<dmkok; j++) { s = s + kektpk[j] * pow(t, g); g = g + hf; } ektpks[k] = s; }
	ektpk = vydelPolkvi(ektpks, eteks, cemk, 0); delete[]etek; 
	etek = vydelPolkvi(ektpks, eteks, cemk, 1);
	double nf = 0.0, *p = vydelPolkvi(ektpks, eteks, cemk, 2), tn;
	nf = p[0]; delete[]p; k = 0; while (nf>0.0) { nf = nf - hf; k++; } 
	nk = k; cemk = nk; tn = etek[0]; delete[]ektpks; delete[]eteks;
	mu = opredTempHolGorkvi(ektpk, etek, cemk, dmkok, y0kvi);
	double *tgks = mu[0], *tkhks = mu[1], *qks = mu[2], *tsks = mu[3]; delete[]mu;
	etek[0] = tn; for (k = 1; k < cemk; k++) etek[k] = etek[k - 1] + detek; //for (k = 0; k < cemk; k++) cout << "t_g = " << tgks[k] << "\tt_h = " << tkhks[k] << "\tqks = " << qks[k] << "\ttsks = " << tsks[k] << "\tete = " << etek[k] << "\tektpk = " << ektpk[k] << endl;
	if (stchsrkvi) delete[]stchsrkvi; stchsrkvi = new double[cemk]; 
	if (!stchsrkvi) { cout << snmk << endl; k = getchar(); exit(1); }
	for (k = 0; k < cemk; k++) { g = epsisredkvi(etek[k], tkusck, kusck, dmkoosck, dkosckt, dkosckm, dkosckl, vybkvi); stchsrkvi[k] = g; } 
	for (k = 0; k < cemk; k++) { cout << "tem = " << etek[k] << "\tst_ch = " << stchsrkvi[k] << "\t"; cout << endl; } //zapisvfile(stchsrkvi, cemk, svfdk); //cout << "por kvi = " << porkvi << "\tce = " << cemk << endl; k=getchar();
	kttkk = opredKTPTverKarkkvi(etek, ektpk, porkvi, wsi, wal, wmg, vybkvi, dmkok, tnosck, dtosck, dmkoosck, kusck, tkusck, stchsrkvi, cemk, vmik, vpv);
}
double RasFracXeffkvi(int v)
{	int l = 2, k; double *mkbr = new double[l], *mv = new double[l], *tol = new double[l];
	k = 0; mkbr[k] = 230.078; k++; mkbr[k] = 231.006; k++;
	k = 0; mv[k] = 0.95; k++; mv[k] = 1.19; k++;
	k = 0; tol[k] = 0.64; k++; tol[k] = 0.64; k++;
	double rokbr = 2.75, rov = 0.44, *xv = EffectTols(mkbr, mv, tol, rov, rokbr, l), xsr = 0.0;
	if ((!v) || (v == 1)) xsr = xv[v];
	delete[]mkbr; delete[]mv; delete[]tol; delete[]xv; return xsr; }
double RaschAlphaTvKarkvi()
{	int l = 2, k=0, d = 30, q=0, cfk = l, cek, w=0;
	long j, jk = 1, pj = 10, h = 6;
	for (j = 0; j < h; j++) jk = jk * pj;
	double x = 0.0, y = 0.0, p = 0.0, xt = 0.0, yp = 0.0, hf = 1e0, epsilon=1e-9, e=1e-1;
	double pork=22.0, tc=pork+tek0, altc=1e-6, *srra=NULL, dvpn = 1448e0*altc / tc / 4e0, **mu=NULL;
	double *pn=new double[d], *alx=new double[d], *rcf=new double[cfk], sp=0.0;
	double *rapo=NULL, srp, ce, cet, *pp=NULL, *legr=NULL, *prgr=NULL, **pt=NULL, lamtem=0.0, *rpr; 
	mu = rasPorpoRazkvi(vybkvi); k=0; q=k; 
	rpr = mu[k]; k++; //0
	rapo = mu[k]; k++; //1
	srra = mu[k]; k++; //2
	prgr = mu[k]; k++; //3
	legr = mu[k]; //4
	pp = mu[k]; srp = pp[q]; k++; if (pp) delete[]pp; //4
	pp = mu[k]; ce = pp[q]; k++; if (pp) delete[]pp; if (mu) delete[]mu; //5
	cet = ce; j = 0; while (cet > e) { cet = cet - hf; j++; } cek = j; 
	x = 0.0; for (j = 0; j < cek; j++) if (j <= 21) x = x + srra[j] * rapo[j]; pork = x*porkvi / srp; 
	double pr = 0.0, rmf = 0.0, prf = 0.0, po = 0.0, *alsf = new double[cfk];
	for (j = 0; j < RAND_MAX; j++) rmf = rmf + hf;
	for (j = 0; j < jk; j++) po = po + hf;
	if ((pn) && (alx) && (rcf) && (alsf))
	{
		for (j = 0; j < d; j++) { pn[j] = 0.0; alx[j] = 0.0; }
		for (j = 0; j < cfk; j++) { rcf[j] = 0.0; alsf[j] = 0.0; }
	}
	else { cout << snmk << endl; j = getchar(); exit(1); }
	long lt; unsigned int st; lt = time(NULL); st = (unsigned int)(lt - (lt % 2)) / 2; srand(st);
		for (k = 0; k < l; k++) {
		x = RasFracXeffkvi(k); //размер частицы
		rcf[k] = x*1e-6; y = x*pork; //размер поры
		for (j = 0; j < d; j++) pn[j] = 0.0;
		for (j = 0; j < jk; j++) {
			pj = rand(); prf = 0.0; 
			for (h = 0; h < pj; h++) prf = prf + hf; pr = prf / rmf;
			yp = y*pr; xt = yp*(hf - pork) / pork;
			p = x / (xt + yp);
			pr = 0.0; for (h = 0; h < d; h++)	
			{ xt = pr + hf; if ((p >= pr) && (p < xt)) pn[h] = pn[h] + hf; pr = xt; } }
		pr = 0.0; for (j = 0; j < d; j++) { pn[j] = pn[j] / po; pr = pr + pn[j]; } //cout << "Summa = " << pr << endl; for (j = 0; j < d; j++) cout << "pn ( " << j << " ) = " << pn[j] << "\t";
		for (j = 0; j < d; j++) pn[j] = pn[j] / pr;
		for (j = 2; j < d; j++) {
			p = 0.0; for (h = 0; h<j; h++) p = p + hf;
			yp = pork * x * 1e-6 / (p - hf); xt = (hf - pork)*x*1e-6 / p; //cout << "x = " << x << "\txp = " << yp << "\txt = " << xt << endl;
			if (yp > dvpn) {
				lamtem = F0_lamT(yp*tc); if (lamtem<epsilon) lamtem = 0.0; 
				if (lamtem>hf) lamtem = hf; //cout << "f(lam*T)_1 = " << lamtem; 
				lamtem = bbfn(yp*tc*1e6); if (lamtem<epsilon) lamtem = 0.0; //cout << "\tf(lam*T)_2 = " << lamtem << endl; 
				if (lamtem>0.0) alx[j] = BolTochRasAlpha(0, j, yp, xt, tc, ktpvokvi, vtekvi, etek, cemk, dmkvokvi, snmk, Tkvi, Rkvi, 2 * N, kttkk, Akvi, 3, Rakvi, Takvi, Aakvi, Rtkvi, Ttkvi, Atkvi)*lamtem; else alx[j]=0.0; }
			else alx[j] = 0.0; 
		} w=0; pt = RaschRTAkvi(d, xt, 0.0, 0.0, 1, w, yp, 1, tc, w, w, w); pp = pt[w]; altc = pp[w]; delete[]pp; delete[]pt;
		alx[0] = 0.0; alx[1] = altc / (hf - pork); //for (j=0; j<d; j++) cout << "j = " << j << "\talx = " << alx[j] << "\t";
		sp=0.0; for (j = 1; j < d; j++) {
			p = 0.0; for (h = 0; h <= j; h++) p = p + hf;
			yp = pork*x*1e-6 / (p - hf); if (j>1) xt = oprProcSoderpoPoris(rapo, legr, yp, cek); else xt = 1e0; //cout << "j = " << j << "\txt = " << xt << "\t";
			alsf[k] = alsf[k] + pn[j] * alx[j] * xt; sp=sp+pn[j]*xt; } alsf[k]=alsf[k]/sp; }
	x = 0.0; yp = 0.0; for (j = 0; j < cfk; j++) { x = x + alsf[j]; yp = yp + hf; } x = x / yp; for (j=0; j<cfk; j++) cout << "j = " << j << "\tal_sr = " << alsf[j] << endl; 
	delete[]rapo; delete[]srra; delete[]rcf; delete[]pn; delete[]alx; delete[]alsf; delete[]legr;
	x = (x / altc); printf("sr_okp = %0.10lf\n", x); //ослабление КП за счет пористой структуры вермикулита
	return x; }
double **RaschRTAkvi(int kost, double htm, double kta, double ktb, int izm, int vyte, double hvo, int v, double ti, int c, int w, int u)
{	int k = kost, j, r; double *te = NULL, e = tocraskvi, **mauk=NULL; //cout << "c = " << c << "\tw = " << w << "\tu = " << u << endl;
	if (!v) { 
		alphakvi = kopokvi(kta, ktb, 2, 0.0, kost, vyte, htm, hvo, w, alphakvi);
		te = kopokvi(kta, ktb, 0, 0.0, kost, vyte, htm, hvo, w, alphakvi); //cout << "htm = " << htm << endl;
		if (!c) {
			k = kost; 
			if (Rakvi) delete[]Rakvi; Rakvi = new double[k];
			if (Rkvi) delete[]Rkvi; Rkvi = new double[k];
			if (Rtkvi) delete[]Rtkvi; Rtkvi = new double[k];
			if ((!Rakvi) || (!Rkvi) || (!Rtkvi)) { cout << snmk << endl; j = getchar(); exit(1); }
			k = kost; j=3; Rakvi = opredKoefOtr(te, Rakvi, k, tocraskvi, j, cemk, mkok, etek, vybkvi);
			for (j = 0; j < k; j++) { Rtkvi[j] = Rakvi[j]; Rkvi[j] = Rakvi[j]; }
			if (Takvi) delete[]Takvi; Takvi = new double[k];
			if (Ttkvi) delete[]Ttkvi; Ttkvi = new double[k];
			if (Tkvi) delete[]Tkvi; Tkvi = new double[k];
			if (Aakvi) delete[]Aakvi; Aakvi = new double[k];
			if (Atkvi) delete[]Atkvi; Atkvi = new double[k];
			if (Akvi) delete[]Akvi; Akvi = new double[k];
			if ((!Takvi) || (!Ttkvi) || (!Tkvi) || (!Aakvi) || (!Atkvi) || (!Akvi)) { cout << snmk << endl; j = getchar(); exit(1); }
			k = kost; for (j = 0; j < k; j++) {
				Takvi[j] = -alphakvi[j] * htm; Takvi[j] = exp(Takvi[j]); Aakvi[j] = 1e0 - Takvi[j] - Rakvi[j];
				if (Aakvi[j] < 0.0) { Aakvi[j] = alphakvi[j] * htm; Takvi[j] = 1e0 - Rakvi[j] - Aakvi[j]; }
				Atkvi[j] = Aakvi[j]; Ttkvi[j] = Takvi[j]; Akvi[j] = Aakvi[j]; Tkvi[j] = Takvi[j]; }
			k = kost; mauk = izmRTAkvi(te, k, 1, Rkvi, Tkvi, Akvi, Rakvi, Takvi, Aakvi, 0);
			r = 6; mauk[r] = Rtkvi; r++; mauk[r] = Ttkvi; r++; mauk[r] = Atkvi; }
		else if (c == 1) { r = 1; mauk = new double*[r]; if (!mauk) { cout << snmk << endl; j = getchar(); exit(1); } j=0; mauk[j] = te; } }
	else if (v == 1) { //cout << "v = " << v << endl;
		te = kopokvi(0.0, 0.0, 1, ti, kost, c, htm, hvo, w, alphakvi); r = 1;
		mauk = new double*[r]; if (!mauk) { cout << snmk << endl; j = getchar(); exit(1); } j=0; mauk[j] = te; }
	j=0; if (!u) mauk[j] = te; return mauk; }
double *kopokvi(double ta, double tb, int vyb, double tl, int kost, int vyte, double htk, double hvo, int w, double *al)
{ 	int k, q = 100, j, r, f;
	double t0 = 2e1, te = tek0 + t0, dte = 1e0, *p=NULL, e = 1e-3, t, a1, a2, t1, t2;
	char *s = new char[q]; KoePog *kp = new KoePog, *ne=NULL, *roo=NULL, *pre=NULL;
	if ((!kp) || (!s)) { cout << snmk << endl; k = getchar(); exit(1); } for (j = 0; j < q; j++) s[j] = '\0';
	ifstream fin; fin.open(sfatk); if (!fin.is_open()) { cout << sfnok << endl; k = getchar(); exit(1); }
	roo = kp; k = 0; while (!fin.eof()) {
		fin.getline(s, q, '\n'); ne = new KoePog; if (!ne) { cout << snmk << endl; j = getchar(); exit(1); }
		kp->alp = atof(s)*dkoalk*dkospk; kp->tem = te; kp->nex = ne; pre = kp; kp = ne; k++; te = te + dte; kp->nex = NULL;
	}
	if (ta>templak) ta=templak; if (tb>templak) tb=templak;
	delete[]ne; kp = NULL; pre->nex = kp; fin.close(); delete[]s; r = k;
	if ((!vyb) || (vyb == 2)) {
		f = 2; double *xi = new double[kost], *teks = new double[kost], *koe = new double[f], knat;
		if ((!teks) || (!koe) || (!xi)) { cout << snmk << endl; k = getchar(); exit(1); }
		teks = oprRasTemNachkvi(cemk, dmkok, teks, koe, kost, xi, htk, hvo, hkkvi, ta, tb, w);
		for (k = 0; k < kost; k++) {
			kp = roo->nex; pre = roo; j = 0; t = teks[k];
			while ((kp) && (j<r)) {
				if (kp->tem>t) {
					a1 = pre->alp; a2 = kp->alp;
					t1 = pre->tem; t2 = kp->tem; knat = (a2 - a1) / (t2 - t1); 
					al[k] = a1 + knat*(t - t1); break;	}
				else { pre = kp; kp = kp->nex; } j++;	} } delete[]xi; delete[]koe;
		if (vyb == 2) { p = al; delete[]teks; }
		else p = teks; }
	else if (vyb == 1) {
		kp = roo->nex; j = 0; pre = roo; t = 0.0;
		while ((kp) && (j<r)) {
			if (kp->tem>tl) {
				a1 = pre->alp; a2 = kp->alp; t1 = pre->tem; t2 = kp->tem;
				t = a1 + (a2 - a1)*(tl - t1) / (t2 - t1); break; }
			else { pre = kp; kp = kp->nex; } j++; }
		j = 1; p = new double[j]; j=0; if (p) p[j] = t; else { cout << snmk << endl; j = getchar(); exit(1); } }
	kp = roo; while (kp) { ne = kp->nex; delete kp; kp = ne; } //удаление списка
	return p; }
double **izmRTAkvi(double *tere, int kost, int izm, double *Ra, double *Ta, double *Aa, double *Rb, double *Tb, double *Ab, int v) //izm = 0 - нет изменений, izm - учитываются изменения //поиск изменения степени черноты или безразмерного коэффициента поглощения
{	double **mu, ko, dkosck; int k, rt = dmkoosck;
	for (k = 0; k < kost; k++) {
		ko = opredKTPTKTochSha(kusck, tkusck, tere[k], rt);
		dkosck = opredKTPTKTochSha(dkosckm, dkosckt, tere[k], dkosckl);
		if ((ko<0.0) || (ko>1e0) || (dkosck<0.0) || (dkosck>1e0) || (!izm)) { dkosck = 1e0; ko = 1e0; }
		Aa[k] = Aa[k] * ko*dkosck; Ta[k] = 1e0 - Aa[k] - Ra[k]; Tb[k] = Ta[k]; Ab[k] = Aa[k]; }
	mu = chaRTAkvi(kost, Ra, Ta, Aa, Rb, Tb, Ab, v);
	return mu; }
double **chaRTAkvi(int kost, double *Ra, double *Ta, double *Aa, double *Rb, double *Tb, double *Ab, int v)
{
	int f = 9, k = kost; double *tmp = new double[k], **mu = new double*[f]; 
	if ((!tmp) || (!mu)) { cout << snmk << endl; k = getchar(); exit(1); }
	for (k = 0; k < kost; k++) { if (Ta[k] * Ra[k] >= 1e0) v = 3; if (v == 3) { Ab[k] = Aa[k]; Rb[k] = Ra[k]; Tb[k] = Ta[k]; } }
	if (v < 3) {
		for (k = 0; k < kost; k++) {
			Ab[k] = (1e0 - Ta[k] + Ra[k] * Ta[k] - Ra[k]); Ab[k] = Ab[k] / (1e0 - Ta[k] * Ra[k]);
			tmp[k] = pow((1e0 - Ra[k])*Ta[k], 2e0)*Ra[k] / (1e0 - pow(Ra[k] * Ta[k], 2e0)) + Ra[k]; Rb[k] = tmp[k];
			Tb[k] = pow(1e0 - Ra[k], 2e0)*Ta[k] / (1e0 - pow(Ra[k] * Ta[k], 2e0)); } }
	delete[]tmp; k = 0; mu[k] = Ra; k++; mu[k] = Ta; k++; mu[k] = Aa; k++; 
	mu[k] = Rb; k++; mu[k] = Tb; k++; mu[k] = Ab; k++; return mu; }
double *oprRasTemNachkvi(int ce, int cee, double *teks, double *koe, int kost, double *xi, double htch, double hvozd, double hkokvi, double a, double b, int w)
{ 	double e = tocraskvi, hkx = 0.0, ht = 0.0; int k, j;
	if (!w) {
		k = 1; for (j = 0; j < ce; j++) if ((tgork[j] < e) || (tholk[j] < e) || (qobk[j] < e)) { k = 0; break; }
		if (!k) opredtemphckvi(ektpk, etek, tgork, tholk, qobk, ce, cee, y0kvi);
		hkx = hkokvi; koe[0] = (tholk[vtk] - tgork[vtk]) / y0kvi; 
		koe[1] = tgork[vtk]; tmak = tgork[vtk]; tmik = tholk[vtk]; }
	else if (w == 1) { hkx = y0kvi / 2e0; koe[0] = a; koe[1] = b; tmak = b; tmik = b - a*y0kvi; }
	xi[0] = hkx + htch / 2e0; ht = hvozd + htch;
	for (k = 1; k < kost; k++) xi[k] = xi[k - 1] + ht; //массив середин каждой из стенок по толщине
	for (k = 0; k < kost; k++) teks[k] = koe[0] * xi[k] + koe[1]; //for (k=0; k<kost; k++) cout << "k = " << k << "\txk = " << xi[k] << "\ttex = " << teks[k] << endl; //линеаризация поля температур
	qobkvi = opredKTPTKTochSha(qobk, etek, (teks[0] + teks[kost - 1]) / 2e0, ce); 
	if (qobkvi < 0.0) qobkvi = 0.0; return teks; }
void opredtemphckvi(double *efktpv, double *eftev, double *tgv, double *thv, double *qon, int dlma, int n, double h) //n=3 - длина массива коэффициентов приближающего многочлена, tgv - температура горячей стенки, thv - температура холодной стенки, dlma - длина массива ЭКТП
{ 	int k = 0, j; double *koeq = new double[n], *kho = new double[n], *kgo = new double[n];
	double *tsv = new double[n], t = tek0, ts, g, p, hf = 1e0, r, s;
	if ((!koeq) || (!kho) || (!kgo) || (!tsv)) { cout << snmk << endl; j = getchar(); exit(1); }
	for (k = 0; k < n; k++) tsv[k] = (tgv[k] + thv[k]) / 2e0;
	koeq = koefPribSha(qon, tsv, n, koeq);
	kgo = koefPribSha(tgv, tsv, n, kgo);
	kho = koefPribSha(thv, tsv, n, kho);
	for (k = 0; k < dlma; k++) {
		ts = eftev[k];
		g = 0.0; p = 0.0; 
		for (j = 0; j < n; j++) { g = g + pow(ts, p)*koeq[j]; p = p + hf; } 
		if (g<0.0) g = 0.0; if (qobk) qobk[k] = g;
		p = efktpv[k]; 
		if (fabs(p)>0.0) { g = fabs(g*h / p / 2e0); tgork[k] = eftev[k] + g; tholk[k] = eftev[k] - g; }
		else {
			r = 0.0; s = 0.0; for (j = 0; j < n; j++) { r = r + pow(ts, s)*kgo[j]; s = s + hf; } 
			if (r < 0.0) r = 0.0; if (tgork) tgork[k] = r;
			r = 0.0; s = 0.0; for (j = 0; j < n; j++) { r = r + pow(ts, s)*kho[j]; s = s + hf; } 
			if (r < 0.0) r = 0.0; if (tholk) tholk[k] = r; } }
delete[]tsv; delete[]koeq; delete[]kho; delete[]kgo; }
double *zadrktkvi(int zf, int kost, double d, int vyte, double htk, double hvo, int prod, int vy, double tc, int c, int u, double *rta)
{ 	int j; double **mu = RaschRTAkvi(kost, htk, 0.0, 0.0, 1, vyte, hvo, 0, tc, 1, 0, 0), *tt = mu[0], *te = new double[kost];
	if (!te) { cout << snmk << endl; j = getchar(); exit(1); }
	else for (j = 0; j < kost; j++) te[j] = tt[j];
	if (!u) {
		int k = 0, kst, q = 0, m = 4, b;
		double *prs = new double[m*kost*kost], *pr = new double[m], Er;
		if ((!prs) || (!pr)) { cout << snmk << endl; k = getchar(); exit(1); }
		for (k = 0; k < (m*kost*kost); k++) prs[k] = 0.0;
		q = 0; for (kst = 1; kst <= kost; kst++)
		{ 	for (b = 0; b < m; b++) pr[b] = 0.0;
			for (k = 1; k <= kost; k++) {
				pr = izstNkvi(k, kst, m, kost);
				for (b = 0; b < m; b++) prs[q + b] = pr[b]; q = q + m; } }
		if (prod == 1) Er = opredLuchSostkvi(prs, te, q, zf, vyte);
		delete[]pr; delete[]te; delete[]tt; delete[]mu; return prs; }
	else { delete[]tt; delete[]mu; return te; } }
double *izstNkvi(int izst, int kst, int l, int ocs)
{ 	int k; double *o;
	if (abs(izst - kst) <= N) 
		o = podschchieleSha(izst, kst, ocs, Rkvi, Tkvi);
	else {
		o = new double[l]; if (!o) { cout << snmk << endl; k = getchar(); exit(1); }
		for (k = 0; k < l; k++) o[k] = 0.0; }
	return o; }
double opredLuchSostkvi(double *prot, double *Tsr, int le, int zf, int vyte)
{
	int k = 0, j = 0, q = 0, m = 0, l = 3, mma = 30000, mm = 0, mt = 0, vy = 9, idr=3;
	double **mauk, *pt, s = 0.0, r = 0.0, p = 0.0, rez = 0.0, rap = 0.0;
	double *ao = new double[2 * N], *po = NULL, d = 0.0, hb = 5e-1, hb0 = hb, hf = 1e0; for (j = 0; j < ks; j++) d = d + hf;
	double de = hf, ta = 0.0, tb = 0.0, tc = 0.0, fa = 0.0, fb = 0.0, fc = 0.0, ba = bnk, bb = bkk, bc = (ba + bb) / 2e0, xk = 0.0, dxk = 1e-3;
	double ca = cnk, cb = ckk, cc = (ca + cb) / 2e0, hc = 1e-1, hc0 = hc, tt = 0.0;
	double mk = 1e2, mkm = mk, c00, b0 = bc, c0 = cc, kdc = 5e0, kdb = 5e0, tminkvi, tmaxkvi;
	if (!ao) { cout << snmk << endl; k = getchar(); exit(1); }
	for (k = 0; k < (2 * N); k++) ao[k] = 0.0;
	qobkvi = opredKTPTKTochSha(qobk, etek, etek[vtk], cemk); if (qobkvi < 0.0) qobkvi = 0.0;
	if ((vtk == vtkn) && (sctxk <= 2) && (!temrask)) {
		j = nxtk - nnxtk; temrask = new double[j];
		if (!temrask) { cout << snmk << endl; k = getchar(); exit(1); }
		ta = tholk[vtk] - tgork[vtk]; ta = ta / y0kvi; tb = tgork[vtk];
		for (k = 0; k < j; k++) { temrask[k] = tb + ta*xk; xk = xk + dxk; } }
	if ((vtk == vtkn) && (sctxk <= 2) && (!ktrk)) {
	k=3; ktrk = KoefPoglRosselNac(etek, idr, cemk, smgok, ssiok, salok, dkosckm, dkosckt, dkosckl, dkospk, dkoalk, kusck, tkusck, dmkoosck, vybkvi, k);
		for (k = 0; k < cemk; k++) cout << "ktr = " << ktrk[k] << "\t"; cout << endl; }
	l = vtk; tminkvi = tholk[l]; tmaxkvi = tgork[l];
	ta = tminkvi; tb = tmaxkvi; tc = temrask[sctxk];
	mauk = RaschRTAkvi(ks, hkvi, 0.0, 0.0, 1, vtk, hvokt, 0, tc, 0, 0, 1);
	k = 0; Rakvi = mauk[k]; k++; Takvi = mauk[k]; k++; Aakvi = mauk[k]; k++; Rkvi = mauk[k]; k++;
	Tkvi = mauk[k]; k++; Akvi = mauk[k]; k++; Rtkvi = mauk[k]; k++; Ttkvi = mauk[k]; k++; Atkvi = mauk[k]; delete[]mauk;
	ca = cnk; cb = ckk; m = 0; q = 0;
	while ((m<mma) && (fabs(cb - ca)>1e-3) && (!q)) {
		cc = (ca + cb) / 2e0;
		pt = FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, cc, vy, qobk, etek, cemk, 2, snmk, dpctk, vtk, ektpk, y0kvi, ktrk, kttkk, ks, ktpvokvi, vtekvi, dmkvokvi, ecktpkvi, hkvi, hvokt, Rkvi, Tkvi, Akvi, Rtkvi, Ttkvi, Atkvi, txkvi, qxkvi, dtxkvi, sctxk, tocraskvi, porkvi, 2 * N, stchsrkvi, Tpctk, b0, c0, mk, svfdk, 0);
		k = 3; fc = pt[k]; delete[]pt;
		if (fabs(fc) < 1e0) { c0 = cc; break; }
		pt = FuncRaschIzl(ta, tb, d, ta, m, prot, de, ao, rez, rap, bc, cc, vy, qobk, etek, cemk, 2, snmk, dpctk, vtk, ektpk, y0kvi, ktrk, kttkk, ks, ktpvokvi, vtekvi, dmkvokvi, ecktpkvi, hkvi, hvokt, Rkvi, Tkvi, Akvi, Rtkvi, Ttkvi, Atkvi, txkvi, qxkvi, dtxkvi, sctxk, tocraskvi, porkvi, 2 * N, stchsrkvi, Tpctk, b0, c0, mk, svfdk, 0);
		k = 3; fa = pt[k]; delete[]pt;
		if (fabs(fa) < 1e0) { c0 = ca; break; }
		pt = FuncRaschIzl(ta, tb, d, tb, m, prot, de, ao, rez, rap, bc, cc, vy, qobk, etek, cemk, 2, snmk, dpctk, vtk, ektpk, y0kvi, ktrk, kttkk, ks, ktpvokvi, vtekvi, dmkvokvi, ecktpkvi, hkvi, hvokt, Rkvi, Tkvi, Akvi, Rtkvi, Ttkvi, Atkvi, txkvi, qxkvi, dtxkvi, sctxk, tocraskvi, porkvi, 2 * N, stchsrkvi, Tpctk, b0, c0, mk, svfdk, 0);
		k = 3; fb = pt[k]; delete[]pt;
		if (fabs(fb) < 1e0) { c0 = cb; break; }
		if ((fa*fc) < 0.0) { cb = cc; q = 0; }
		if ((fb*fc) < 0.0) { ca = cc; q = 0; }
		if (((fb*fc)<0.0) && ((fa*fc)<0.0)) { cb = cc; q = 0; }
		if (((fb*fc)>0.0) && ((fa*fc)>0.0)) q = 1;
		m++; }
	c0 = cc; c00 = c0; cok = c00;
	m = 0; ba = bnk; bb = bkk; l = vtk; ta = tholk[l]; tb = tgork[l]; tc = temrask[sctxk];
	while ((m<mma) && ((bb - ba)>1e-1))
	{ 	bc = (ba + bb) / 2e0; b0 = bc; c0 = c00; cc = c0;
		pt = FuncRaschIzl(ta, tb, d, tb, m, prot, de, ao, rez, rap, bc, cc, vy, qobk, etek, cemk, 2, snmk, dpctk, vtk, ektpk, y0kvi, ktrk, kttkk, ks, ktpvokvi, vtekvi, dmkvokvi, ecktpkvi, hkvi, hvokt, Rkvi, Tkvi, Akvi, Rtkvi, Ttkvi, Atkvi, txkvi, qxkvi, dtxkvi, sctxk, tocraskvi, porkvi, 2 * N, stchsrkvi, Tpctk, b0, c0, mk, svfdk, 0);
		k = 9; s = pt[k]; k = 7; mk = pt[k]; k = 5; b0 = pt[k]; k = 4; tt = pt[k]; k = 10; r = pt[k]; delete[]pt;
		if ((mk<mkm) && (mk>0.0)) {
			mkm = mk; bok = bc; b0 = bc;
			lxkvi[sctxk] = r; ooxkvi[sctxk] = s; txkvi[sctxk] = tt; qxkvi[sctxk] = qobkvi; }
		m++; if (bc < 1e0) hb = hb0 / kdb; else hb = hb0; bb = bb - hb; }
	if (fabs(tnrk - tt)>1e-1)
		tnrk = tt; else {
			cnk = cok - hc0; if (cnk <= 0.0) cnk = 1e-2;
			ckk = cok + 2e0*hc0;
			bnk = bok - hb0; if (bnk <= 0.0) bnk = 1e-2;
			bkk = bok + 2e0*hb0;
			tnrk = tt;
		}	cout << "qo = " << qobkvi << "\ttc = " << tc << endl;
	delete[]ao; return rez; }
double KorrZnachVozdProskvi(double hps, double ksf, double por, int vy)
{	int j = 0, k = 1000; double pa = 1e-2, pb = 1e0, *po, pc, ra = fabs(pa - pb), e = tocraskvi, hf = 1e0;
	double fa, fb, fc, ta, tb, tc, tca, tcb, tcc, ka = hps*pa / por, kb = pb*hps / por, kc; //cout << "hps = " << hps << "\tksuf = " << ksf << "\tpor = " << por << endl;
	while ((ra > e) && (j < k)) { //подтягиваем пористость к значению, которое задали изначально, во время подстройки ЭКТП
		pc = (pa + pb) / 2e0;
		kc = hps*pc / por;
		po = oprEffDoliTepPerenkvi(kc, ksf, pc); tc = po[0]; delete[]po; //при 373 К
		tcc = kc*(hf - pc) / pc;
		fc = (hf - tc)*kc / (kc + tcc) - por;
		ka = hps*pa / por;
		po = oprEffDoliTepPerenkvi(ka, ksf, pa); ta = po[0]; delete[]po; //определяем долю площади сечения перемычки
		tca = ka*(hf - pa) / pa;
		fa = (hf - ta)*ka / (ka + tca) - por;
		kb = hps*pb / por;
		po = oprEffDoliTepPerenkvi(kb, ksf, pb); tb = po[0]; delete[]po; //через перемычку тепло распространяется чистой теплопроводностью
		tcb = kb*(hf - pb) / pb;
		fb = (hf - tb)*kb / (kb + tcb) - por;
		if ((fc*fb > 0.0) && (fa*fc<0.0)) pb = pc; if ((fc*fa>0.0) && (fb*fc < 0.0)) pa = pc;
		j++; ra = fabs(pa - pb); }
	dpctk = tc; //cout << "Dol Plo CTP = " << tc << endl;
	if (!vy) return kc; //скорректированное значение толщины воздушной прослойки (размер поры), когда ввели перемычку
	else if (vy == 1) return tc;
} //доля площади, через которую происходит перенос тепла чистой теплопроводностью
double *oprEffDoliTepPerenkvi(double ko, double d, double por)
{ 	int k, j, f, kost = ks, ksu = 2 * ks; double hvozd = ko, htch, hf = 1e0, e = tocraskvi;
	ecktpkvi = new double[kost]; if (!ecktpkvi) { cout << snmk << endl; k = getchar(); exit(1); }
	double *tepo = new double[ksu], *dol = new double[cemk], r, p, sa, sb, sc, t, fa, fb, fc;
	if ((!tepo) || (!dol)) { cout << snmk << endl; k = getchar(); exit(1); }
	for (j = 0; j < ksu; j++) tepo[j] = 0.0;
	hvozd = ko; htch = hvozd*(hf - por) / por;
	k = 1; for (j = 0; j < cemk; j++) if (qobk[j] < e) k = 0;
	if (!k) opredtemphckvi(ektpk, etek, tgork, tholk, qobk, cemk, dmkok, y0kvi); //cout << "hvoz = " << hvozd << "\tht = " << htch << endl;
	for (j = 0; j < cemk; j++) { //определяем ЭКТП многослойной стенки
		tepo = opredTempStenShaFragm(tepo, 2 * kost, ktpvokvi, vtekvi, etek, kttkk, cemk, dmkvokvi, htch, hvozd, qobk[j], etek[j], -hf); //в середине слоя
		r = 0.0; for (k = 0; k < ksu; k++) for (f = 0; f<ksu; f++) { p = tepo[k] - tepo[f]; if (p>r) r = p; }
		p = hvozd*(d - hf) + htch*d; r = r / p; t = qobk[j] / r; ecktpkvi[j] = t;
	} //for (j=0; j<cemi; j++) cout << "SKTPTK ( " << j << " ) = " << ecktpitom[j] << endl;
	delete[]tepo;
	f = 1000; for (j = 0; j<cemk; j++)
	{	sa = 0.0; sb = 1e0; k = 0;
		do { sc = (sa + sb) / 2e0;
			fa = kttkk[j] * sa + ecktpkvi[j] * (hf - sa) - ektpk[j]; //эффективные КТП многослойной стенки и  перемычки должны сравняться
			fb = kttkk[j] * sb + ecktpkvi[j] * (hf - sb) - ektpk[j]; //чтобы найти относительные доли площадей сечения переноса общей ПТП
			fc = kttkk[j] * sc + ecktpkvi[j] * (hf - sc) - ektpk[j];
			if ((fc*fb>0.0) && (fa*fc<0.0)) sb = sc;
			if ((fc*fa>0.0) && (fb*fc<0.0)) sa = sc;
			r = fabs(sa - sb); k++;
		} while ((r>e) && (k<f));
		dol[j] = sc; }
	return dol; }
double KorrZnachVozdProskviKon(double hps, double ksf, double por)
{ 	int j = 0, k = 1000; double pa = 1e-3, pb = 1e0, pc, ra = fabs(pa - pb);
	double fa, fb, fc, tca, tcb, tcc, ka = hps*pa / por, kb = pb*hps / por, kc, e = tocraskvi; //cout << "hps = " << hps << "\tksuf = " << ksf << "\tpor = " << por << endl;
	while ((ra>e) && (j < k)) { //подтягиваем пористость к значению, которое задали изначально, во время подстройки ЭКТП
		pc = (pa + pb) / 2e0;
		kc = hps*pc / por;
		tcc = kc*(1e0 - pc) / pc;
		fc = (1e0 - dpctk)*kc / (kc + tcc) - por;
		ka = hps*pa / por;
		tca = ka*(1e0 - pa) / pa;
		fa = (1e0 - dpctk)*ka / (ka + tca) - por;
		kb = hps*pb / por;
		tcb = kb*(1e0 - pb) / pb;
		fb = (1e0 - dpctk)*kb / (kb + tcb) - por;
		if ((fc*fb > 0.0) && (fa*fc<0.0)) pb = pc;
		if ((fc*fa>0.0) && (fb*fc < 0.0)) pa = pc;
		j++; ra = fabs(pa - pb);
	} //скорректированное значение толщины воздушной прослойки (размер поры), когда ввели перемычку
	return kc; }
void oprsodoxkvi()
{ int k=0; double wal, wsi, wmg, salok, smgok, ssiok, ko=1e-2;
double tnd = 6e2, dtd = 2e2, tm; dkosckt[0] = tnd; 
	for (k = 1; k < dkosckl; k++) dkosckt[k] = dkosckt[k - 1] + dtd;
	if (vybkvi == 4) {
		salok = 33e0; smgok = 15e0; ssiok = 52e0; porkvi = porkvi400; 
		wal = 25e0; wsi = 11e0; wmg = 4e1;
	} //КВИ-400
	else if (vybkvi == 5) {
		salok = 34e0; smgok = 11e0; ssiok = 54e0;
		wal = 28e0; wmg = 8e0; wsi = 44e0; porkvi = porkvi500;
	} //КВИ-500
	else if (vybkvi == 6) {
		salok = 36e0; smgok = 9e0; ssiok = 55e0;
		wal = 3e1; wsi = 7e0; wmg = 45e0; porkvi = porkvi600;
	} //КВИ-600
	else if (vybkvi == 7) {
		salok = 37e0; smgok = 8e0; ssiok = 55e0; 
		wal = 31e0; wmg = 6e0; wsi = 45e0; porkvi = porkvi700;
	} //КВИ-700
	else if (vybkvi == 8) {
		salok = 38e0; smgok = 7e0; ssiok = 55e0; 
		wal = 3e1; wmg = 5e0; wsi = 45e0; porkvi = porkvi800;
	} //КВИ-800
	else if (vybkvi == 9) {
		salok = 39e0; smgok = 6e0; ssiok = 55e0; 
		wal = 32e0; wmg = 5e0; wsi = 45e0; porkvi = porkvi900;
	} //КВИ-900
	else if (vybkvi == 10) {
		salok = 39e0; smgok = 6e0; ssiok = 55e0; 
		wal = 32e0; wmg = 4e0; wsi = 45e0; porkvi = porkvi1000;
	} //КВИ-1000
	else { cout << "Net takoy marki KVI!" << endl; k = getchar(); exit(1); }
	salok = salok*ko; smgok = smgok*ko; ssiok = ssiok*ko;
	wal = wal*ko; wsi = wsi*ko; wmg = wmg*ko;
		k = 0; dkosckm[k] = 3e0; k++; dkosckm[k] = 6e0; k++; dkosckm[k] = 7e0; k++; 
		dkosckm[k] = 15e0; k++; dkosckm[k] = 2e1; k++; dkosckm[k] = 26e0;
		if (vybkvi>6) { k=1; dkosckm[k]=7e0; } if (vybkvi>9) { k=2; dkosckm[k]=8e0; }
	for (k = 0; k < dkosckl; k++) {
		tm = dkosckm[k]*ko; dkosckm[k] = 1e0 - tm; //cout << "dko = " <<dkosckm[k] << "\t";
	}
	sodoxkvi[0]=wal; sodoxkvi[1]=wsi; sodoxkvi[2]=wmg; 
	sodoxkvi[3]=salok; sodoxkvi[4]=ssiok; sodoxkvi[5]=smgok; 
}