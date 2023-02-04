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
const double pori440 = (8e1+82e0)*1e-2/2e0, pori620 = (75e0+78e0)*1e-2/2e0, pori860 = (65e0+68e0)*1e-2/2, pori1000 = (62e0+65e0)*1e-2/2e0;
const double tei0 = 273.15, tnosci = 3e2, dtosci = 1e2, tnaci = 2e2;
const double detei = 1e2, y0itom = 3e1*1e-3, tocrasitom = 1e-4, templai = 1750.0 + tei0, dkoali = 0.639201;
const int dsi = 60, ks = 10, dmkoosci = 14, dmkvoitom = 28, dmkoi = 3;
const int N = 13, vtin = 0, nnxti = 0, nxti = 30, vybitom=3, vmik=1;
const char nfi = 'B';
struct derevo {
	int otre; //1 - îòðàæåíèå èëè 0 - ïðîïóñêàíèå
	int ste; //íîìåð ñòåíêè
	int vis; //âèäèìîñòü: 1 - âèäåí, 0 - íåò
	int lev; //íîìåð óðîâíÿ
	struct derevo *back; //óêàçàòåëü íàçàä
	struct derevo *next; //óêàçàòåëü âïåðåä
};
class KoePog { public: double alp, tem; KoePog *nex; };
void RasKTPItom();
double *poisMasKoefItom(int);
void napstritom();
void initarritom(int);
double *koefoslab(double, double, double, double *, int, double *);
void NapMasVozdSha(double *, double *, int);
double epsisreditom(double, double *, double *, int, double *, double *, int, int);
double *opredKTPTverKarkitom(double *, double *, double, double, double, double, int, int, double, double, int, double *, double *, double *, int, int, int);
double **opredTempHolGor(double *, double *, int, int, double, int, double **, int, double *, double *);
double *NapMasKTPItom(double *, int, int);
double **PoiskZavVelTemVer(int, double **, int, double *, double *, int);
double *koefPribSha(double *, double *, int, double *);
void zadrktITOMNac();
void napMasEKTPitomNac(double, double, double);
double *vydelPolitom(double *, double *, int, int);
double RasFracXeffitom63(int);
double *EffectTols(double *, double *, double *, double, double, int);
double RaschAlphaTvKaritom();
double **rasPorpoRazitom(int);
double BolTochRasAlpha(int, int, double, double, double, double *, double *, double *, int, int, char *, double *, double *, int, double *, double *, int, double *, double *, double *, double *, double *, double *);
double **RaschRTAitom(int, double, double, double, int, int, double, int, double, int, int, int);
double *kopoitom(double, double, int, double, int, int, double, double, int, double *);
double **izmRTAitom(double *, int, int, double *, double *, double *, double *, double *, double *, int);
double **chaRTAitom(int, double *, double *, double *, double *, double *, double *, int);
double *oprRasTemNachitom(int, int, double *, double *, int, double *, double, double, double, double, double, int);
double opredKTPTKTochSha(double *, double *, double, int);
double *zadrktitom(int, int, double, int, double, double, int, int, double, int, int, double *);
double *izstNitom(int, int, int, int);
double F0_lamT(double);
double *KoefPoglRosselNac(double *, int, int, double, double, double, double *, double *, int, double, double, double *, double *, int, int, int);
double **opredTempLPStenSha(double *, double *, double *, double *, double, int, double *, double *, double *, double *, int, int, double, double, double, int, char *);
double *FuncRaschIzl(double, double, double, double, int, double *, double, double *, double, double, double, double, int, double *, double *, int, int, char *, double, int, double *, double, double *, double *, int, double *, double *, int, double *, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, int, double, double, int, double *, double *, double, double, double, char *, int);
double *podschchieleSha(int, int, int, double *, double *);
double oprProcSoderpoPoris(double *, double *, double, int);
double opredLuchSostitom(double *, double *, int, int, int);
void opredtemphcitom(double *, double *, double *, double *, double *, int, int, double);
double *opredKoefOtr(double *, double *, int, double, int, int, double *, double *, int);
double KorrZnachVozdPrositomKon(double, double, double);
double *oprEffDoliTepPerenitom(double, double, double);
double KorrZnachVozdPrositom(double, double, double, int);
double *opredTempStenShaFragm(double *, int, double *, double *, double *, double *, int, int, double, double, double, double, double);
void vyvodfile(double *, int, int, double, char *);
double bbfn(double);
//-----------
char *sfati = NULL, *sfnoi = NULL, *sni = NULL, *ssi = NULL, *skpti = NULL, *svsi = NULL, *snmi = NULL;
char *sfoi = NULL, *svfdiu = NULL, *svfdi = NULL;
double *dkoscim = NULL, *dkoscit = NULL, *qobi = NULL, *etei = NULL, *kttki = NULL, *mkoi = NULL, *ektpi = NULL;
double *tgori = NULL, *tholi = NULL, *kektpi = NULL, *tsreditom = NULL, *stchsritom = NULL, *Titom = NULL;
double *Aitom = NULL, *Ritom = NULL, *Raitom = NULL, *Taitom = NULL, tnai = tnaci + tei0, *Aaitom = NULL;
double *Rtitom = NULL, *Ttitom = NULL, *Atitom = NULL, *alphaitom = NULL, *Tpcti = NULL, poritom, smgoi, saloi;
double ssioi, *txitom = NULL, *qxitom = NULL, *lxitom = NULL, *kxitom = NULL, *dtxitom = NULL, *tkusci = NULL;
double *kusci = NULL, *ktpvoitom = NULL, *vteitom = NULL, *ktri = NULL, qobitom, dkospi = 1e0, hkitom, hitom;
double hvoit, tmiitom, tmaitom, dpcti = 1e0, *ecktpitom = NULL, tmai, tmii, coi, boi;
double cni = 1e-2, cki = 1e2, bni = 1e-2, bki = 1e2, *temrasi = NULL, tnri, *ooxitom = NULL, *mtsi = NULL;
int dkoscil = 6, cemi = 11, vti = 0, sctxi = 0, vtik=cemi;
//------------
void zadrktITOMNac()
{
	int j, jk = nxti, jn = nnxti, k, q, f = 6, qk; initarritom(jk); 
	double hf = 1e0, nf = 0.0, ka, kb;
	double hnitom = 0.0, dhk = 1e-3, hvh = 1e-6, hvko = (13e1)*hvh, hvna=0.0, p=0.0, r=0.0, d = 0.0, *atr = NULL, t=0.0, ko=0.0, *po=NULL, **mu=NULL; 
	for (j = 0; j < ks; j++) d = d + hf; //cout << "d = " << d << "\t";
	//--------------
	int cei=0;
	double x = 0.0, *srra=NULL;
	double *rapo, *rpr, srp=0.0, ce=0.0, cet=0.0, *pp=NULL, *legr=NULL, *prgr=NULL;
	mu = rasPorpoRazitom(vybitom);
	k = 0; j=0; rapo = mu[k]; k++; //0
	srra = mu[k]; k++; //1
	prgr = mu[k]; k++; //2
	legr = mu[k]; k++; //3
	pp = mu[k]; srp=pp[j]; k++; if (pp) delete[]pp; //4
	pp = mu[k]; ce=pp[j]; k++; if (pp) delete[]pp; if (mu) delete[]mu; //5
	cet = ce; j = 0; while (cet > 0.0) { cet = cet - hf; j++; } cei = j; delete[]mu;
	k=31; x = 0.0; for (j = 0; j < cei; j++) if (j <= k) x = x + srra[j] * rapo[j]; else break; //for (j=0; j<k; j++) cout << "sr_ra = " << srra[j] << "\tlegr = " << legr[j] << endl;
	//--------------
	dkospi = RaschAlphaTvKaritom(); vtik=cemi;
	for (vti = vtin; vti<vtik; vti++) { //ïðîáåãàåì ïî òåìïåðàòóðå
		ko = srp; cout << "Sred razm por = " << ko << "\t";
		r = d; t = ko; ko = KorrZnachVozdPrositom(ko, r, poritom, 0); cout << "Korr sred razm por = " << ko << "\t";
		r = d; ko = t; p = KorrZnachVozdPrositom(ko, r, poritom, 1); cout << "Dol Plo = " << p << "\t";
		r = d; ko = t; ko = KorrZnachVozdPrositomKon(ko, r, poritom); hvna = ko; hvko = ko;
		if (vti>vtin) {
			q = jk - jn; ka = (tholi[vti] - tgori[vti]) / y0itom; kb = tgori[vti]; hkitom = hnitom;
			for (k = 0; k < q; k++) { temrasi[k] = kb + ka*hkitom; hkitom = hkitom + dhk; }
		}
		hvoit = hvna; hitom = hvoit*(1e0 - poritom) / poritom;
		while (hvoit <= hvko) {
			j = jn; hkitom = hnitom*dhk; //ïðîáåãàåì ïî ðàçìåðàì ïîð
			while ((j < jk) && (hkitom < y0itom)) {
				cout << "hk = " << hkitom << "\tte = " << etei[vti] << endl; sctxi = j; //ïðîáåãàåì ïî êîîðäèíàòå
				mu = RaschRTAitom(ks, hitom, 0.0, 0.0, 1, vti, hvoit, 0, etei[vti], 0, 0, 1);
				qk=0;	Raitom = mu[qk]; qk++; Taitom = mu[qk]; qk++; Aaitom = mu[qk]; qk++; 
						Ritom = mu[qk];  qk++; Titom = mu[qk];  qk++; Aitom = mu[qk];  qk++;
						Rtitom = mu[qk]; qk++; Ttitom = mu[qk]; qk++; Atitom = mu[qk]; delete[]mu;
				atr = zadrktitom(j, ks, d, vti, hitom, hvoit, 1, 0, etei[vti], 0, 0, atr);
				kxitom[j - jn] = hkitom;
				hkitom = hkitom + dhk;
				j++; delete[]atr;
			}
			vyvodfile(lxitom, jk - jn, 0, hvoit, sfoi);
			for (j = 0; j<jk - jn - 1; j++) {
				p = kxitom[j + 1] - kxitom[j];
				r = txitom[j] - txitom[j + 1]; if (fabs(p)>0.0) r = fabs(r / p); else r = 0.0;
				p = fabs(qxitom[j + 1] + qxitom[j]) / 2e0; lxitom[j] = p / r;
			}
			vyvodfile(ooxitom, jk - jn, 2, hvoit, sfoi);
			vyvodfile(txitom, jk - jn, 2, hvoit, sfoi);
			vyvodfile(lxitom, jk - jn - 1, 2, hvoit, sfoi);
			hvoit = hvoit + hvh;
		}
	}
}
void napstritom()
{
	if ((!sfati) || (!sfoi) || (!sni) || (!ssi) || (!skpti) || (!svsi) || (!snmi) || (!sfnoi)) { cout << "No_memory!" << endl; getchar(); exit(1); }
	int k = 0;
	sni[k] = 'D'; k++; sni[k] = ':';   k++; sni[k] = '\\'; k++; sni[k] = '\\'; k++; sni[k] = '_';  k++; sni[k] = 'À';  k++;
	sni[k] = 'ñ'; k++; sni[k] = 'ï';   k++; sni[k] = 'è';  k++; sni[k] = 'ð';  k++; sni[k] = 'à';  k++; sni[k] = 'í';  k++;
	sni[k] = 'ò';  k++; sni[k] = 'ó';  k++; sni[k] = 'ð';  k++; sni[k] = 'à';  k++; sni[k] = '\\'; k++; sni[k] = '\\'; k++;
	sni[k] = 't';  k++; sni[k] = 'm';  k++; sni[k] = 'p';  k++; sni[k] = '\\'; k++; sni[k] = '\\'; k++; sni[k] = '\0';
	k = 0;
	ssi[k] = 'C';  k++; ssi[k] = ':';  k++; ssi[k] = '\\'; k++; ssi[k] = '\\'; k++; ssi[k] = 'U';   k++; ssi[k] = 's';  k++;
	ssi[k] = 'e';  k++; ssi[k] = 'r';  k++; ssi[k] = 's';  k++; ssi[k] = '\\'; k++; ssi[k] = '\\';  k++; ssi[k] = 'À';  k++;
	ssi[k] = 'í';  k++; ssi[k] = 'ä';  k++; ssi[k] = 'ð';  k++; ssi[k] = 'å';  k++; ssi[k] = 'é';   k++; ssi[k] = '\\'; k++;
	ssi[k] = '\\'; k++; ssi[k] = 'D';  k++; ssi[k] = 'o';  k++; ssi[k] = 'c';  k++; ssi[k] = 'u';   k++; ssi[k] = 'm';  k++;
	ssi[k] = 'e';  k++; ssi[k] = 'n';  k++; ssi[k] = 't';  k++; ssi[k] = 's';  k++; ssi[k] = '\\';  k++; ssi[k] = '\\'; k++;
	ssi[k] = '_';  k++; ssi[k] = 'À';  k++; ssi[k] = 'ñ';  k++; ssi[k] = 'ï';  k++; ssi[k] = 'è';   k++; ssi[k] = 'ð';  k++;
	ssi[k] = 'à';  k++; ssi[k] = 'í';  k++; ssi[k] = 'ò';  k++; ssi[k] = 'ó';  k++; ssi[k] = 'ð';   k++; ssi[k] = 'à';  k++;
	ssi[k] = '\\'; k++; ssi[k] = '\\'; k++; ssi[k] = 't';  k++; ssi[k] = 'm';  k++; ssi[k] = 'p';   k++; ssi[k] = '\\'; k++;
	ssi[k] = '\\'; k++; ssi[k] = '\0';
	k = 0;
	skpti[k] = 'K';  k++; skpti[k] = 'o'; k++; skpti[k] = 'e'; k++; skpti[k] = 'f'; k++; skpti[k] = 'f'; k++; skpti[k] = 'i'; k++;
	skpti[k] = 'c';  k++; skpti[k] = 'i'; k++; skpti[k] = 'e'; k++; skpti[k] = 'n'; k++; skpti[k] = 't'; k++; skpti[k] = '_'; k++;
	skpti[k] = 'p';  k++; skpti[k] = 'o'; k++; skpti[k] = 'g'; k++; skpti[k] = 'l'; k++; skpti[k] = 'o'; k++; skpti[k] = 's'; k++;
	skpti[k] = 'c';  k++; skpti[k] = 'h'; k++; skpti[k] = 'e'; k++; skpti[k] = 'n'; k++; skpti[k] = 'i'; k++; skpti[k] = 'y'; k++;
	skpti[k] = 'a';  k++; skpti[k] = '_'; k++; skpti[k] = 'i'; k++; skpti[k] = 't'; k++; skpti[k] = 'o'; k++; skpti[k] = 'm'; k++;
	if (!vybitom) { skpti[k] = '4'; k++; skpti[k] = '4'; k++; skpti[k] = '0'; k++; }
	else if (vybitom == 1) { skpti[k] = '6'; k++; skpti[k] = '2'; k++; skpti[k] = '0'; k++; }
	else if (vybitom == 2) { skpti[k] = '8'; k++; skpti[k] = '6'; k++; skpti[k] = '0'; k++; }
	else if (vybitom == 3) { skpti[k] = '1'; k++; skpti[k] = '0'; k++; skpti[k] = '0'; k++; skpti[k] = '0'; k++; }
	skpti[k] = '_'; k++; skpti[k] = 'T'; k++; skpti[k] = '.'; k++; skpti[k] = 't'; k++; skpti[k] = 'x'; k++; skpti[k] = 't'; k++; skpti[k] = '\0'; k++;
	k = 0;
	svsi[k] = 'D'; k++; svsi[k] = 'o'; k++; svsi[k] = 'l'; k++; svsi[k] = 'i'; k++; svsi[k] = '_'; k++; svsi[k] = 'p'; k++; svsi[k] = 'r'; k++;
	svsi[k] = 'o'; k++; svsi[k] = 'p'; k++; svsi[k] = '_'; k++; svsi[k] = 'i'; k++; svsi[k] = 't'; k++; svsi[k] = 'o'; k++; svsi[k] = 'm'; k++;
	svsi[k] = '-'; k++; svsi[k] = nfi; k++; svsi[k] = '.'; k++; svsi[k] = 't'; k++; svsi[k] = 'x'; k++; svsi[k] = 't'; k++; svsi[k] = '\0';
	k = 0;
	snmi[k] = 'N'; k++; snmi[k] = 'o'; k++; snmi[k] = '_'; k++; snmi[k] = 'm'; k++; snmi[k] = 'e'; k++; snmi[k] = 'm'; k++; snmi[k] = 'o'; k++;
	snmi[k] = 'r'; k++; snmi[k] = 'y'; k++; snmi[k] = '!'; k++; snmi[k] = '\0';
	k = 0;
	sfnoi[k] = 'F'; k++; sfnoi[k] = 'i'; k++; sfnoi[k] = 'l'; k++; sfnoi[k] = 'e'; k++; sfnoi[k] = '_'; k++; sfnoi[k] = 'i'; k++;
	sfnoi[k] = 's'; k++; sfnoi[k] = '_'; k++; sfnoi[k] = 'n'; k++; sfnoi[k] = 'o'; k++; sfnoi[k] = 't'; k++; sfnoi[k] = '_'; k++;
	sfnoi[k] = 'o'; k++; sfnoi[k] = 'p'; k++; sfnoi[k] = 'e'; k++; sfnoi[k] = 'n'; k++; sfnoi[k] = '!'; k++; sfnoi[k] = '\0';
	k = 0;
	svfdiu[k] = 'V'; k++; svfdiu[k] = 'y'; k++; svfdiu[k] = 'v'; k++; svfdiu[k] = 'o'; k++; svfdiu[k] = 'd'; k++; svfdiu[k] = 'v'; k++;
	svfdiu[k] = 'F'; k++; svfdiu[k] = 'i'; k++; svfdiu[k] = 'l'; k++; svfdiu[k] = 'e'; k++; svfdiu[k] = '.'; k++; svfdiu[k] = 't'; k++;
	svfdiu[k] = 'x'; k++; svfdiu[k] = 't'; k++; svfdiu[k] = '\0';
	for (k = 0; k < (2 * dsi); k++) { sfati[k] = '\0'; sfoi[k] = '\0'; }
	strcpy(sfati, ssi); strcat(sfati, skpti); sfati[strlen(sfati) + 1] = '\0';
	strcpy(sfoi, ssi); strcat(sfoi, svsi); sfoi[strlen(sfoi) + 1] = '\0';
	strcpy(svfdi, ssi); strcat(svfdi, svfdiu); svfdi[strlen(svfdi) + 1] = '\0';
}
void initarritom(int koel)
{
	double wmg, wsi, wal, ko=1e-2; int k;
	dkoscim = new double[dkoscil]; dkoscit = new double[dkoscil];
	if ((!dkoscim) || (!dkoscit)) { cout << snmi << endl; k = getchar(); exit(1); }
	double tnd = 6e2, dtd = 2e2, tm; k=0; dkoscit[k] = tnd; for (k = 1; k < dkoscil; k++) dkoscit[k] = dkoscit[k - 1] + dtd;
	if (!vybitom) {
		saloi = 26e0*ko; smgoi = 22e0*ko; ssioi = 52e0*ko; poritom = pori440; //äëÿ òðåõêîìïîíåíòíîé ñìåñè
		wal = 23e0*ko; wmg = 19e0*ko; wsi = 49e0*ko; //äëÿ ìíîãîêîìïîíåíòíîé ñìåñè
		k = 0; dkoscim[k] = 6.07; k++; dkoscim[k] = 5.36; k++; dkoscim[k] = 6.19; k++; 
		dkoscim[k] = 13.48; k++; dkoscim[k] = 19.93; k++; dkoscim[k] = 27.69;
	} //ÈÒÎÌ-440
	else if (vybitom == 1) {
		saloi = 29e0*ko; smgoi = 16e0*ko; ssioi = 55e0*ko;
		wal = 26e0*ko; wmg = 15e0*ko; wsi = 5e1*ko; poritom = pori620;
		k = 0; dkoscim[k] = 6.41; k++; dkoscim[k] = 5.53; k++; dkoscim[k] = 6.32; k++; 
		dkoscim[k] = 13.56; k++; dkoscim[k] = 19.86; k++; dkoscim[k] = 27.66;
	} //ÈÒÎÌ-620
	else if (vybitom == 2) {
		saloi = 33e0*ko; smgoi = 11e0*ko; ssioi = 56e0*ko;
		wal = 3e1*ko; wmg = 1e1*ko; wsi = 52e0*ko; poritom = pori860;
		k = 0; dkoscim[k] = 7.28; k++; dkoscim[k] = 6.1; k++; dkoscim[k] = 6.83; k++; 
		dkoscim[k] = 14.02; k++; dkoscim[k] = 19.86; k++; dkoscim[k] = 27.22;
	} //ÈÒÎÌ-860
	else if (vybitom == 3) {
		saloi = 35e0*ko; smgoi = 9e0*ko; ssioi = 56e0*ko; 
		wal = 33e0*ko; wmg = 7e0*ko; wsi = 53e0*ko; poritom = pori1000;
		k = 0; dkoscim[k] = 7.45; k++; dkoscim[k] = 6.24; k++; dkoscim[k] = 6.95; k++; 
		dkoscim[k] = 14.15; k++; dkoscim[k] = 19.89; k++; dkoscim[k] = 27.09;
	} //ÈÒÎÌ-1000
	else { cout << "Net takoy marki ITOM! (initarr)" << endl; k = getchar(); exit(1); }
	for (k = 0; k < dkoscil; k++) {
		tm = dkoscim[k]*ko; dkoscim[k] = 1e0 - tm; //cout << "dko = " <<dkoscim[k] << "\t";
	} sni = new char[dsi]; ssi = new char[dsi]; skpti = new char[dsi]; svsi = new char[dsi];
	snmi = new char[dsi]; sfnoi = new char[dsi]; svfdiu = new char[dsi];
	int j; if ((!sni) || (!ssi) || (!skpti) || (!svsi) || (!snmi) || (!sfnoi) || (!svfdiu)) 
	{ cout << "No memory!" << endl; j = getchar(); exit(1); }
	for (j = 0; j < dsi; j++) { sni[j] = '\0'; ssi[j] = '\0'; skpti[j] = '\0'; svsi[j] = '\0'; snmi[j] = '\0'; 
	sfnoi[j] = '\0'; svfdiu[j] = '\0'; }
	j = 2 * dsi; sfati = new char[j]; sfoi = new char[j]; svfdi = new char[j];
	if ((!sfati) || (!sfoi) || (!svfdi)) { cout << "No memory!" << endl; j = getchar(); exit(1); }
	for (j = 0; j < (2 * dsi); j++) { sfati[j] = '\0'; sfoi[j] = '\0'; svfdi[j] = '\0'; }
	napstritom();
	txitom = new double[koel]; qxitom = new double[koel]; lxitom = new double[koel]; kxitom = new double[koel]; dtxitom = new double[koel]; ooxitom = new double[koel];
	if ((!txitom) || (!qxitom) || (!lxitom) || (!kxitom) || (!dtxitom) || (!ooxitom)) { cout << snmi << endl; j = getchar(); exit(1); }
	for (j = 0; j < koel; j++) { txitom[j] = 0.0; qxitom[j] = 0.0; lxitom[j] = 0.0; kxitom[j] = 0.0; dtxitom[j] = 0.0; ooxitom[j] = 0.0; }
	qobi = new double[cemi]; etei = new double[cemi]; kttki = new double[cemi]; mkoi = new double[cemi]; ektpi = new double[cemi];
	tgori = new double[cemi]; tholi = new double[cemi]; tsreditom = new double[cemi]; stchsritom = new double[cemi];
	Titom = new double[ks]; Aitom = new double[ks]; Ritom = new double[ks]; Raitom = new double[ks];
	Taitom = new double[ks]; Aaitom = new double[ks]; Rtitom = new double[ks]; Ttitom = new double[ks];
	Atitom = new double[ks]; alphaitom = new double[ks]; Tpcti = new double[ks];
	if ((!Titom) || (!Aitom) || (!Ritom) || (!Raitom) || (!Taitom) || (!Aaitom) || (!Ttitom) || (!Rtitom) || (!Atitom) || (!alphaitom) || (!Tpcti)) { cout << snmi << endl; j = getchar(); exit(1); }
	if ((!qobi) || (!etei) || (!kttki) || (!mkoi) || (!ektpi) || (!tgori) || (!tholi) || (!stchsritom) || (!tsreditom)) { cout << snmi << endl; j = getchar(); exit(1); }
	for (j = 0; j < cemi; j++) {
		qobi[j] = 0.0; etei[j] = 0.0; kttki[j] = 0.0; mkoi[j] = 0.0;
		ektpi[j] = 0.0; stchsritom[j] = 0.0; tgori[j] = 0.0; tholi[j] = 0.0; tsreditom[j] = 0.0;
	}
	for (j = 0; j < ks; j++) {
		Aitom[j] = 0.0; Ritom[j] = 0.0; Raitom[j] = 0.0; Taitom[j] = 0.0;
		Aaitom[j] = 0.0; Rtitom[j] = 0.0; Ttitom[j] = 0.0; Atitom[j] = 0.0; Titom[j] = 0.0;
	}
	txitom = new double[koel]; qxitom = new double[koel]; lxitom = new double[koel]; kxitom = new double[koel]; dtxitom = new double[koel];
	if ((!txitom) || (!qxitom) || (!lxitom) || (!kxitom) || (!dtxitom)) { cout << snmi << endl; j = getchar(); exit(1); }
	for (j = 0; j < koel; j++) { txitom[j] = 0.0; qxitom[j] = 0.0; lxitom[j] = 0.0; kxitom[j] = 0.0; dtxitom[j] = 0.0; }
	etei[0] = tnai; for (j = 1; j < cemi; j++) etei[j] = etei[j - 1] + detei;
	tkusci = new double[dmkoosci]; kusci = new double[dmkoosci];
	if ((!tkusci) || (!kusci)) { cout << snmi << endl; j = getchar(); exit(1); }
	k=0; tkusci[k] = tnosci; for (k = 1; k < dmkoosci; k++) tkusci[k] = tkusci[k - 1] + dtosci;
	kusci = koefoslab(wmg, wsi, wal, tkusci, dmkoosci, kusci);
	ktpvoitom = new double[dmkvoitom]; vteitom = new double[dmkvoitom];
	for (j = 0; j < dmkvoitom; j++) { ktpvoitom[j] = 0.0; vteitom[j] = 0.0; }
	NapMasVozdSha(ktpvoitom, vteitom, dmkvoitom); //for (k=0; k<dmkvoitom; k++) cout << "te = " << vteitom[k] << "\tktp vozd ( " << k << " ) = " << ktpvoitom[k] << "\t";
	double s; for (k = 0; k < cemi; k++) { s = epsisreditom(etei[k], tkusci, kusci, dmkoosci, dkoscit, dkoscim, dkoscil, vybitom); stchsritom[k] = s; } //for (k=0; k<cemi; k++) cout << "Step cher ( " << k << " ) = " << stchsritom[k] << "\t"; cout << endl; 
	napMasEKTPitomNac(wmg, wsi, wal);
}
void napMasEKTPitomNac(double wmg, double wsi, double wal)
{
	int vpv=vybitom; 
	kektpi = poisMasKoefItom(vybitom); 
	int k, nk, j, f=6; 
	double hf = 1e0, g=0.0, s=0.0, t=0.0, *p=NULL, nf = 0.0, tn=0.0;
	double *ektpi=new double[cemi], *eteis=new double[cemi], **mu=new double*[f];
	if ((!ektpi) || (!eteis) || (!mu)) { cout << snmi << endl; k = getchar(); exit(1); }
	for (k = 0; k < cemi; k++) { t = etei[k]; eteis[k] = t; s = 0.0; g = 0.0; 
	for (j = 0; j<dmkoi; j++) { s = s + kektpi[j] * pow(t, g); g = g + hf; } ektpi[k] = s; }
	mu = opredTempHolGor(ektpi, etei, cemi, dmkoi, y0itom, 0, mu, etei, qobi);
	mu=vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv, mu, f);
	k=0; ektpi=mu[k]; if (etei) delete[]etei; etei=mu[k]; k++; p=mu[k];
	k = 0; nf = p[k]; if (p) delete[]p; while (nf>0.0) { nf = nf - hf; k++; } nk = k; 
	k = 0; tn = etei[0]; cemi = nk; if (ektpi) delete[]ektpi; if (eteis) delete[]eteis;
	double *tgis = mu[0], *tkhis = mu[1], *qis = mu[2], *tsis = mu[3]; delete[]mu;
	etei[0] = tn; for (k = 1; k < cemi; k++) etei[k] = etei[k - 1] + detei; //for (k = 0; k < cemi; k++) cout << "t_g = " << tgis[k] << "\tt_h = " << tkhis[k] << "\tqis = " << qis[k] << "\ttsis = " << tsis[k] << "\tete = " << etei[k] << "\tektpi = " << ektpi[k] << endl;
	if (stchsritom) delete[]stchsritom; stchsritom = new double[cemi]; if (!stchsritom) { cout << snmi << endl; k = getchar(); exit(1); }
	for (k = 0; k < cemi; k++) { g = epsisreditom(etei[k], tkusci, kusci, dmkoosci, dkoscit, dkoscim, dkoscil, vybitom); stchsritom[k] = g; } 
	for (k = 0; k < cemi; k++) cout << "tem = " << etei[k] << "\tst_ch = " << stchsritom[k] << "\t"; cout << endl;
	kttki = opredKTPTverKarkitom(etei, ektpi, poritom, wsi, wal, wmg, vybitom, dmkoi, tnosci, dtosci, dmkoosci, kusci, tkusci, stchsritom, cemi, vmik, vpv); 
}
void RasKTPItom()
{
	int k, no = 0, f = 4;
	double *sv = new double[f], *kktpi, kk=1e-2;
	k=0; sv[k] = 6e1; k++; sv[k] = 5e1; k++; sv[k] = 25.0; k++; sv[k] = 2e1; 
	for (k = 0; k < f; k++) sv[k] = sv[k] * kk; //ñîäåðæàíèå âåðìèêóëèòà
	kktpi = poisMasKoefItom(no);
	delete[]sv;
}
double *poisMasKoefItom(int no)
{
	int f = 2, k; double *kti, t1=2e2, t2=38e1, kk=1e-2, te200=t1+tei0, te380=t2+tei0;
	if (!no) {
		double *kitom440 = new double[f], *ktpit = new double[f];
		if ((!kitom440) || (!ktpit)) { cout << snmi << endl; getchar(); exit(1); }
		k=0; ktpit[k] = 9e0*kk; k++; ktpit[k] = 12e0*kk;
		kitom440[1] = (ktpit[1] - ktpit[0]) / (te380 - te200); 
		kitom440[0] = ktpit[0] - kitom440[1] * te200;
		kti = kitom440; delete[]ktpit;
	}
	else if (no == 1) {
		double *kitom620 = new double[f], *ktpit = new double[f];
		if ((!kitom620) || (!ktpit)) { cout << snmi << endl; getchar(); exit(1); } 
		//ktpit[0] = 12.0*1e-2; ktpit[1] = 139.0*1e-3; //èç Äèññåðòàöèè
		k=0; ktpit[k] = 18.0*kk; k++; ktpit[k] = 19.0*kk; //Äàííûå 2017 ãîäà
		kitom620[1] = (ktpit[1] - ktpit[0]) / (te380 - te200); 
		kitom620[0] = ktpit[1] - kitom620[1] * te380;
		kti = kitom620; delete[]ktpit;
	}
	else if (no == 2) {
		double *kitom860 = new double[f], *ktpit = new double[f];
		if ((!kitom860) || (!ktpit)) { cout << snmi << endl; getchar(); exit(1); } 
		//ktpit[0] = 18.3*kk; ktpit[1] = 19.4*kk; //èç Äèññåðòàöèè
		ktpit[0] = 26.0*kk; ktpit[1] = 37.0*kk; //Äàííûå 2017 ãîäà
		kitom860[1] = (ktpit[1] - ktpit[0]) / (te380 - te200); kitom860[0] = ktpit[1] - kitom860[1] * te380;
		kti = kitom860; delete[]ktpit;
	}
	else if (no == 3) {
		double *kitom1000 = new double[f], *ktpit = new double[f];
		if ((!kitom1000) || (!ktpit)) { cout << snmi << endl; getchar(); exit(1); } 
		//ktpit[0] = 23.0*kk; ktpit[1] = 25.0*kk; //èç Äèññåðòàöèè
		ktpit[0] = 42.0*kk; ktpit[1] = 52.0*kk; //Äàííûå 2017 ãîäà
		kitom1000[1] = (ktpit[1] - ktpit[0]) / (te380 - te200); kitom1000[0] = ktpit[1] - kitom1000[1] * te380;
		kti = kitom1000; delete[]ktpit;
	}
	else { cout << "Net takoy marki ITOM!" << endl; k = getchar(); exit(1); }
	return kti;
}
double RasFracXeffitom63(int v)
{
	int l = 2, k; double *mkbr = new double[l], *mv = new double[l], *tol = new double[l];
	k = 0; mkbr[k] = 230.078; k++; mkbr[k] = 231.006; k++;
	k = 0; mv[k] = 0.95; k++; mv[k] = 1.19; k++;
	k = 0; tol[k] = 0.64; k++; tol[k] = 0.64; k++;
	double rokbr = 2.75, rov = 0.44, *xv = EffectTols(mkbr, mv, tol, rov, rokbr, l), xsr = 0.0;
	if ((!v) || (v==1)) xsr = xv[v]; 
	delete[]mkbr; delete[]mv; delete[]tol; delete[]xv; return xsr;
}
double RaschAlphaTvKaritom()
{
	int l = 2, k, d = 30, q, cfk = l, cei;
	long j, jk = 100 * 100 * 100, pj, h;
	double x = 0.0, y = 0.0, p = 0.0, xt = 0.0, yp = 0.0, hf = 1e0, tc = 22e0 + tei0;
	double e=1e-1, altc=0.0, *srra=NULL, pori, dvpn = 1448e0*1e-6/tc/4e0, **mu=NULL;
	double *pn = new double[d], *alx = new double[d], *rcf = new double[cfk], *rapo=NULL;
	double srp, ce, cet, *pp=NULL, *legr=NULL, **pt=NULL, lamtem, *prgr=NULL;
	mu=rasPorpoRazitom(vybitom);
	k=0; rapo=mu[k]; k++; srra=mu[k]; k++; prgr=mu[k]; k++; legr=mu[k]; k++; //0, 1, 2, 3
	k=0; q=0; pp=mu[k]; srp=pp[q]; k++; if (pp) delete[]pp; //4
	pp=mu[k]; ce=pp[q]; k++; if (pp) delete[]pp; if (mu) delete[]mu; //5
	cet=ce; j=0; while (cet>e) { cet = cet - hf; j++; } cei = j; delete[]mu;
	x = 0.0; for (j = 0; j < cei; j++) if (j <= 31) x = x + srra[j] * rapo[j]; pori = x*poritom / srp;
	double pr = 0.0, rmf = 0.0, prf = 0.0, po = 0.0, *alsf = new double[cfk];
	for (j = 0; j < RAND_MAX; j++) rmf = rmf + hf;
	for (j = 0; j < jk; j++) po = po + hf;
	if ((pn) && (alx) && (rcf) && (alsf))
	{
		for (j = 0; j < d; j++) { pn[j] = 0.0; alx[j] = 0.0; }
		for (j = 0; j < cfk; j++) { rcf[j] = 0.0; alsf[j] = 0.0; }
	}
	else { cout << snmi << endl; j = getchar(); exit(1); }
	long lt; unsigned int st; lt = time(NULL); st = (unsigned int)(lt - (lt % 2)) / 2; srand(st);
	for (k = 0; k < l; k++) {
		x = RasFracXeffitom63(k); //ðàçìåð ÷àñòèöû
		rcf[k] = x*1e-6; y = x*pori; //ðàçìåð ïîðû
		for (j = 0; j < d; j++) pn[j] = 0.0;
		for (j = 0; j < jk; j++) {
			pj = rand(); prf = 0.0; for (h = 0; h < pj; h++) prf = prf + hf; pr = prf / rmf;
			yp = y*pr; xt = yp*(hf - pori) / pori;
			p = x / (xt + yp);
			pr = 0.0; for (h = 0; h < d; h++)	{ xt = pr + hf; if ((p >= pr) && (p < xt)) pn[h] = pn[h] + hf; pr = xt; }
		}
		pr = 0.0; for (j = 0; j < d; j++) { pn[j] = pn[j] / po; pr = pr + pn[j]; } //cout << "Summa = " << pr << endl; for (j=0; j<d; j++) cout << "pn ( " << j << " ) = " << pn[j] << "\t";
		for (j = 0; j < d; j++) pn[j] = pn[j] / pr;
		for (j = 2; j < d; j++) {
			p = 0.0; for (h = 0; h<j; h++) p = p + hf;
			yp = pori*x*1e-6 / (p - hf); xt = (hf - pori)*x*1e-6 / p; //cout << "x = " << x << "\txp = " << yp << "\txt = " << xt << "\n";
			if (yp>dvpn) {
				lamtem = F0_lamT(yp*tc); if (lamtem<0.0) lamtem = 0.0; if (lamtem>hf) lamtem = hf; //cout << "f(lam*T)_1 = " << lamtem; 
				lamtem = bbfn(yp*tc*1e6);
				alx[j] = BolTochRasAlpha(0, j, yp, xt, tc, ktpvoitom, vteitom, etei, cemi, dmkvoitom, snmi, Titom, Ritom, 2 * N, kttki, Aitom, 2, Raitom, Taitom, Aaitom, Rtitom, Ttitom, Atitom)*lamtem;
			}
			else alx[j] = 0.0; //cout << "\tf(lam*T)_2 = " << lamtem << endl; 
		} pt = RaschRTAitom(d, xt, 0.0, 0.0, 1, 0, yp, 1, tc, 0, 0, 0); pp = pt[0]; altc = pp[0]; delete[]pp; delete[]pt;
		alx[0] = 0.0; alx[1] = altc / (hf - pori); //for (j=0; j<d; j++) cout << "j = " << j << "\talx = " << alx[j] << "\t";
		for (j = 1; j<d; j++) {
			p = 0.0; for (h = 0; h <= j; h++) p = p + hf;
			yp = pori*x*1e-6 / (p - hf); if (j>1) xt = oprProcSoderpoPoris(rapo, legr, yp, cei); else xt = 1e0; //cout << "j = " << j << "\txt = " << xt << "\t";
			alsf[k] = alsf[k] + pn[j] * alx[j] * xt;
		}
	}
	x = 0.0; yp = 0.0; for (j = 0; j < cfk; j++) { x = x + alsf[j]; yp = yp + hf; } 
	x = x / yp; for (j=0; j<cfk; j++) cout << "j = " << j << "\tal_sr = " << alsf[j] << endl; 
	delete[]rapo; delete[]srra; delete[]rcf; delete[]pn; delete[]alx; delete[]alsf; delete[]legr; delete[]prgr;
	x = (x / altc); printf("sr_okp = %0.20lf\n", x); //îñëàáëåíèå ÊÏ çà ñ÷åò ïîðèñòîé ñòðóêòóðû âåðìèêóëèòà
	return x;
}
double **RaschRTAitom(int kost, double htm, double kta, double ktb, int izm, int vyte, double hvo, int v, double ti, int c, int w, int u)
{
	int k = kost, j, r; double *te = NULL, e = tocrasitom, **mauk=NULL;
	if (!v) {
		alphaitom = kopoitom(kta, ktb, 2, 0.0, kost, vyte, htm, hvo, w, alphaitom);
		te = kopoitom(kta, ktb, 0, 0.0, kost, vyte, htm, hvo, w, alphaitom);
		if (!c) {
			k = kost;
			if (Raitom) delete[]Raitom; Raitom = new double[k];
			if (Ritom) delete[]Ritom; Ritom = new double[k];
			if (Rtitom) delete[]Rtitom; Rtitom = new double[k];
			k = kost; Raitom = opredKoefOtr(te, Raitom, k, tocrasitom, 2, cemi, mkoi, etei, vybitom);
			for (j = 0; j < k; j++) { Rtitom[j] = Raitom[j]; Ritom[j] = Raitom[j]; }
			if (Taitom) delete[]Taitom; Taitom = new double[k];
			if (Ttitom) delete[]Ttitom; Ttitom = new double[k];
			if (Titom) delete[]Titom; Titom = new double[k];
			if (Aaitom) delete[]Aaitom; Aaitom = new double[k];
			if (Atitom) delete[]Atitom; Atitom = new double[k];
			if (Aitom) delete[]Aitom; Aitom = new double[k];
			k = kost; for (j = 0; j < k; j++) {
				Taitom[j] = -alphaitom[j] * htm; Taitom[j] = exp(Taitom[j]); Aaitom[j] = 1e0 - Taitom[j] - Raitom[j];
				if (Aaitom[j] < 0.0) { Aaitom[j] = alphaitom[j] * htm; Taitom[j] = 1e0 - Raitom[j] - Aaitom[j]; }
				Atitom[j] = Aaitom[j]; Ttitom[j] = Taitom[j]; Aitom[j] = Aaitom[j]; Titom[j] = Taitom[j];
			}
			k = kost; mauk = izmRTAitom(te, k, 1, Ritom, Titom, Aitom, Raitom, Taitom, Aaitom, 0);
			r = 6; mauk[r] = Rtitom; r++; mauk[r] = Ttitom; r++; mauk[r] = Atitom;
		}
		else if (c == 1) { r = 1; mauk = new double*[r]; if (!mauk) { cout << snmi << endl; j = getchar(); exit(1); } mauk[0] = te; }
	}
	else if (v == 1) {
		te = kopoitom(0.0, 0.0, 1, ti, kost, c, htm, hvo, w, alphaitom); r = 1;
		mauk = new double*[r]; if (!mauk) { cout << snmi << endl; j = getchar(); exit(1); } mauk[0] = te;
	}
	if (!u) mauk[0] = te; return mauk;
}
double *kopoitom(double ta, double tb, int vyb, double tl, int kost, int vyte, double htk, double hvo, int w, double *al)
{
	int k, q = 100, j, r, f;
	double t0 = 2e1, te = tei0 + t0, dte = 1e0, *p=NULL, e = 1e-3, t, a1, a2, t1, t2;
	char *s = new char[q]; KoePog *kp = new KoePog, *ne=NULL, *roo=NULL, *pre=NULL;
	if ((!kp) || (!s)) { cout << snmi << endl; k = getchar(); exit(1); } for (j = 0; j < q; j++) s[j] = '\0';
	ifstream fin; fin.open(sfati); if (!fin.is_open()) { cout << sfnoi << endl; k = getchar(); exit(1); }
	roo = kp; k = 0; while (!fin.eof()) {
		fin.getline(s, q, '\n'); ne = new KoePog; if (!ne) { cout << snmi << endl; j = getchar(); exit(1); }
		kp->alp = atof(s)*dkoali*dkospi; kp->tem = te; kp->nex = ne; pre = kp; kp = ne; k++; te = te + dte; kp->nex = NULL;
	}
	delete[]ne; kp = NULL; pre->nex = kp; fin.close(); delete[]s; r = k;
	if ((!vyb) || (vyb == 2)) {
		f = 2; double *xi = new double[kost], *teks = new double[kost], *koe = new double[f], knat;
		if ((!teks) || (!koe) || (!xi)) { cout << snmi << endl; k = getchar(); exit(1); }
		teks = oprRasTemNachitom(cemi, dmkoi, teks, koe, kost, xi, htk, hvo, hkitom, ta, tb, w);
		for (k = 0; k < kost; k++) {
			kp = roo->nex; pre = roo; j = 0; t = teks[k];
			while ((kp) && (j<r)) {
				if (kp->tem>t) {
					a1 = pre->alp; a2 = kp->alp;
					t1 = pre->tem; t2 = kp->tem; knat = (a2 - a1) / (t2 - t1); al[k] = a1 + knat*(t - t1); break;
				}
				else { pre = kp; kp = kp->nex; } j++;
			}
		} delete[]xi; delete[]koe;
		if (vyb == 2) { p = al; delete[]teks; }
		else p = teks;
	}
	else if (vyb == 1) {
		kp = roo->nex; j = 0; pre = roo; t = 0.0;
		while ((kp) && (j<r)) {
			if (kp->tem>tl)
			{
				a1 = pre->alp; a2 = kp->alp; t1 = pre->tem; t2 = kp->tem;
				t = a1 + (a2 - a1)*(tl - t1) / (t2 - t1); break;
			}
			else { pre = kp; kp = kp->nex; } j++;
		}
		j = 1; p = new double[j]; if (p) p[0] = t; else { cout << snmi << endl; j = getchar(); exit(1); }
	}
	kp = roo; while (kp) { ne = kp->nex; delete kp; kp = ne; } //óäàëåíèå ñïèñêà
	return p;
}
double **izmRTAitom(double *tere, int kost, int izm, double *Ra, double *Ta, double *Aa, double *Rb, double *Tb, double *Ab, int v) //izm = 0 - íåò èçìåíåíèé, izm - ó÷èòûâàþòñÿ èçìåíåíèÿ //ïîèñê èçìåíåíèÿ ñòåïåíè ÷åðíîòû èëè áåçðàçìåðíîãî êîýôôèöèåíòà ïîãëîùåíèÿ
{
	double **mu, ko, dkosci; int k, rt = dmkoosci;
	for (k = 0; k < kost; k++) {
		ko = opredKTPTKTochSha(kusci, tkusci, tere[k], rt);
		dkosci = opredKTPTKTochSha(dkoscim, dkoscit, tere[k], dkoscil);
		if ((ko<0.0) || (ko>1e0) || (dkosci<0.0) || (dkosci>1e0) || (!izm)) { dkosci = 1e0; ko = 1e0; }
		Aa[k] = Aa[k] * ko*dkosci; Ta[k] = 1e0 - Aa[k] - Ra[k]; Tb[k] = Ta[k]; Ab[k] = Aa[k];
	}
	mu = chaRTAitom(kost, Ra, Ta, Aa, Rb, Tb, Ab, v);
	return mu;
}
double **chaRTAitom(int kost, double *Ra, double *Ta, double *Aa, double *Rb, double *Tb, double *Ab, int v)
{
	int f = 9, k = kost; double *tmp = new double[k], **mu = new double*[f]; if ((!tmp) || (!mu)) { cout << snmi << endl; k = getchar(); exit(1); }
	for (k = 0; k < kost; k++) { if (Ta[k] * Ra[k] >= 1e0) v = 3; if (v == 3) { Ab[k] = Aa[k]; Rb[k] = Ra[k]; Tb[k] = Ta[k]; } }
	if (v < 3) {
		for (k = 0; k < kost; k++) {
			Ab[k] = (1e0 - Ta[k] + Ra[k] * Ta[k] - Ra[k]); Ab[k] = Ab[k] / (1e0 - Ta[k] * Ra[k]);
			tmp[k] = pow((1e0 - Ra[k])*Ta[k], 2e0)*Ra[k] / (1e0 - pow(Ra[k] * Ta[k], 2e0)) + Ra[k]; Rb[k] = tmp[k];
			Tb[k] = pow(1e0 - Ra[k], 2e0)*Ta[k] / (1e0 - pow(Ra[k] * Ta[k], 2e0));
		}
	}
	delete[]tmp; k = 0; mu[k] = Ra; k++; mu[k] = Ta; k++; mu[k] = Aa; k++; mu[k] = Rb; k++; mu[k] = Tb; k++; mu[k] = Ab; k++; return mu;
}
double *oprRasTemNachitom(int ce, int cee, double *teks, double *koe, int kost, double *xi, double htch, double hvozd, double hkoitom, double a, double b, int w)
{
	double e = tocrasitom, hkx = 0.0, ht = 0.0; int k, j;
	if (!w) {
		k = 1; for (j = 0; j < ce; j++) if ((tgori[j] < e) || (tholi[j] < e) || (qobi[j] < e)) { k = 0; break; }
		if (!k) opredtemphcitom(ektpi, etei, tgori, tholi, qobi, ce, cee, y0itom);
		hkx = hkoitom; koe[0] = (tholi[vti] - tgori[vti]) / y0itom; koe[1] = tgori[vti]; tmai = tgori[vti]; tmii = tholi[vti];
	}
	else if (w == 1) { hkx = y0itom / 2e0; koe[0] = a; koe[1] = b; tmai = b; tmii = b - a*y0itom; }
	xi[0] = hkx + htch / 2e0; ht = hvozd + htch;
	for (k = 1; k < kost; k++) xi[k] = xi[k - 1] + ht; //ìàññèâ ñåðåäèí êàæäîé èç ñòåíîê ïî òîëùèíå
	for (k = 0; k < kost; k++) teks[k] = koe[0] * xi[k] + koe[1]; //for (k=0; k<kost; k++) cout << "k = " << k << "\txk = " << xi[k] << "\ttex = " << teks[k] << endl; //ëèíåàðèçàöèÿ ïîëÿ òåìïåðàòóð
	qobitom = opredKTPTKTochSha(qobi, etei, (teks[0] + teks[kost - 1]) / 2e0, ce); if (qobitom < 0.0) qobitom = 0.0; return teks;
}
void opredtemphcitom(double *efktpv, double *eftev, double *tgv, double *thv, double *qon, int dlma, int n, double h) //n=3 - äëèíà ìàññèâà êîýôôèöèåíòîâ ïðèáëèæàþùåãî ìíîãî÷ëåíà, tgv - òåìïåðàòóðà ãîðÿ÷åé ñòåíêè, thv - òåìïåðàòóðà õîëîäíîé ñòåíêè, dlma - äëèíà ìàññèâà ÝÊÒÏ
{
	int k = 0, j; double *koeq = new double[n], *kho = new double[n], *kgo = new double[n];
	double *tsv = new double[n], t = tei0, ts, g, p, hf = 1e0, r, s;
	if ((!koeq) || (!kho) || (!kgo) || (!tsv)) { cout << snmi << endl; j = getchar(); exit(1); }
	for (k = 0; k < n; k++) tsv[k] = (tgv[k] + thv[k]) / 2e0;
	koeq = koefPribSha(qon, tsv, n, koeq);
	kgo = koefPribSha(tgv, tsv, n, kgo);
	kho = koefPribSha(thv, tsv, n, kho);
	for (k = 0; k < dlma; k++) {
		ts = eftev[k];
		g = 0.0; p = 0.0; for (j = 0; j < n; j++) { g = g + pow(ts, p)*koeq[j]; p = p + hf; } if (g<0.0) g = 0.0; if (qobi) qobi[k] = g;
		p = efktpv[k]; if (fabs(p)>0.0) { g = fabs(g*h / p / 2e0); tgori[k] = eftev[k] + g; tholi[k] = eftev[k] - g; }
		else {
			r = 0.0; s = 0.0; for (j = 0; j < n; j++) { r = r + pow(ts, s)*kgo[j]; s = s + hf; } if (r < 0.0) r = 0.0; if (tgori) tgori[k] = r;
			r = 0.0; s = 0.0; for (j = 0; j < n; j++) { r = r + pow(ts, s)*kho[j]; s = s + hf; } if (r < 0.0) r = 0.0; if (tholi) tholi[k] = r;
		}
	}
	delete[]tsv; delete[]koeq; delete[]kho; delete[]kgo;
}
double *zadrktitom(int zf, int kost, double d, int vyte, double htk, double hvo, int prod, int vy, double tc, int c, int u, double *rta)
{
	int j; double **mu = RaschRTAitom(kost, htk, 0.0, 0.0, 1, vyte, hvo, 0, tc, 1, 0, 0), *tt = mu[0], *te = new double[kost];
	if (!te) { cout << snmi << endl; j = getchar(); exit(1); }
	else for (j = 0; j < kost; j++) te[j] = tt[j];
	if (!u) {
		int k = 0, kst, q = 0, m = 4, b;
		double *prs = new double[m*kost*kost], *pr = new double[m], Er;
		if ((!prs) || (!pr)) { cout << snmi << endl; k = getchar(); exit(1); }
		for (k = 0; k < (m*kost*kost); k++) prs[k] = 0.0;
		q = 0; for (kst = 1; kst <= kost; kst++)
		{
			for (b = 0; b < m; b++) pr[b] = 0.0;
			for (k = 1; k <= kost; k++) {
				pr = izstNitom(k, kst, m, kost);
				for (b = 0; b < m; b++) prs[q + b] = pr[b]; q = q + m;
			}
		}
		if (prod == 1) Er = opredLuchSostitom(prs, te, q, zf, vyte);
		delete[]pr; delete[]te; delete[]tt; delete[]mu; return prs;
	}
	else { delete[]tt; delete[]mu; return te; }
}
double *izstNitom(int izst, int kst, int l, int ocs)
{
	int k; double *o=NULL;
	if (abs(izst - kst) <= N) 
		o = podschchieleSha(izst, kst, ocs, Ritom, Titom);
	else {
		o = new double[l]; if (!o) { cout << snmi << endl; k = getchar(); exit(1); }
		for (k = 0; k < l; k++) o[k] = 0.0;
	}
	return o;
}
double opredLuchSostitom(double *prot, double *Tsr, int le, int zf, int vyte)
{
	int k = 0, j = 0, q = 0, m = 0, l = 3, mma = 30000, mm = 0, mt = 0, vy = 9, idr=2;
	double **mauk, *pt, s = 0.0, r = 0.0, p = 0.0, rez = 0.0, rap = 0.0;
	double *ao = new double[2 * N], *po = NULL, d = 0.0, hb = 5e-1, hb0 = hb, hf = 1e0; for (j = 0; j < ks; j++) d = d + hf;
	double de = hf, ta = 0.0, tb = 0.0, tc = 0.0, fa = 0.0, fb = 0.0, fc = 0.0, ba = bni, bb = bki, bc = (ba + bb) / 2e0, xk = 0.0, dxk = 1e-3;
	double ca = cni, cb = cki, cc = (ca + cb) / 2e0, hc = 1e-1, hc0 = hc, tt = 0.0;
	double mk = 1e2, mkm = mk, c00, b0 = bc, c0 = cc, kdc = 5e0, kdb = 5e0, tminitom, tmaxitom;
	if (!ao) { cout << snmi << endl; k = getchar(); exit(1); }
	for (k = 0; k < (2 * N); k++) ao[k] = 0.0;
	qobitom = opredKTPTKTochSha(qobi, etei, etei[vti], cemi); if (qobitom < 0.0) qobitom = 0.0;
	if ((vti == vtin) && (sctxi <= 2) && (!temrasi)) {
		j = nxti - nnxti; temrasi = new double[j];
		if (!temrasi) { cout << snmi << endl; k = getchar(); exit(1); }
		ta = tholi[vti] - tgori[vti]; ta = ta / y0itom; tb = tgori[vti];
		for (k = 0; k < j; k++) { temrasi[k] = tb + ta*xk; xk = xk + dxk; }
	}
	if ((vti == vtin) && (sctxi <= 2) && (!ktri)) {
		ktri = KoefPoglRosselNac(etei, idr, cemi, smgoi, ssioi, saloi, dkoscim, dkoscit, dkoscil, dkospi, dkoali, kusci, tkusci, dmkoosci, vybitom, 2);
		for (k = 0; k < cemi; k++) cout << "ktr = " << ktri[k] << "\t"; cout << endl;
	}
	l = vti; tminitom = tholi[l]; tmaxitom = tgori[l];
	ta = tminitom; tb = tmaxitom; tc = temrasi[sctxi];
	mauk = RaschRTAitom(ks, hitom, 0.0, 0.0, 1, vti, hvoit, 0, tc, 0, 0, 1);
	k = 0; Raitom = mauk[k]; k++; Taitom = mauk[k]; k++; Aaitom = mauk[k]; k++; Ritom = mauk[k]; k++;
	Titom = mauk[k]; k++; Aitom = mauk[k]; k++; Rtitom = mauk[k]; k++; Ttitom = mauk[k]; k++; Atitom = mauk[k]; delete[]mauk;
	ca = cni; cb = cki; m = 0; q = 0;
	while ((m<mma) && (fabs(cb - ca)>1e-3) && (!q)) {
		cc = (ca + cb) / 2e0;
		pt = FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, cc, vy, qobi, etei, cemi, 2, snmi, dpcti, vti, ektpi, y0itom, ktri, kttki, ks, ktpvoitom, vteitom, dmkvoitom, ecktpitom, hitom, hvoit, Ritom, Titom, Aitom, Rtitom, Ttitom, Atitom, txitom, qxitom, dtxitom, sctxi, tocrasitom, poritom, 2 * N, stchsritom, Tpcti, b0, c0, mk, svfdi, 0);
		k = 3; fc = pt[k]; delete[]pt;
		if (fabs(fc) < 1e0) { c0 = cc; break; }
		pt = FuncRaschIzl(ta, tb, d, ta, m, prot, de, ao, rez, rap, bc, cc, vy, qobi, etei, cemi, 2, snmi, dpcti, vti, ektpi, y0itom, ktri, kttki, ks, ktpvoitom, vteitom, dmkvoitom, ecktpitom, hitom, hvoit, Ritom, Titom, Aitom, Rtitom, Ttitom, Atitom, txitom, qxitom, dtxitom, sctxi, tocrasitom, poritom, 2 * N, stchsritom, Tpcti, b0, c0, mk, svfdi, 0);
		k = 3; fa = pt[k]; delete[]pt;
		if (fabs(fa) < 1e0) { c0 = ca; break; }
		pt = FuncRaschIzl(ta, tb, d, tb, m, prot, de, ao, rez, rap, bc, cc, vy, qobi, etei, cemi, 2, snmi, dpcti, vti, ektpi, y0itom, ktri, kttki, ks, ktpvoitom, vteitom, dmkvoitom, ecktpitom, hitom, hvoit, Ritom, Titom, Aitom, Rtitom, Ttitom, Atitom, txitom, qxitom, dtxitom, sctxi, tocrasitom, poritom, 2 * N, stchsritom, Tpcti, b0, c0, mk, svfdi, 0);
		k = 3; fb = pt[k]; delete[]pt;
		if (fabs(fb) < 1e0) { c0 = cb; break; }
		if ((fa*fc) < 0.0) { cb = cc; q = 0; }
		if ((fb*fc) < 0.0) { ca = cc; q = 0; }
		if (((fb*fc)<0.0) && ((fa*fc)<0.0)) { cb = cc; q = 0; }
		if (((fb*fc)>0.0) && ((fa*fc)>0.0)) q = 1;
		m++;
	}
	c0 = cc; c00 = c0; coi = c00;
	m = 0; ba = bni; bb = bki; l = vti; ta = tholi[l]; tb = tgori[l]; tc = temrasi[sctxi];
	while ((m<mma) && ((bb - ba)>1e-1))
	{
		bc = (ba + bb) / 2e0; b0 = bc; c0 = c00; cc = c0;
		pt = FuncRaschIzl(ta, tb, d, tb, m, prot, de, ao, rez, rap, bc, cc, vy, qobi, etei, cemi, 2, snmi, dpcti, vti, ektpi, y0itom, ktri, kttki, ks, ktpvoitom, vteitom, dmkvoitom, ecktpitom, hitom, hvoit, Ritom, Titom, Aitom, Rtitom, Ttitom, Atitom, txitom, qxitom, dtxitom, sctxi, tocrasitom, poritom, 2 * N, stchsritom, Tpcti, b0, c0, mk, svfdi, 0);
		k = 9; s = pt[k]; k = 7; mk = pt[k]; k = 5; b0 = pt[k]; k = 4; tt = pt[k]; k = 10; r = pt[k]; delete[]pt;
		if ((mk<mkm) && (mk>0.0)) {
			mkm = mk; boi = bc; b0 = bc;
			lxitom[sctxi] = r; ooxitom[sctxi] = s; txitom[sctxi] = tt; qxitom[sctxi] = qobitom;
		}
		m++; if (bc < 1e0) hb = hb0 / kdb; else hb = hb0; bb = bb - hb;
	}
	if (fabs(tnri - tt)>1e-1)
		tnri = tt; else {
			cni = coi - hc0; if (cni <= 0.0) cni = 1e-2;
			cki = coi + 2e0*hc0;
			bni = boi - hb0; if (bni <= 0.0) bni = 1e-2;
			bki = boi + 2e0*hb0;
			tnri = tt;
		}	cout << "qo = " << qobitom << "\ttc = " << tc << endl;
	delete[]ao; return rez;
}
double KorrZnachVozdPrositom(double hps, double ksf, double por, int vy)
{
	int j = 0, k = 1000; double pa = 1e-2, pb = 1e0, *po, pc, ra = fabs(pa - pb), e = tocrasitom, hf = 1e0;
	double fa, fb, fc, ta, tb, tc, tca, tcb, tcc, ka = hps*pa / por, kb = pb*hps / por, kc; //cout << "hps = " << hps << "\tksuf = " << ksf << "\tpor = " << por << endl;
	while ((ra > e) && (j < k)) { //ïîäòÿãèâàåì ïîðèñòîñòü ê çíà÷åíèþ, êîòîðîå çàäàëè èçíà÷àëüíî, âî âðåìÿ ïîäñòðîéêè ÝÊÒÏ
		pc = (pa + pb) / 2e0;
		kc = hps*pc / por;
		po = oprEffDoliTepPerenitom(kc, ksf, pc); tc = po[0]; delete[]po; //ïðè 373 Ê
		tcc = kc*(hf - pc) / pc;
		fc = (hf - tc)*kc / (kc + tcc) - por;
		ka = hps*pa / por;
		po = oprEffDoliTepPerenitom(ka, ksf, pa); ta = po[0]; delete[]po; //îïðåäåëÿåì äîëþ ïëîùàäè ñå÷åíèÿ ïåðåìû÷êè
		tca = ka*(hf - pa) / pa;
		fa = (hf - ta)*ka / (ka + tca) - por;
		kb = hps*pb / por;
		po = oprEffDoliTepPerenitom(kb, ksf, pb); tb = po[0]; delete[]po; //÷åðåç ïåðåìû÷êó òåïëî ðàñïðîñòðàíÿåòñÿ ÷èñòîé òåïëîïðîâîäíîñòüþ
		tcb = kb*(hf - pb) / pb;
		fb = (hf - tb)*kb / (kb + tcb) - por;
		if ((fc*fb > 0.0) && (fa*fc<0.0)) pb = pc; if ((fc*fa>0.0) && (fb*fc < 0.0)) pa = pc;
		j++; ra = fabs(pa - pb);
	}
	dpcti = tc; //cout << "Dol Plo CTP = " << tc << endl;
	if (!vy) return kc; //ñêîððåêòèðîâàííîå çíà÷åíèå òîëùèíû âîçäóøíîé ïðîñëîéêè (ðàçìåð ïîðû), êîãäà ââåëè ïåðåìû÷êó
	else if (vy == 1) return tc;
} //äîëÿ ïëîùàäè, ÷åðåç êîòîðóþ ïðîèñõîäèò ïåðåíîñ òåïëà ÷èñòîé òåïëîïðîâîäíîñòüþ
double *oprEffDoliTepPerenitom(double ko, double d, double por)
{
	int k, j, f, kost = ks, ksu = 2 * ks; double hvozd = ko, htch, hf = 1e0, e = tocrasitom;
	ecktpitom = new double[kost]; if (!ecktpitom) { cout << snmi << endl; k = getchar(); exit(1); }
	double *tepo = new double[ksu], *dol = new double[cemi], r, p, sa, sb, sc, t, fa, fb, fc;
	if ((!tepo) || (!dol)) { cout << snmi << endl; k = getchar(); exit(1); }
	for (j = 0; j < ksu; j++) tepo[j] = 0.0;
	hvozd = ko; htch = hvozd*(hf - por) / por;
	k = 1; for (j = 0; j < cemi; j++) if (qobi[j] < e) k = 0;
	if (!k) opredtemphcitom(ektpi, etei, tgori, tholi, qobi, cemi, dmkoi, y0itom); //cout << "hvoz = " << hvozd << "\tht = " << htch << endl;
	for (j = 0; j < cemi; j++) { //îïðåäåëÿåì ÝÊÒÏ ìíîãîñëîéíîé ñòåíêè
		tepo = opredTempStenShaFragm(tepo, 2 * kost, ktpvoitom, vteitom, etei, kttki, cemi, dmkvoitom, htch, hvozd, qobi[j], etei[j], -hf); //â ñåðåäèíå ñëîÿ
		r = 0.0; for (k = 0; k < ksu; k++) for (f = 0; f<ksu; f++) { p = tepo[k] - tepo[f]; if (p>r) r = p; }
		p = hvozd*(d - hf) + htch*d; r = r / p; t = qobi[j] / r; ecktpitom[j] = t;
	} //for (j=0; j<cemi; j++) cout << "SKTPTK ( " << j << " ) = " << ecktpitom[j] << endl;
	delete[]tepo;
	f = 1000; for (j = 0; j<cemi; j++)
	{
		sa = 0.0; sb = 1e0; k = 0;
		do {
			sc = (sa + sb) / 2e0;
			fa = kttki[j] * sa + ecktpitom[j] * (hf - sa) - ektpi[j]; //ýôôåêòèâíûå ÊÒÏ ìíîãîñëîéíîé ñòåíêè è  ïåðåìû÷êè äîëæíû ñðàâíÿòüñÿ
			fb = kttki[j] * sb + ecktpitom[j] * (hf - sb) - ektpi[j]; //÷òîáû íàéòè îòíîñèòåëüíûå äîëè ïëîùàäåé ñå÷åíèÿ ïåðåíîñà îáùåé ÏÒÏ
			fc = kttki[j] * sc + ecktpitom[j] * (hf - sc) - ektpi[j];
			if ((fc*fb>0.0) && (fa*fc<0.0)) sb = sc;
			if ((fc*fa>0.0) && (fb*fc<0.0)) sa = sc;
			r = fabs(sa - sb); k++;
		} while ((r>e) && (k<f));
		dol[j] = sc;
	}
	return dol;
}
double KorrZnachVozdPrositomKon(double hps, double ksf, double por)
{
	int j = 0, k = 1000; double pa = 1e-3, pb = 1e0, pc, ra = fabs(pa - pb);
	double fa, fb, fc, tca, tcb, tcc, ka = hps*pa / por, kb = pb*hps / por, kc, e = tocrasitom; //cout << "hps = " << hps << "\tksuf = " << ksf << "\tpor = " << por << endl;
	while ((ra>e) && (j < k)) { //ïîäòÿãèâàåì ïîðèñòîñòü ê çíà÷åíèþ, êîòîðîå çàäàëè èçíà÷àëüíî, âî âðåìÿ ïîäñòðîéêè ÝÊÒÏ
		pc = (pa + pb) / 2e0;
		kc = hps*pc / por;
		tcc = kc*(1e0 - pc) / pc;
		fc = (1e0 - dpcti)*kc / (kc + tcc) - por;
		ka = hps*pa / por;
		tca = ka*(1e0 - pa) / pa;
		fa = (1e0 - dpcti)*ka / (ka + tca) - por;
		kb = hps*pb / por;
		tcb = kb*(1e0 - pb) / pb;
		fb = (1e0 - dpcti)*kb / (kb + tcb) - por;
		if ((fc*fb > 0.0) && (fa*fc<0.0)) pb = pc;
		if ((fc*fa>0.0) && (fb*fc < 0.0)) pa = pc;
		j++; ra = fabs(pa - pb);
	} //ñêîððåêòèðîâàííîå çíà÷åíèå òîëùèíû âîçäóøíîé ïðîñëîéêè (ðàçìåð ïîðû), êîãäà ââåëè ïåðåìû÷êó
	return kc;
}
