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
const double tei0=273.15, tnosci=3e2, dtosci=1e2, tnaci=2e2, tki=22.0+tei0;
const double detei=1e2, y0itom = 3e1*1e-3, tocrasitom = 1e-4, templai = 1750.0 + tei0, dkoali = 0.639201;
const int dsi = 60, ks = 10, dmkoosci = 14, dmkvoitom = 28, dmkoi = 3, cemni=11, cemdu=6;
const int N = 13, vtin = 0, nnxti = 0, nxti = 30, vybitom=0, vmik=1, vybves=2;
const char nfi = 'B';
struct derevo {
	int otre; //1 - отражение или 0 - пропускание
	int ste; //номер стенки
	int vis; //видимость: 1 - виден, 0 - нет
	int lev; //номер уровня
	struct derevo *back; //указатель назад
	struct derevo *next; //указатель вперед
};
void RasKTPItom();
double *poisMasKoefItom(int, double *, int);
void initarritom(int);
double *koefoslab(double, double, double, double *, int, double *);
double epsisred(double, double *, double *, int, double *, double *, int, int);
double *opredKTPTverKark(double *, double *, double, double, double, double, int, int, double, double, int, double *, double *, double *, int, int, int);
double **opredTempHolGor(double *, double *, int, int, double, int, double **, int, int, double *, int, double);
double *koefPribSha(double *, double *, int, double *, char *);
void zadrktNac();
void napMasEKTPitomNac(double, double, double);
double RasFracXeffitom63(int);
double *EffectTols(double *, double *, double *, double, double, int);
double **rasPorpoRazitom(int);
double **RaschRTAitom(int, double, double, double, int, int, double, int, double, int, int, int);
double *oprRasTemNach(int, int, double *, double *, int, double *, double, double, double, double, double, int);
double opredKTPTKTochSha(double *, double *, double, int);
double *zadrkt(int, int, double, int, double, double, int, int, double, int, int, double *);
double oprProcSoderpoPoris(double *, double *, double, int);
double opredLuchSost(double *, double *, int, int, int);
double *opredKoefOtr(double *, double *, int, double, int, int, double *, double *, int);
void vyvodfile(double *, int, int, double, char *);
double **vydelPol(int, int, double **, double **, int, int);
double novNapMas(int, int, int, int, int, int, int, int, int);
double KorrZnachVozdPros(double, double, double, int);
double KorrZnachVozdProsKon(double, double, double);
char **napolStrok();
char **napNazvFile(int, int, char *, char);
//-----------
double *dkoscim = NULL, *dkoscit = NULL, *kttki = NULL, *mkoi = NULL;
double *kektpi = NULL;
double poritom=0.0, smgoi=0.0, saloi=0.0;
double ssioi=0.0, *txitom = NULL, *qxitom = NULL, *lxitom = NULL, *kxitom = NULL, *dtxitom = NULL, *tkusci = NULL;
double *kusci = NULL, dkospi = 1e0, hkitom=0.0, hitom=0.0;
double hvoit=0.0, dpcti=1e0;
double *temrasi=NULL, *ooxitom=NULL;
int dkoscil = 6, cemi = cemni, vti = 0, sctxi = 0, vtik=cemi;
//------------
void zadrktNac()
{
int j, jk = nxti, jn = nnxti, k=0, q=0, f = 6, qk=0, vv=2, vm=0, cei=0;
char *snmi=NULL, *sfnoi=NULL, *svfdi=NULL, *s=NULL;
char **ms=napolStrok();
k=0; snmi=ms[k]; k++; sfnoi=ms[k];
ms=napNazvFile(vv, vm, snmi, nfi); if (ms) delete[]ms;
k=0; svfdi=ms[k]; k++; j=4; for (q=k; q<j; q++) { s=ms[q]; if (s) delete[]s; } if (ms) delete[]ms;
	double hf=1e0, nf=0.0, ka=0.0, kb=0.0, hnitom = 0.0, dhk=1e-3, hvh=1e-6, hvko=(13e1)*hvh;
	double hvna=0.0, p=0.0, r=0.0, d=0.0, *atr=NULL, t=0.0, ko=0.0, *po=NULL, **mu=NULL, e=hvh; 
	for (j = 0; j < ks; j++) d = d + hf; 
	//--------------
	double x=0.0, *srra=NULL, *rapo=NULL, *rpr=NULL, srp=0.0, ce=0.0, cet=0.0, *pp=NULL, *legr=NULL, *prgr=NULL;
if (!vybves) mu=rasPorpoRazSha(vybsha, vybstasha);
else if (vybves==1) mu=rasPorpoRazVer(vyfr, vysv);
else if (vybves==2) mu = rasPorpoRazitom(vybitom);
else if (vybves==3) mu = rasPorpoRazkvi(vybkvi);
	k = 0; j=0; rapo = mu[k]; k++; //0
	srra = mu[k]; k++; //1
	prgr = mu[k]; k++; //2
	legr = mu[k]; k++; //3
	pp = mu[k]; srp=pp[j]; k++; if (pp) delete[]pp; //4
	pp = mu[k]; ce=pp[j]; k++; if (pp) delete[]pp; if (mu) delete[]mu; //5
	cet = ce; j = 0; while (cet > e) { cet = cet - hf; j++; } cei = j; delete[]mu;
	k=31; x = 0.0; for (j = 0; j < cei; j++) if (j <= k) x = x + srra[j] * rapo[j]; else break; //for (j=0; j<k; j++) cout << "sr_ra = " << srra[j] << "\tlegr = " << legr[j] << endl;
	//-------------//dkospi = RaschAlphaTvKar(); 
vtik=cemi;
	for (vti = vtin; vti<vtik; vti++) { //пробегаем по температуре
		ko = srp; cout << "Sred razm por = " << ko << "\t";
		r = d; t = ko; ko = KorrZnachVozdPros(ko, r, por, 0); cout << "Korr sred razm por = " << ko << "\t";
		r = d; ko = t; p = KorrZnachVozdPros(ko, r, por, 1); cout << "Dol Plo = " << p << "\t";
		r = d; ko = t; ko = KorrZnachVozdProsKon(ko, r, por); hvna = ko; hvko = ko;
		if (vti>vtin) {
			q = jk - jn; ka = (tholi[vti] - tgori[vti]) / y0itom; kb = tgori[vti]; hkitom = hnitom;
			for (k = 0; k < q; k++) { temrasi[k] = kb + ka*hkitom; hkitom = hkitom + dhk; }
		}
		hvoit = hvna; hitom = hvoit*(hf - poritom) / poritom;
		while (hvoit <= hvko) {
			j = jn; hkitom = hnitom*dhk; //пробегаем по размерам пор
			while ((j < jk) && (hkitom < y0itom)) {
				cout << "hk = " << hkitom << "\tte = " << etei[vti] << endl; sctxi = j; //пробегаем по координате
				atr = zadrkt(j, ks, d, vti, hitom, hvoit, 1, 0, etei[vti], 0, 0, atr);
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
void initarritom(int koel)
{
	double wmg, wsi, wal, ko=1e-2; int k;
	dkoscim = new double[dkoscil]; dkoscit = new double[dkoscil];
	if ((!dkoscim) || (!dkoscit)) { cout << snmi << endl; k = getchar(); exit(1); }
	double tnd = 6e2, dtd = 2e2, tm; k=0; dkoscit[k] = tnd; for (k = 1; k < dkoscil; k++) dkoscit[k] = dkoscit[k - 1] + dtd;
	if (!vybitom) {
		saloi = 26e0*ko; smgoi = 22e0*ko; ssioi = 52e0*ko; //для трехкомпонентной смеси
		wal = 23e0*ko; wmg = 19e0*ko; wsi = 49e0*ko; //для многокомпонентной смеси
		k = 0; dkoscim[k] = 6.07; k++; dkoscim[k] = 5.36; k++; dkoscim[k] = 6.19; k++; 
		dkoscim[k] = 13.48; k++; dkoscim[k] = 19.93; k++; dkoscim[k] = 27.69;
	} //ИТОМ-440
	else if (vybitom == 1) {
		saloi = 29e0*ko; smgoi = 16e0*ko; ssioi = 55e0*ko;
		wal = 26e0*ko; wmg = 15e0*ko; wsi = 5e1*ko; 
		k = 0; dkoscim[k] = 6.41; k++; dkoscim[k] = 5.53; k++; dkoscim[k] = 6.32; k++; 
		dkoscim[k] = 13.56; k++; dkoscim[k] = 19.86; k++; dkoscim[k] = 27.66;
	} //ИТОМ-620
	else if (vybitom == 2) {
		saloi = 33e0*ko; smgoi = 11e0*ko; ssioi = 56e0*ko;
		wal = 3e1*ko; wmg = 1e1*ko; wsi = 52e0*ko; 
		k = 0; dkoscim[k] = 7.28; k++; dkoscim[k] = 6.1; k++; dkoscim[k] = 6.83; k++; 
		dkoscim[k] = 14.02; k++; dkoscim[k] = 19.86; k++; dkoscim[k] = 27.22;
	} //ИТОМ-860
	else if (vybitom == 3) {
		saloi = 35e0*ko; smgoi = 9e0*ko; ssioi = 56e0*ko; 
		wal = 33e0*ko; wmg = 7e0*ko; wsi = 53e0*ko; 
		k = 0; dkoscim[k] = 7.45; k++; dkoscim[k] = 6.24; k++; dkoscim[k] = 6.95; k++; 
		dkoscim[k] = 14.15; k++; dkoscim[k] = 19.89; k++; dkoscim[k] = 27.09;
	} //ИТОМ-1000
	else { cout << "Takoy marki ITOM net!" << endl; k = getchar(); exit(1); }
	for (k = 0; k < dkoscil; k++) { tm = dkoscim[k]*ko; dkoscim[k] = 1e0 - tm; } 
	txitom = new double[koel]; qxitom = new double[koel]; lxitom = new double[koel]; kxitom = new double[koel]; dtxitom = new double[koel]; ooxitom = new double[koel];
	if ((!txitom) || (!qxitom) || (!lxitom) || (!kxitom) || (!dtxitom) || (!ooxitom)) { cout << snmi << endl; j = getchar(); exit(1); }
	for (j = 0; j < koel; j++) { txitom[j] = 0.0; qxitom[j] = 0.0; lxitom[j] = 0.0; kxitom[j] = 0.0; dtxitom[j] = 0.0; ooxitom[j] = 0.0; }
	qobi = new double[cemi]; etei = new double[cemi]; kttki = new double[cemi]; mkoi = new double[cemi]; ektpi = new double[cemi];
	tgori = new double[cemi]; tholi = new double[cemi]; tsreditom = new double[cemi]; stchsritom = new double[cemi];
	if ((!qobi) || (!etei) || (!kttki) || (!mkoi) || (!ektpi) || (!tgori) || (!tholi) || (!stchsritom) || (!tsreditom)) { cout << snmi << endl; j = getchar(); exit(1); }
	for (j = 0; j < cemi; j++) {
		qobi[j] = 0.0; etei[j] = 0.0; kttki[j] = 0.0; mkoi[j] = 0.0;
		ektpi[j] = 0.0; stchsritom[j] = 0.0; tgori[j] = 0.0; tholi[j] = 0.0; tsreditom[j] = 0.0;
	}
	txitom = new double[koel]; qxitom = new double[koel]; lxitom = new double[koel]; kxitom = new double[koel]; dtxitom = new double[koel];
	if ((!txitom) || (!qxitom) || (!lxitom) || (!kxitom) || (!dtxitom)) { cout << snmi << endl; j = getchar(); exit(1); }
	for (j = 0; j < koel; j++) { txitom[j] = 0.0; qxitom[j] = 0.0; lxitom[j] = 0.0; kxitom[j] = 0.0; dtxitom[j] = 0.0; }
	j=0; etei[j] = tnai; for (j = 1; j < cemi; j++) etei[j] = etei[j - 1] + detei;
	tkusci = new double[dmkoosci]; kusci = new double[dmkoosci];
	if ((!tkusci) || (!kusci)) { cout << snmi << endl; j = getchar(); exit(1); }
	k=0; tkusci[k] = tnosci; for (k = 1; k < dmkoosci; k++) tkusci[k] = tkusci[k - 1] + dtosci;
	kusci = koefoslab(wmg, wsi, wal, tkusci, dmkoosci, kusci);
	double s; for (k = 0; k < cemi; k++) { s=epsisreditom(etei[k], tkusci, kusci, dmkoosci, dkoscit, dkoscim, dkoscil, vybitom); stchsritom[k]=s; } //for (k=0; k<cemi; k++) cout << "Step cher ( " << k << " ) = " << stchsritom[k] << "\t"; cout << endl; 
	napMasEKTPitomNac(wmg, wsi, wal);
}
void napMasEKTPitomNac(double wmg, double wsi, double wal)
{
	int vpi=vybitom, k=0, j=0, f=cemdu, vyve=2;
	double *kektpi = new double[dmkoi]; if (!kektpi) { cout << snmi << endl; k=getchar(); exit(1); }
	kektpi=poisMasKoefItom(vpi, kektpi, dmkoi);
	poritom=novNapMas(vybves, vpi, k, k, k, k, k, k, cemi); cout << "por = " << poritom << "\t";	
	double hf = 1e0, g=0.0, s=0.0, t=0.0, *po=NULL, nf = 0.0, tn=0.0, ys0=y0itom;
	double **mu=new double*[f], **muv=new double*[f];
	double *qob=NULL, *thol=NULL, *tgor=NULL, *tsred=NULL;
	if ((!mu) || (!muv)) { cout << snmi << endl; k = getchar(); exit(1); }
	for (k=0; k<f; k++) { po=new double[cemi]; if (!po) { cout << snmi << endl; k=getchar(); exit(1); } muv[k]=po; } 
	for (k = 0; k < cemi; k++) { t=etei[k]; s=0.0; g=0.0; 
	for (j = 0; j<dmkoi; j++) { s = s + kektpi[j] * pow(t, g); g = g + hf; } ektpi[k] = s; } //for (k=0; k<cem; k++) cout << "tem = " << ete[k] << "\tktp = " << ektpi[k] << "\t"; cout << endl;
	muv=opredTempHolGor(ektpi, etei, cemi, k, ys0, k, muv, f, cemi, etei, dmkoi, tki); //cem - длина массива efte
	k=0; j=1; mu=vydelPol(k, cemi, muv, mu, f, j);
	k=0; thol=mu[k]; k++; tgor=mu[k]; k++; qob=mu[k]; k++; 
	ektpi=mu[k]; k++; tsred=mu[k]; etei=tsred; k++; po=mu[k]; 
	k=0; nf=po[k]; t=0.0; k=0; while (t<nf) { t=t+hf; k++; } cemi=k; //cout << "cem = " << cem << "\tnf = " << nf << "\t";	
	for (k=0; k<cemi; k++) cout << "tem = " << etei[k] << "\tktp = " << ektpi[k] << "\tth = " << thol[k] << "\ttg = " << tgor[k] << "\tqo = " << qob[k] << "\t"; cout << endl;
	for (k=0; k<f; k++) { po=muv[k]; if (po) delete[]po; } for (k=0; k<f; k++) { po=mu[k]; if (po) delete[]po; }
	if (muv) delete[]muv; if (mu) delete[]mu; if (kektpi) delete[]kektpi; k=getchar(); //for (k = 0; k < cemi; k++) cout << "t_g = " << tgis[k] << "\tt_h = " << tkhis[k] << "\tqis = " << qis[k] << "\ttsis = " << tsis[k] << "\tete = " << etei[k] << "\tektpi = " << ektpi[k] << endl;
	if (stchsritom) delete[]stchsritom; stchsritom = new double[cemi]; if (!stchsritom) { cout << snmi << endl; k = getchar(); exit(1); }
	for (k = 0; k < cemi; k++) { g = epsisred(etei[k], tkusci, kusci, dmkoosci, dkoscit, dkoscim, dkoscil, vybitom); stchsritom[k] = g; } 
	for (k = 0; k < cemi; k++) cout << "tem = " << etei[k] << "\tst_ch = " << stchsritom[k] << "\t"; cout << endl;
	kttki = opredKTPTverKarkitom(etei, ektpi, poritom, wsi, wal, wmg, vybitom, dmkoi, tnosci, dtosci, dmkoosci, kusci, tkusci, stchsritom, cemi, vmik, vpi); 
}
void RasKTPItom()
{
	int k, f = 4;
	double *sv = new double[f], kk=1e-2;
	k=0; sv[k] = 6e1; k++; sv[k] = 5e1; k++; sv[k] = 25.0; k++; sv[k] = 2e1; 
	for (k = 0; k < f; k++) sv[k] = sv[k] * kk; //содержание вермикулита
	delete[]sv;
}
double *poisMasKoefItom(int no, double *kti, int n)
{
	int f=n, k=0; double t1=2e2, t2=38e1, kk=1e-2, te200=t1+tei0, te380=t2+tei0, *ktpit=NULL;
	if (!no) {
		double *kitom440 = new double[f]; ktpit = new double[f];
		if ((!kitom440) || (!ktpit)) { cout << snmi << endl; getchar(); exit(1); }
		for (k=0; k<f; k++) { ktpit[k]=0.0; kitom440[k]=0.0; }
		k=0; ktpit[k] = 9e0*kk; k++; ktpit[k] = 12e0*kk;
		kitom440[1] = (ktpit[1] - ktpit[0]) / (te380 - te200); 
		kitom440[0] = ktpit[0] - kitom440[1] * te200;
		for (k=0; k<f; k++) kti[k] = kitom440[k]; if (kitom440) delete[]kitom440;
	}
	else if (no == 1) {
		double *kitom620 = new double[f]; ktpit = new double[f];
		if ((!kitom620) || (!ktpit)) { cout << snmi << endl; getchar(); exit(1); } 
		for (k=0; k<f; k++) { ktpit[k]=0.0; kitom620[k]=0.0; }
		//ktpit[0] = 12.0*1e-2; ktpit[1] = 139.0*1e-3; //из Диссертации
		k=0; ktpit[k] = 18.0*kk; k++; ktpit[k] = 19.0*kk; //Данные 2017 года
		kitom620[1] = (ktpit[1] - ktpit[0]) / (te380 - te200); 
		kitom620[0] = ktpit[1] - kitom620[1] * te380;
		for (k=0; k<f; k++) kti[k] = kitom620[k]; if (kitom620) delete[]kitom620;
	}
	else if (no == 2) {
		double *kitom860 = new double[f]; ktpit = new double[f];
		if ((!kitom860) || (!ktpit)) { cout << snmi << endl; getchar(); exit(1); } 
		for (k=0; k<f; k++) { ktpit[k]=0.0; kitom860[k]=0.0; }
		//ktpit[0] = 18.3*kk; ktpit[1] = 19.4*kk; //из Диссертации
		ktpit[0] = 26.0*kk; ktpit[1] = 37.0*kk; //Данные 2017 года
		kitom860[1] = (ktpit[1] - ktpit[0]) / (te380 - te200); kitom860[0] = ktpit[1] - kitom860[1] * te380;
		for (k=0; k<f; k++) kti[k] = kitom860[k]; if (kitom860) delete[]kitom860;
	}
	else if (no == 3) {
		double *kitom1000 = new double[f]; ktpit = new double[f];
		if ((!kitom1000) || (!ktpit)) { cout << snmi << endl; getchar(); exit(1); } 
		for (k=0; k<f; k++) { ktpit[k]=0.0; kitom1000[k]=0.0; }
		//ktpit[0] = 23.0*kk; ktpit[1] = 25.0*kk; //из Диссертации
		ktpit[0] = 42.0*kk; ktpit[1] = 52.0*kk; //Данные 2017 года
		kitom1000[1] = (ktpit[1] - ktpit[0]) / (te380 - te200); kitom1000[0] = ktpit[1] - kitom1000[1] * te380;
		for (k=0; k<f; k++) kti[k] = kitom1000[k]; if (kitom1000) delete[]kitom1000;
	}
	else { cout << "Net takoy marki ITOM!" << endl; k = getchar(); exit(1); }
	if (ktpit) delete[]ktpit;
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
	if (mkbr) delete[]mkbr; if (mv) delete[]mv; if (tol) delete[]tol; if (xv) delete[]xv; return xsr;
}