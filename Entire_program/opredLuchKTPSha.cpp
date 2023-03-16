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
const int ceektp = 4, ds = 50, N = 13, ks = 10, dmkooscs = 14, nxts=65;
const int nnxts=0, vtn=0, vrsh=3, cemns=11, cemdu=6, vybves=0;
const double tnoscs = 3e2, tnacs=2e2, detes = 1e2, nnxtfs = 0.0, nxtfs = 65e0, y0sha = nxtfs*1e-3;
const double dtoscs = 1e2, pksvs = 0.0, tocrassha = 1e-7, epsi = 1e-15;
const double te0 = 273.15, dkoals = 0.647528, tks=22.0+te0, templas=175.0*1e1+te0; //dkoals - получение КП из КО
const char nf = 'F';
struct derevo {
	int otre; //1 - отражение или 0 - пропускание
	int ste; //номер стенки
	int vis; //видимость: 1 - виден, 0 - нет
	int lev; //номер уровня
	struct derevo *back; //указатель назад
	struct derevo *next; //указатель вперед
};
double *lxsha = NULL, *qxsha = NULL, *tkuscs = NULL, *kuscs = NULL, *qob = NULL, *tems = NULL, *etesha = NULL, *dtxsha = NULL;
double *dkosct = NULL, *kttk = NULL, *stchsrsha = NULL, *ktr = NULL, *ecktpsha = NULL, *Tpct = NULL, *dkoscm = NULL; //dkoscm - отклонение СЧ (или ПС) от литературных данных
double hsha, hvsha, hksha = 0.0, qobsha = 0.0, porsha = 16.5*1e-2, tos = 1750.0 + te0, tnas = tnacs+te0, dkosps = 0.665; //dkosps - дополнительный коэффициент для КО с учетом пористой структуры;
double *mko=NULL, *kektp=NULL, *ektpsha=NULL; //dpcts - доля площади перемычки, через которую переносится чистой теплопроводностью ТК
double *alphaSha=NULL, *txsha=NULL, *kxsha=NULL, tnr=0.0, dpcts=0.314;
double cn=1e-2, ck=4e0, *temrass=NULL, *ooxsha=NULL, sfeosha=1.24e-2, smgosha=0.29e-2, salosha=41.9e-2, ssiosha=54e-2, *dpctsm=NULL;
int vtsh = 0, sctx = 0, dkoscl = 6, vystsha = 1, vybsha = 3, cem = cemns;
//----------------------------------------
double *podschchiele(int, int, int, double *, double *);
double *izstN(int, int, int, int);
double rasotprpovsSha(int *, int *, int, int, int, int, double *, double *);
double otrprovsSha(int, int, double *, double *);
double opredLuchSostSha(double *, double *, int, int);
double *koefPribSha(double *, double *, int, double *);
double *zadrkt(int, int, double, int, double, double, int, double, int, int, double *);
double **chaRTASha(int, double *, double *, double *, double *, double *, double *, int);
double **izmRTASha(double *, int, double *, double *, double *, double *, double *, double *, int, int);
double *opredKoefOtr(double *, double *, int, double, int, int, double *, double *, int);
void zadrktNac();
double **NapMasVozdSha();
double *opredKTPTverKarkSha(double *, double *, double, double, double, double, double, double, int, double *, double *, double *, int, int, int);
double opredKTPTKTochSha(double *, double *, double, int);
double *usrednen(double *, double *, int, int);
void vyvodfile(double *, int, int, double, char *);
void napstrsha();
void PrintMatr(double **, int);
void osvopam();
void iniarrele(int, double, double, double);
double **opredTempLPStenSha(double *, double *, double *, double *, double, int, double *, double *, double *, double *, int, int, double, double, double, int, char *);
double **RaschVneshIzluchSha(double *, double *, double *, double *, double *, int, double, double, char *);
double **RasLuchPloTepPot(int, double *, double *, double *, double *, double *, double *, int, double *, double *, double *, char *);
double **RaschSobLuchPlotTepPot(int, double *, double *, double *, double *, double *, double *, double, double, double *, int, char *);
double opredUrovPodderM03(double);
double *opredTempStenShaFragm(double *, int, double *, double *, double *, double *, int, int, double, double, double, double, double);
double *KoefPoglRosselNac(double *, int, int, double, double, double, double *, double *, int, double, double, double *, double *, int, int, int);
double LuchKTPChudnovsky(double *, double, int, double);
double epsisredsha(double, double *, double *, int, double *, double *, int);
double **RaschRTASha(int, double, double, double, double, double, int, int, int, int);
double *oprRasTemNach(int, int, double *, double *, int, double *, double, double, double, double, double, int);
double **vybFunRasPorpoRazSha(double, int, int);
double *oprEffDoliTepPeren(double, double, double);
double *KorrZnachVozdPros(double, double, double, int);
double *koefoslab(double, double, double, double *, int, double *);
double vychNevyaz(double **, double **, double *, double *, int);
double *opredPolTempTvKarShaFragm(double *, int, double *, double *, int, double, double, double, double, double);
double *EffectTols(double *, double *, double *, double, double, int);
double SeryeStenkiRasIzl(double, double, double, double *, double *, double *, double *, double *, int *, int);
double BolTochRasAlpha(int, int, double, double, double, double *, double *, double *, int, int, char *, double *, double *, int, double *, double *, int, double *, double *, double *, double *, double *, double *);
double RasFracXeffSha60(int);
double RaschAlphaTvKarSha();
double RasIzlSerStenNac(double *, double *, double *, double *, double *, double *, double *, double, double, double, double, int, double, double, double, double, int, double *, double *, int, int, double *, double *, double *, char *, int, double *, double *, int);
double KorrZnachVozdProsSham(double, double, double);
double *rasMetNewtRafs(double *, double **, double **, double *, double, double, int);
double *FuncRaschIzl(double, double, double, double, int, double *, double, double *, double, double, double, double, int, double *, double *, int, int, char *, double, int, double *, double, double *, double *, int, double *, double *, int, double *, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, int, double, double, int, double *, double *, double, double, double, char *, int);
double **PoiskZavVelTemVer(int, double **, int, double *, double *, int);
void zapisvfile(double *, int, char *);
double *poisMasKoefsha(double *, int);
void napMasEKTPshaNac();
double **vydelPol(int, int, double **, double **, int, int);
double novNapMas(int, int, int, int, int, int, int, int, int);
double **opredTempHolGor(double *, double *, int, int, double, int, double **, int, int, double *, int, double);
double *reshnewtrafs(double *, double *, double *, double *, int, double *, double *, int, double *, double *, int, double *, double *, int, double, double *, double *, int, double, double, double, double *, int, double, double, double, char *);
//----------------------------------------
double *poisMasKoefsha(double *kektp, int n)
{
	double ko = 1e-2, p20 = 2e1*ko, p24 = 24.0*ko, p30 = 3e1*ko, p33 = 33.0*ko, p16 = 16.0*ko, p10 = 1e1*ko; int k=0;
	double s28=28e-2, s38=38e-2, s45=45e-2;
	if (!vrsh) { sfeosha = 1.24*ko; smgosha = 0.29*ko, salosha = 41.9*ko, ssiosha = 54.0*ko; porsha=21.8*ko; } //Роучка
	else if (vrsh==1) { sfeosha = 1.21*ko; smgosha = 0.29*ko; salosha = 42.6*ko; ssiosha = 1e0-salosha-sfeosha-smgosha; porsha=11.0144*ko; } //ШПД-41
	else if (vrsh==2) { sfeosha = 1.64*ko; smgosha = 0.36*ko; salosha = 35.9*ko; ssiosha = 59.1e-2;  porsha=25.2*ko; } //ШБ-1 2-1
	else if (vrsh==3) { sfeosha = 1.66*ko; smgosha = 0.4*ko; salosha = 37.3*ko; ssiosha = 57.4*ko;  porsha=26.5*ko; } //ШВ-1 1-1
	else if (vrsh==4) { sfeosha = 1.24*ko; smgosha = 0.29*ko; salosha = 41.9*ko; ssiosha = 54*ko;  porsha=11.5*ko; } //ШПД
	else if (vrsh==5) { sfeosha = 1.54*ko; smgosha = 0.3*ko; salosha = 38.6*ko; ssiosha = 56.5*ko;  porsha=16.5*ko; } //ШКУ-32 3-1
	for (k=0; k<n; k++) kektp[k]=0.0;
	if ((salosha >= s28) && (salosha <= s38)) {
		if ((porsha >= p20) && (porsha < p24)) { vybsha = 0; vystsha = 0; kektp[3] = -0.435e-9; kektp[2] = 0.685e-6; kektp[1] = 0.134e-3; kektp[0] = 0.725; }
		else if ((porsha >= p24) && (porsha < p30)) { vybsha = 0; vystsha = 1; kektp[3] = -0.867e-9; kektp[2] = 1.77e-6; kektp[1] = -0.523e-3; kektp[0] = 0.806; }  //задание коэффициентов - шамот средней пористости
		else if ((porsha >= p16) && (porsha < p20)) { vybsha = 1; kektp[3] = -0.397e-9; kektp[2] = 0.71e-6; kektp[1] = 0.011e-3; kektp[0] = 0.851; } //уплотненный шамот
		else if ((porsha >= p30) && (porsha <= p33)) { vybsha = 2; kektp[3] = -0.377e-9; kektp[2] = 0.918e-6; kektp[1] = -0.338e-3; kektp[0] = 0.77; }  //низкоплотный шамот
		else if ((porsha >= p10) && (porsha<p16)) { vybsha = 3; kektp[3] = 0.0; kektp[2] = -0.607e-6; kektp[1] = 1.14e-3; kektp[0] = 0.641; }
	} //повышенной плотности шамот
	if ((salosha>s38) && (salosha <= s45)) {
		if ((porsha >= p20) && (porsha < p24)) { vybsha = 0; vystsha = 0; kektp[3] = -0.124e-9; kektp[2] = 0.215e-6; kektp[1] = 0.125e-3; kektp[0] = 1.01; }
		else if ((porsha >= p24) && (porsha < p30)) { vybsha = 0; vystsha = 1; kektp[3] = -0.333e-9; kektp[2] = 0.805e-6; kektp[1] = -0.289e-3; kektp[0] = 0.903; } //задание коэффициентов - шамот средней пористости
		else if ((porsha >= p16) && (porsha < p20)) { vybsha = 1; kektp[3] = 0.0; kektp[2] = -0.154e-6; kektp[1] = 0.369e-3; kektp[0] = 1.03; } //уплотненный шамот
		else if ((porsha >= p30) && (porsha < p33)) { vybsha = 2; kektp[3] = -0.377e-9; kektp[2] = 0.918e-6; kektp[1] = -0.338e-3; kektp[0] = 0.77; }  //низкоплотный шамот
		else if ((porsha >= p10) && (porsha < p16)) { vybsha = 3; kektp[3] = 0.0; kektp[2] = -0.141e-6; kektp[1] = 0.437e-3; kektp[0] = 1.32; }
	} //повышенной плотности шамот
	return kektp;
}
void iniarrele(int n, double wmg, double wsi, double wal)
{
	int j, k; tems = new  double[cem]; kttk = new double[cem]; stchsrsha = new double[cem];
	mko = new double[cem]; ektpsha = new double[cem]; etesha = new double[cem]; ecktpsha = new double[cem]; 
	if ((!tems) || (!kttk) || (!mko) || (!ektpsha) || (!etesha) || (!stchsrsha) || (!ecktpsha))
	{ cout << snms << endl; k = getchar(); exit(1); }
	for (j = 0; j < cem; j++) { qob[j] = 0.0; tems[j] = 0.0; kttk[j] = 0.0; 
	stchsrsha[j] = 0.0; mko[j] = 0.0; ektpsha[j] = 0.0; etesha[j] = 0.0; ecktpsha[j] = 0.0; }
	j=0; etesha[j]=tnas; for (j = 1; j < cem; j++) etesha[j] = etesha[j - 1] + detes;
	dkoscm = new double[dkoscl]; dkosct = new double[dkoscl];
	if ((!dkoscm) || (!dkosct)) { cout << snms << endl; k = getchar(); exit(1); }
	double tnd=6e2, dtd=2e2, s=0.0, tm=0.0; k=0; dkosct[k] = tnd; for (k = 1; k < dkoscl; k++) dkosct[k] = dkosct[k - 1] + dtd;
	k = 0; dkoscm[k] = 8.15; k++; dkoscm[k] = 6.87; k++; dkoscm[k] = 7.55; k++; dkoscm[k] = 14.84; k++; dkoscm[k] = 20.09; k++; dkoscm[k] = 26.3;
	for (k = 0; k < dkoscl; k++) { tm = dkoscm[k] / 1e2; dkoscm[k] = 1e0 - tm; }
	alphaSha = new double[ks]; Tpct = new double[ks];
	if ((!alphaSha) || (!Tpct)) { cout << snms << endl; k = getchar(); exit(1); }
	for (j = 0; j < ks; j++) { alphaSha[j] = 0.0; Tpct[j] = 0.0; }
	txsha = new double[n]; kxsha = new double[n]; lxsha = new double[n]; qxsha = new double[n]; dtxsha = new double[n]; ooxsha = new double[n];
	if ((!txsha) || (!kxsha) || (!lxsha) || (!qxsha) || (!dtxsha) || (!ooxsha)) { cout << snms << endl; k = getchar(); exit(1); }
	for (j = 0; j < n; j++) { txsha[j] = 0.0; kxsha[j] = 0.0; lxsha[j] = 0.0; qxsha[j] = 0.0; dtxsha[j] = 0.0; ooxsha[j] = 0.0; }
	tkuscs = new double[dmkooscs]; kuscs = new double[dmkooscs];
	if ((!tkuscs) || (!kuscs)) { cout << snms << endl; k = getchar(); exit(1); }
	k=0; tkuscs[k] = tnoscs; for (k = 1; k < dmkooscs; k++) tkuscs[k] = tkuscs[k - 1] + detes;
	kuscs = koefoslab(wmg, wsi, wal, tkuscs, dmkooscs, kuscs);
	for (k = 0; k < cem; k++) { s = epsisredsha(etesha[k], tkuscs, kuscs, dmkooscs, dkosct, dkoscm, dkoscl); stchsrsha[k] = s; }
	kttk = opredKTPTverKarkSha(etesha, ektpsha, porsha, wsi, wal, wmg, tnoscs, dtoscs, dmkooscs, kuscs, tkuscs, stchsrsha, cem, vybsha, vystsha);
	NapMasVozd(); //for (k = 0; k < cem; k++) cout << "te = " << etesha[k] << "\tektp = " << ektpsha[k] << endl;
}
void napMasEKTPshaNac()
{
	int vps=vybsha, k=0, j=0, f=cemdu, cems=cem, dmkos=ceektp, vyve=0;
char *snms=NULL, *sfnos=NULL, *svfds=NULL, *s=NULL;
char **ms=napolStrok();
k=0; snms=ms[k]; k++; sfnos=ms[k];
ms=napNazvFile(vyve, vm, snmi, nfs); if (ms) delete[]ms;
k=0; svfds=ms[k]; k++; j=4; for (q=k; q<j; q++) { s=ms[q]; if (s) delete[]s; } if (ms) delete[]ms;
	double *kektps=new double[dmkos], ys0=y0sha; if (!kektps) { cout << snms << endl; k=getchar(); exit(1); }
	kektps=poisMasKoefsha(kektps, dmkos);
	porsha=novNapMas(vybves, vps, k, k, k, k, k, k, cem);
	double hf=1e0, g=0.0, s=0.0, t=0.0, *po=NULL, nf = 0.0, tn=0.0;
	double **mu=new double*[f], **muv=new double*[f];
	if ((!mu) || (!muv)) { cout << snms << endl; k = getchar(); exit(1); }
	for (k=0; k<f; k++) { po=new double[cems]; if (!po) { cout << snms << endl; k=getchar(); exit(1); } muv[k]=po; } 
	for (k=0; k<cems; k++) { t=etesha[k]-te0; s=0.0; g=0.0; 
	for (j=0; j<dmkos; j++) { s=s+kektps[j]*pow(t, g); g=g+hf; } ektpsha[k]=s; } //for (k=0; k<cem; k++) cout << "tem = " << ete[k] << "\tktp = " << ektpi[k] << "\t"; cout << endl;
	muv=opredTempHolGor(ektpsha, etesha, cems, k, ys0, k, muv, f, cems, etesha, dmkos, tks); //cem - длина массива efte
	k=0; j=1; mu=vydelPol(k, cems, muv, mu, f, j);
	k=0; tholsha=mu[k]; k++; tgorsha=mu[k]; k++; qob=mu[k]; k++; 
	ektpsha=mu[k]; k++; etesha=mu[k]; k++; po=mu[k]; 
	k=0; nf=po[k]; t=0.0; k=0; while (t<nf) { t=t+hf; k++; } cems=k; //cout << "cem = " << cem << "\tnf = " << nf << "\t";	
	for (k=0; k<cems; k++) { cout << "tem = " << etesha[k] << "\tktp = " << ektpsha[k] << "\tth = " << tholsha[k] << "\ttg = " << tgorsha[k] << "\tqo = " << qob[k] << "\t"; cout << endl; }
	for (k=0; k<f; k++) { po=muv[k]; delete[]po; } for (k=0; k<f; k++) { po=mu[k]; delete[]po; }
	if (muv) delete[]muv; if (mu) delete[]mu; if (kektps) delete[]kektps; k=getchar();
}
void osvopam()
{
	delete[]Tsh; delete[]Ash; delete[]Rsh; delete[]Tash; delete[]Aash; delete[]Rash; delete[]Ttsh; delete[]Atsh;
	delete[]Rtsh; delete[]alphaSha; delete[]txsha; delete[]kxsha; delete[]lxsha; delete[]qxsha; delete[]dtxsha; 
	delete[]qob; delete[]tems; delete[]stchsrsha; delete[]kttk; delete[]tholsha; delete[]tgorsha; delete[]mko; delete[]kektp; 
	delete[]etesha; delete[]ektpsha; delete[]ktpvosha; delete[]vtesha; delete[]ecktpsha; delete[]Tpct; delete[]dkoscm; 
	delete[]dkosct; delete[]kuscs; delete[]tkuscs; delete[]ktr; delete[]sns; delete[]sss, delete[]skpts; delete[]skptsk;
	delete[]svsk; delete[]svs; delete[]snms; delete[]sfnos; delete[]temrass; delete[]ooxsha;
}
void zadrktNac()
{
	int j, jk = nxts, jn = nnxts, k, f = 6, q; napstrsha(); iniarrele(jk-jn, smgosha, ssiosha, salosha);
	double dhk = y0sha/fabs(nxtfs-nnxtfs), hnsha = nnxtfs*dhk, hvko = 130e-6;
	double hvh = 1e-6, hvna = 1e-6, p, r, s; dkosps = RaschAlphaTvKarSha(); //dkosps - дополнительный коэффициент для КО с учетом пористой структуры
	double d = 0.0, *mw = NULL, *po, ko, **mauk=NULL, hf = 1e0, ka, kb, *ras;
	double *srra, *legr, *prgr, marp, srp; for (j = 0; j < ks; j++) d = d + hf; //cout << "Chislo stenok = " << d << endl;
	for (vtsh = vtn; vtsh < cem; vtsh++) { //пробегаем по температуре
		if (vtsh == vtn) {
			mauk=vybFunRasPorpoRazSha(porsha, vrsh, vystsha);
			k=0; ras=mauk[k]; k++; srra=mauk[k]; k++; prgr=mauk[k]; k++; legr=mauk[k]; k++; 
			po=mauk[k]; srp=po[0]; k++; if (po) delete[]po; po=mauk[k]; marp=po[0]; ko=srp; 
			if (po) delete[]po; delete[]mauk; cout << "Sred razm por = " << ko << "\t"; 
			if (ras) delete[]ras; if (srra) delete[]srra; if (legr) delete[]legr; if (prgr) delete[]prgr; 
			r = d; s = ko; po = KorrZnachVozdPros(ko, r, porsha, 0); ko=po[0]; cout << "Korr sred razm por = " << ko << "\t";
			r = d; ko = s; po = KorrZnachVozdPros(ko, r, porsha, 1); p=po[0]; cout << "Dol Plo = " << p << "\t";
			r = d; ko = s; ko = KorrZnachVozdProsSham(ko, r, porsha);
			bk = 2e0 / pow(hf - dpcts, 2e0); hvna = ko; hvko = ko;  }
		else {
			q = jk - jn; ka = (tholsha[vtsh]-tgorsha[vtsh]) / y0sha; kb=tgorsha[vtsh]; hksha = hnsha;
			for (k = 0; k < q; k++) { temrass[k] = kb + ka*hksha; hksha = hksha + dhk; } }
		hvsha = hvna; hsha = hvsha*(hf - porsha) / porsha;
		while (hvsha <= hvko) {
			j = jn; hksha = hnsha; //пробегаем по размерам пор
			while (j < jk) {
				cout << "hk = " << hksha << "\tte = " << etesha[vtsh] << endl; sctx = j-jn; p = d; //пробегаем по координате
				mw = zadrktSha(j, ks, p, vtsh, hsha, hvsha, 1, etesha[vtsh], 0, 0, mw); kxsha[j-jn] = hksha;
				hksha = hksha + dhk; delete[]mw; j++;
			} vyvodfile(lxsha, jk - jn, 0, hvsha, svs); 
			vyvodfile(ooxsha, jk - jn, 2, 0.0, svs); 
			for (j = 0; j < (jk-jn-1); j++) { r = fabs((txsha[j] - txsha[j + 1]) / (kxsha[j + 1] - kxsha[j])); 
			p = fabs(qxsha[j + 1] + qxsha[j]) / 2e0; 
			lxsha[j] = p / r; //cout << "kx = " << kxsha[j] << "\ttx = " << txsha[j] << "\tqx = " << qxsha[j] << "\tlx = " << lxsha[j] << endl; 
			} vyvodfile(lxsha, jk-jn-1, 0, porsha, svs); 
			vyvodfile(qxsha, jk - jn, 2, 0.0, svs); 
			vyvodfile(txsha, jk - jn, 2, 0.0, svs); 
			hvsha = hvsha + hvh;
		}
	}
	osvopam();
}
double *podschchiele(int no, int kste, int ocs, double *Ref, double *Tra)
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
double *koefPribSha(double *ktp, double *te, int le, double *ko, char *snm)
{
	int k=0, kem=3; double **A=new double*[kem], *AA=NULL, *b=new double[kem], yx2=0.0, yx=0.0, p=0.0, hf=1e0;
	if ((!A) || (!b) || (!ko)) { cout << snm << endl; k=getchar(); exit(1); }
	double x4 = 0.0, x3 = 0.0, x2 = 0.0, x = 0.0, y = 0.0, de = 0.0, de1 = 0.0, de2 = 0.0, de3 = 0.0;
	if (A) { for (k = 0; k < kem; k++) { AA = new double[kem]; 
	if (AA) A[k] = AA; else { cout << snm << endl; k = getchar(); exit(1); } } }
	else { cout << snm << endl; k = getchar(); exit(1); }
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
	if (b) delete[]b; for (k = 0; k < kem; k++) { AA = A[k]; if (AA) delete[]AA; } 
	return ko;
}
double **NapMasVozd(char *snm)
{
	int k = 0, dlma=28; double *ktpvoz, *tvoz;
	ktpvoz = new double[dlma]; tvoz = new double[dlma]; 
	if ((!ktpvoz) || (!tvoz)) { cout << snm << endl; k = getchar(); exit(1); }
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
	k=2; double **unau=new double*[k];
	k=0; unau[k]=ktpvoz; k++; unau[k]=tvoz;
	return unau;
}
double *usrednen(double *ao, double *usr, int sao, int kost)
{
	double s = 0.0, k = 0.0, minz = tocrassha, hf = 1e0; int j;
	for (j = 0; j<kost; j++) if (fabs(usr[j])>minz) { s = s + usr[j]; k = k + hf; }
	s = s / k; ao[sao] = s; //средняя величина
	return ao;
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
double **RaschVneshIzluch(double *Tra, double *Ref, double *prot, double *hlr1, double *hrl1, int kost, double hlr0, double hrl0, char *snm)
{
	int k = 3, j, m, q = 2; double *khrl = new double[kost], *khlr = new double[kost], *ri = new double[q], **mu = new double*[k], e=1e-15;
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
			if (fabs(wk)>e) {
				tta = tta + tt*khlr[k] * w0 / wk; //доля энергии, отраженная от k-ой стенки, и прошедшая через нее, пришла слева к j-ой стенке
				ttb = ttb + tt*khlr[k] * w1 / wk;
			} //для учета при подсчете лучистой энергии, распространяющейся справа налево
			q = q + 4;
		} q = m + 4 * kost - 1; m = q; //считаем, что стенки, до которых дошло излучение, как бы "сами" излучают, из-за того, что на них попало внешнее излучение
		for (k = (kost - 1); k >= 0; k--) {
			w0 = prot[q - 3]; w1 = prot[q - 2]; 
			w2 = prot[q - 1]; w3 = prot[q]; 
			wk = w2 + w3; q = q - 4; //для внешнего излучения, распространяющегося справа налево
			tt = Ref[k] * w3 + Tra[k] * w2; //деление на два потока
			if (fabs(wk) > e) {
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
	int k=0; double lamt, lamvo, lam, ht, tt; Tna[k] = tenac;
	for (k = 1; k < n; k++)  {
		tt = Tna[k - 1];
		lamt = opredKTPTKTochSha(ktptk, ete, tt, ce); //ищем КТП твердого скелета при данной температуре
		lamvo = opredKTPTKTochSha(ktpvo, vte, tt, dmv); //ищем КТП воздуха при данной температуре
		if (k % 2) { ht = htk; lam = lamt; }
		else { ht = hnf; lam = lamvo; }
		Tna[k] = tt + zn*ptpob*ht / lam;  //рассчитываем температуру на границах стенок 
	} return Tna;
}
double *oprEffDoliTepPeren(double ko, double d, double por)
{
	int k=0, j=0, f=cemdu, kost=ks, ksu=2*ks, q=0; double hvozd=0.0, htch=0.0, hf=1e0, e=tocrassha;
	double *tepo=new double[ksu], *dol=new double[cem], r=0.0, p=0.0, sa=0.0, sb=0.0, sc=0.0;
	double t=0.0, fa=0.0, fb=0.0, fc=0.0, **mu=new double*[f];
	if ((!tepo) || (!dol)) { cout << snms << endl; k = getchar(); exit(1); }
	for (j = 0; j < ksu; j++) tepo[j] = 0.0;
	hvozd = ko; htch = hvozd*(hf - por) / por;
	for (j = 0; j < cem; j++) { //определяем ЭКТП многослойной стенки
		tepo = opredTempStenShaFragm(tepo, ksu, ktpvosha, vtesha, etesha, kttk, cem, dmkv, htch, hvozd, qob[j], etesha[j], -hf); //в середине слоя
		r = 0.0; for (k=0; k<ksu; k++) for (q=0; q<ksu; q++) { p=tepo[k]-tepo[q]; if (p>r) r = p; }
		p = hvozd*(d - hf) + htch*d; r = r / p; t = qob[j] / r; ecktpsha[j] = t;
	} //for (j=0; j<cem; j++) cout << "SKTPTK ( " << j << " ) = " << ecktpsha[j] << endl;
	f = 1000; for (j = 0; j<cem; j++)
	{
		sa = 0.0; sb = 1e0; k = 0;
		do {
			sc = (sa + sb) / 2e0;
			fa = kttk[j] * sa + ecktpsha[j] * (hf - sa) - ektpsha[j]; //эффективные КТП многослойной стенки и перемычки должны сравняться
			fb = kttk[j] * sb + ecktpsha[j] * (hf - sb) - ektpsha[j]; //чтобы найти относительные доли площадей сечения переноса общей ПТП
			fc = kttk[j] * sc + ecktpsha[j] * (hf - sc) - ektpsha[j];
			if ((fc*fb>0.0) && (fa*fc<0.0)) sb = sc;
			if ((fc*fa>0.0) && (fb*fc<0.0)) sa = sc;
			r = fabs(sa - sb); k++;
		} while ((r>tocrassha) && (k<f));
		dol[j] = sc;
	}
	if (mu) delete[]mu; if (tepo) delete[]tepo;
	return dol;
}
double *KorrZnachVozdPros(double hps, double ksf, double por, int vy)
{
	int j = 0, k = 1000, q=0; double pa = 1e-2, pb = 1e0, *po, pc, ra = fabs(pa - pb), *pkc=new double[1];
	double fa, fb, fc, ta, tb, tc, tca, tcb, tcc, ka = hps*pa / por, kb = pb*hps / por, kc; //cout << "hps = " << hps << "\tksuf = " << ksf << "\tpor = " << por << endl;
	dpctsm=new double[cem]; if (!dpctsm) { cout << snms << endl; q=getchar(); exit(1); }
	for (q=0; q<cem; q++) { j=0; pa = 1e-2; pb = 1e0; ra = fabs(pa - pb); ka = hps*pa / por; kb = pb*hps / por;
	while ((ra>tocrassha) && (j < k)) { //подтягиваем пористость к значению, которое задали изначально, во время подстройки ЭКТП
		pc = (pa + pb) / 2e0; kc = hps*pc / por;
		po = oprEffDoliTepPeren(kc, ksf, pc); tc = po[q]; delete[]po; //при 373 К
		tcc = kc*(1e0 - pc) / pc; fc = (1e0 - tc)*kc / (kc + tcc) - por;
		ka = hps*pa / por;
		po = oprEffDoliTepPeren(ka, ksf, pa); ta = po[q]; delete[]po; //определяем долю площади сечения перемычки
		tca = ka*(1e0 - pa) / pa; fa = (1e0 - ta)*ka / (ka + tca) - por;
		kb = hps*pb / por;
		po = oprEffDoliTepPeren(kb, ksf, pb); tb = po[q]; delete[]po; //через перемычку тепло распространяется чистой теплопроводностью
		tcb = kb*(1e0 - pb) / pb; fb = (1e0 - tb)*kb / (kb + tcb) - por;
		if ((fc*fb > 0.0) && (fa*fc<0.0)) pb = pc; if ((fc*fa>0.0) && (fb*fc < 0.0)) pa = pc;
		j++; ra = fabs(pa - pb); } 
	dpctsm[q]=tc; cout << "Dol Plo CTP = " << tc << endl; }
	for (j=0; j<cem; j++) cout << "EC KTP TK ( " << j << " ) = " << ecktpsha[j] << endl;
	if (!vy) { pkc[0]=kc; return pkc; } //скорректированное значение толщины воздушной прослойки (размер поры), когда ввели перемычку
	else if (vy == 1) return dpctsm;
} //доля площади, через которую происходит перенос тепла чистой теплопроводностью
double *opredPolTempTvKarShaFragm(double *Tna, int n, double *ete, double *ktptk, int ce, double htk, double hvp, double ptpob, double tenac, double zn)
{
	int k=0; double lamt, hh = hvp + htk, ht = htk + hvp / 2e0, tt; Tna[k] = tenac;
	for (k = 1; k < n; k++)  {
		tt = Tna[k - 1];
		lamt = opredKTPTKTochSha(ktptk, ete, tt, ce); //ищем КТП твердого скелета при данной температуре
		Tna[k] = tt + zn*ptpob*ht / lamt;  ht = ht + hh;
	} return Tna;
}  //рассчитываем температуру на границах стенок
double RasIzlSerStenNac(double *Refl, double *Tran, double *Abso, double *mpo, double *Reflt, double *Trant, double *Absot, double temc, double hlr0, double hrl0, double b, int vyve, double poristo, double dpct, double htch, double hvoz, int kost, double *ete, double *kttkve, int cemve, int vtve, double *qobve, double *ktpvove, double *tevove, char *snm, int dmao, double *stchsrve, double *Tpctve, int dmkvove)
{
	int j=2, k=0, f=6, q=0, *nost=new int[j];
	double p=0.0, r=p, hf=1e0, ra=p, qs=p, d=p, **mu=NULL, qobvee=p, *Elr=NULL, p1=hf, p2=hf, *Evsv=NULL;
	for (j=0; j<kost; j++) d=d+hf; //cout << "d = " << d << endl;
	q=2*kost; double *tvs=new double[f], *tbs=new double[q], *epst=new double[f], *hvi=new double[f];
	double *Rs=new double[f], *Ap=new double[f], rsf=0.0;
	if ((!tvs) || (!tbs) || (!epst) || (!hvi) || (!Rs) || (!Ap)) { cout << snm << endl; k = getchar(); exit(1); }
	for (j = 0; j < f; j++) { tvs[j] = 0.0; epst[j] = 0.0; hvi[j] = 0.0; Rs[j] = 0.0; Ap[j] = 0.0; }
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
		for (k = 0; k<f; k++) tvs[k] = Tpctve[j]; tvs[0] = tbs[j]; tvs[2] = tbs[j + 1]; //for (k = 0; k < f; k++) Rs[k] = Reflt[j]; Rs[0] = Refl[j]; Rs[2] = Refl[j + 1]; 
		for (k = 0; k<f; k++) epst[k]=Absot[k]; epst[0]=Abso[j]+Tran[j]; epst[2]=Abso[j+1]+Tran[j+1]; //for (k = 0; k<f; k++) epst[k] = opredKTPTKTochSha(stchsrve, ete, tvs[k], cemve); epst[0] = Abso[j]; epst[2] = Abso[j + 1];
		for (k = 0; k<f; k++) Rs[k] = 1e0-epst[k]; //for (k = 0; k<f; k++) epst[k] = 1e0-Rs[k];
		for (k = 0; k<f; k++) hvi[k] = 0.0; hvi[0] = Elr[j]; //hvi[2]=-Elr[j];
		ra = SeryeStenkiRasIzl(w, h, l, tvs, epst, hvi, Rs, Ap, nost, f); //cout << "ra ( " << j << " ) = " << ra << endl; 
		if (fabs(ra)>0.0) { rsf = ra + rsf; r = r + hf; } q = q + 2; 
	}
	delete[]Ts; delete[]Tss; delete[]Tna; delete[]Tsr; delete[]hlr1;
	delete[]hrl1; delete[]reiz; delete[]sislr; delete[]sisrl;
	delete[]Elr; delete[]Evsv; delete[]ao; delete[]tvs; delete[]tbs; delete[]epst;
	delete[]hvi; delete[]Rs; delete[]nost; delete[]Ap; 
	if (fabs(r)>0.0) rsf = rsf / r; else rsf = 0.0; return rsf;
}