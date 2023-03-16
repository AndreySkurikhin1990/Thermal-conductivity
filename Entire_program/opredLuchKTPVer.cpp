#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <iostream>
#include <time.h>
#define vmivmf 0 //выбор метода измерений для фракции 2-0,7 мм: 0 - нестационарный, 1 - стационарный
#define vyfv 1 //выбор фракции: 0 - фракция 2-0,7 мм, 1 - фракция 8-4 мм, 2 - фракция 1,6-0,35 мм
#define vyuv 2 //выбор укладки: 1 - плоскопараллельная, 2 - вертикальная
#define vysv 2 //выбор состояния: 0 - исходное, 1 - после повторных измерений, 2 - после прокаливания при 1000 град С
using namespace std;
const int dmkvover=28, dmkooscv=14, N=13, ks=10, dmkov=4, qtn=8, nxtv=30, nnxtv = 0;
const int vtvn=0, isrp=0, vpkf=0, vpmf=0, cemdu=6, cemdum=2, vybves=1; //0 - старые, 1 - новые
const double detev = 1e2, tev0 = 273.15, ssiv84 = 467.0*1e-3, salv84 = 129.0*1e-3, smgv84 = 282.0*1e-3;
const double ssiv207=458.0*1e-3, salv207 = 14.0*1e-2, smgv207 = 29.0*1e-2, tkv=22.0+tev0, tnoscv=3e2, tnacv=2e2; //dkoscv - максимальное отклонение от экспериментальных данных в литературных источниках;
const double dtoscv = 1e2, tocrasver = 1e-8, templa = 134.0*1e1 + tev0;
const double epsi = 1e-15, y0ver = 3e1*1e-3, pksvv = 0.0, dkoalv = 0.444905; //dkoalv - определение КП из КО
const char nfv = 'A';
//----------------------------------------
double *mgover = NULL, *siover = NULL, *alover = NULL; 
int vtv = 0, sctxv = 0, cel = 4, dkoscvl = 6, cemv = 11, vtvk=cemv; /*выбор температуры вермикулита*/
double *qobv = NULL, *etev = NULL, *txver = NULL, *qxver = NULL, *lxver = NULL, *kxver = NULL;
double *dtxver = NULL, *stchsrver = NULL;
double *ktrv = NULL, *dkoscvm = NULL, *dkoscvt = NULL, hver = 0.0, hvove = 0.0, hkver = 0.0, qobver = 0.0, *ktpvover;
double dkospv = 1.3002913762; //dkospv - дополнительный коэффициент для КО с учетом пористой структуры;
double *kttkv = NULL, *mkov = NULL, *kektpv = NULL, *mtsv = NULL, dpctv = 1e0;
double *Tpctv, tnav = tnacv + tev0, *tkuscv, *kuscv;
double *temrasv = NULL, tnrv, *ooxver = NULL, *dpctvm, ssiv = 368.0*1e-3, salv = 132.0*1e-3, smgv = 217.0*1e-3;
//----------------------------------------
void zadrktVerNac();
double rasotprpovsSha(int *, int *, int, int, int, int, double *, double *);
double otrprovsSha(int, int, double *, double *);
double opredLuchSostVer(double *, double *, int, int, int);
double *koefPribSha(double *, double *, int, double *, char *);
double *zadrkt(int, int, double, int, double, double, int, int, double, int, int, double *);
double opredKTPTKTochSha(double *, double *, double, int);
double *opredKTPTverKarkVerm(double *, double *, double, double, double, double, int, int, double, double, int, double *, double *, double *, int, int, int, int, int);
double **initarrver(int, int, int, int, int);
void vyvodfile(double *, int, int, double, char *);
void zapisvfile(double *, int, char *);
void osvpam();
double **RaschRTAVer(int, double, double, double, int, int, double, int, double, int, int, int);
double epsisredver(double, double *, double *, int, double *, double *, int);
double *KoefPoglRosselNac(double *, int, int, double, double, double, double *, double *, int, double, double, double *, double *, int, int, int);
double **opredTempLPStenSha(double *, double *, double *, double *, double, int, double *, double *, double *, double *, int, int, double, double, double, int, char *);
double **RaschVneshIzluchSha(double *, double *, double *, double *, double *, int, double, double, char *);
double **RasLuchPloTepPot(int, double *, double *, double *, double *, double *, double *, int, double *, double *, double *, char *);
double **RaschSobLuchPlotTepPot(int, double *, double *, double *, double *, double *, double *, double, double, double *, int, char *);
double *opredTempStenShaFragm(double *, int, double *, double *, double *, double *, int, int, double, double, double, double, double);
double LuchKTPChudnovsky(double *, double, int, double);
double **rasPorpoRazVer(double, int, int, int, int, int);
double RaschAlphaTvKarVer();
double *koefoslab(double, double, double, double *, int, double *);
double ReflSredVer(double);
double ReflSredkvi(double, int);
double ReflSredSha(double);
double ReflSreditom(double, int);
double reflsha(double);
double reflver(double);
double oprProcSoderpoPoris(double *, double *, double, int);
double F0_lamT(double);
double *oprRasTemNach(int, int, double *, double *, int, double *, double, double, double, double, double, int);
double *oprEffDoliTepPerenVer(double, double, double);
double *KorrZnachVozdPros(double, double, double, int);
double *opredPolTempTvKarShaFragm(double *, int, double *, double *, int, double, double, double, double, double);
double *FuncRaschIzl(double, double, double, double, int, double *, double, double *, double, double, double, double, int, double *, double *, int, int, char *, double, int, double *, double, double *, double *, int, double *, double *, int, double *, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, int, double, double, int, double *, double *, double, double, double, char *, int);
double KorrZnachVozdProsVermik(double, double, double);
double **vydelPol(int, int, double **, double **, int, int);
double **RaschRTASha(int, double, double, double, double, double, int, int, int, int);
double **opredTempHolGor(double *, double *, int, int, double, int, double **, int, int, double *, int, double);
double novNapMas(int, int, int, int, int, int, int, int, int);
double **oprkoefKTPiskhchao(int, int, double *, double, int, double **, int, int, int);
double **napMasEKTP(int, int, int, double *, int, double, int, int, int, double **, int, int);
double *danPoTemTepl840(double *, double *, int);
double *danPoTemTepl841(double *, double *, int);
double *danPoTemTepl842(double *, double *, int);
double *danPoTemTepl843(double *, double *, int);
double *danPoTemTepl844(double *, double *, int);
double *danPoTemTepl845(double *, double *, int);
double *danPoTemTepl2071(double *, double *, int);
double *danPoTemTepl2072(double *, double *, int);
double *danPoTemTepl2073(double *, double *, int);
double *danPoTemTepl2074(double *, double *, int);
double **arrTemCold84(double **);
double *arrTemHigh84();
double *arrTepPot84();
double **arrTemCold207(double **);
double *arrTemHigh207();
double *arrTepPot207();
double *danIskh207(double *, double *, int, int, int);
double **PoiskZavVelTemVer(int, double **, int, double *, double *, int);
double **napMasEKTPVerNac(double, double, double, int, int, int, int, int);
double **RaschRTAitom(int, double, double, double, int, int, double, int, double, int, int, int);
double bbfn(double);
double *arrKTP_2020();
double *arrKTP_Netzsch();
double **arrTem_Netzsch(double **);
double *arrTem1_2020();
double **arrTem2_2020(double **);
double *arrTem3_2020();
double **opredSredKTPTK();
double **gotSredKTPTK();
char **napolStrok();
//----------------
double **initarrver(int fl, int vfv, int vuv, int vmi, int vsv)
{
	int j=0, k=0, koel=nxtv, vv=vybves, f=3;
char *snmv=NULL, *sfnov=NULL, *svfdv=NULL, *s=NULL;
char **ms=napolStrok();
k=0; snmv=ms[k]; k++; sfnov=ms[k];
ms=napNazvFile(vyve, vm, snmi, nfi); if (ms) delete[]ms;
k=0; svfdi=ms[k]; k++; j=4; for (q=k; q<j; q++) { s=ms[q]; if (s) delete[]s; } if (ms) delete[]ms;
	double wmg=smgv, wsi=ssiv-pksvv, wal=salv+pksvv, **mu=NULL;
	if ((!vfv) || (vfv==2)) { ssiv=ssiv207; salv=salv207; smgv=smgv207; } else if (vfv==1) { ssiv=ssiv84; salv=salv84; smgv = smgv84; } wsi=ssiv; wal=salv; wmg=smgv;
	dkoscvm = new double[dkoscvl]; dkoscvt = new double[dkoscvl]; 
	if ((!dkoscvm) || (!dkoscvt)) { cout << snmv << endl; k = getchar(); exit(1); }
	double tnd = 6e2, dtd = 2e2, tm=0.0, hf=1e0, ko=1e-2; k=0; dkoscvt[k] = tnd; for (k = 1; k < dkoscvl; k++) dkoscvt[k] = dkoscvt[k - 1] + dtd;
	k = 0; dkoscvm[k] = 4.68; k++; dkoscvm[k] = 4.69; k++; dkoscvm[k] = 5.65; k++; 
	dkoscvm[k] = 13.17; k++; dkoscvm[k] = 20.2; k++; dkoscvm[k] = 27.81;
	for (k = 0; k < dkoscvl; k++) {
		tm = dkoscvm[k]*ko; dkoscvm[k] = hf - tm; //cout << "tem = " << dkoscvt[k] << "\tdkosc = " << dkoscvm[k] << endl; 
	} 
	etev = new double[cemv]; mkov = new double[cemv]; ektpv = new double[cemv]; stchsrver = new double[cemv];
	if ((!alphaVer) || (!Tpctv)) { cout << snmv << endl; j = getchar(); exit(1); }
	if ((!etev) || (!mkov) || (!ektpv) || (!stchsrver)) { cout << snmv << endl; j = getchar(); exit(1); }
	for (j = 0; j < cemv; j++) { etev[j] = 0.0; mkov[j] = 0.0; ektpv[j] = 0.0; stchsrver[j] = 0.0; }
	for (j = 0; j < ks; j++) {
		Aver[j] = 0.0; Rver[j] = 0.0; Raver[j] = 0.0; Taver[j] = 0.0;
		Aaver[j] = 0.0; Rtver[j] = 0.0; Ttver[j] = 0.0; Atver[j] = 0.0; Tver[j] = 0.0; }
	txver = new double[koel]; qxver = new double[koel]; lxver = new double[koel]; 
	kxver = new double[koel]; dtxver = new double[koel]; ooxver = new double[koel];
	if ((!txver) || (!qxver) || (!lxver) || (!kxver) || (!dtxver) || (!ooxver)) { cout << snmv << endl; j = getchar(); exit(1); }
	for (j = 0; j < koel; j++) { txver[j] = 0.0; qxver[j] = 0.0; lxver[j] = 0.0; kxver[j] = 0.0; dtxver[j] = 0.0; ooxver[j] = 0.0; }
	tkuscv = new double[dmkooscv]; kuscv = new double[dmkooscv];
	if ((!tkuscv) || (!kuscv)) { cout << snmv << endl; k = getchar(); exit(1); }
	k=0; tkuscv[k] = tnoscv; for (k = 1; k < dmkooscv; k++) tkuscv[k] = tkuscv[k - 1] + dtoscv;
	kuscv = koefoslab(wmg, wsi, wal, tkuscv, dmkooscv, kuscv); //for (k=0; k<dmkooscv; k++) cout << "te = " << tkuscv[k] << "\tku = " << kuscv[k] << "\t"; cout << endl;
	ktpvover = new double[dmkvover]; vtever = new double[dmkvover];
	for (j=0; j<dmkvover; j++) { ktpvover[j]=0.0; vtever[j]=0.0; } if ((!ktpvover) || (!vtever)) { cout << snmv << endl; j = getchar(); exit(1); }
	NapMasVozdSha(ktpvover, vtever, dmkvover);
	j=0; etev[j]=tnav; for (j=1; j<cemv; j++) etev[j]=etev[j-1]+detev;
	porver=novNapMas(vv, j, vmi, vfv, vyuv, vsv, vpmf, vpkf, cemv);
	cout << "por = " << porver << "\twmg = " << wmg << "\twsi = " << wsi << "\twal = " << wal << "\tvyfv = " << vfv << "\tvyuv = " << vuv << "\tvmi " << vmi << "\tvysv = " << vsv << endl;
	for (k = 0; k<cemv; k++) stchsrver[k] = hf; //for (k = 0; k<cemv; k++) { tm = epsisredver(etev[k], tkuscv, kuscv, dmkooscv, dkoscvt, dkoscvm, dkoscvl); stchsrver[k] = tm; } for (k=0; k<cemv; k++) cout << "tem = " << etev[k] << "\tst_ch = " << stchsrver[k] << "\t"; cout << endl;
	mu=napMasEKTPVerNac(wmg, wsi, wal, vfv, vuv, vmi, vsv, fl); //for (k=0; k<cemv; k++) cout << "te = " << etev[k] << "\tktp = " << ektpv[k] << endl; 
	if (fl>0) { if (mkov) delete[]mkov; if (ektpv) delete[]ektpv; if (stchsrver) delete[]stchsrver; if (Tver) delete[]Tver; if (Aver) delete[]Aver; if (Rver) delete[]Rver; if (Raver) delete[]Raver; if (Taver) delete[]Taver; if (Aaver) delete[]Aaver; if (Rtver) delete[]Rtver; if (Ttver) delete[]Ttver; if (Atver) delete[]Atver; if (alphaVer) delete[]alphaVer; if (Tpctv) delete[]Tpctv; if (dkoscvm) delete[]dkoscvm; if (dkoscvt) delete[]dkoscvt; if (txver) delete[]txver; if (qxver) delete[]qxver; if (lxver) delete[]lxver; if (kxver) delete[]kxver; if (dtxver) delete[]dtxver; if (ooxver) delete[]ooxver; if (tkuscv) delete[]tkuscv; if (kuscv) delete[]kuscv; if (ktpvover) delete[]ktpvover; if (vtever) delete[]vtever; }
	return mu;
}
void osvpamver()
{
	delete[]kektpv; delete[]kttkv; delete[]mkov; delete[]qobv; delete[]snv; delete[]ssv;
	delete[]skptv; delete[]svsv; delete[]snmv; delete[]sfnov; delete[]svfdvu; delete[]sfatv; delete[]sfov;
	delete[]svfdv; delete[]etev; delete[]ektpv; delete[]txver; delete[]qxver; delete[]lxver; delete[]kxver; delete[]ooxver;
	delete[]dtxver; delete[]ktpvover; delete[]vtever; delete[]Tpctv; delete[]alphaVer; delete[]tgorv; delete[]tholv;
	delete[]Tver; delete[]Aver; delete[]Rver; delete[]Raver; delete[]Taver; delete[]Aaver; delete[]dkoscvm; delete[]dkoscvt;
	delete[]Rtver; delete[]Ttver; delete[]Atver; delete[]stchsrver; 
}
void zadrktVerNac()
{
	int j=-1, jk = nxtv, jn = nnxtv, k=0, q=0, f = 6, w=0; 
	char **unau=napolStrok();
	double **mu=NULL, *po=NULL, *temtk=NULL, nf=0.0; 
	mu=initarrver(j, vyfv, vyuv, vmivmf, vysv); //if (mu) delete[]mu; mu=gotSredKTPTK(); //mu=opredSredKTPTK();
	k=0; kttkv=mu[k]; k++; temtk=mu[k]; k++; po=mu[k]; k=0; nf=po[k]; if (mu) delete[]mu; //k=getchar(); 
	double hf=1e0, ka=0.0, kb=0.0, *ras=NULL, *srra=NULL, *legr=NULL, *prgr=NULL, nnxtfv=0.0, nxtfv=0.0, e=1e-10;
	for (j=0; j<jn; j++) nnxtfv=nnxtfv+hf; for (j=0; j<jk; j++) nxtfv=nxtfv+hf; 
	double dhk = y0ver / fabs(nxtfv - nnxtfv), hnver = nnxtfv*dhk, ko=1e-6, hvko = (13e1)*ko, srp, marp;
	double hvh = ko, hvna=0.0, p=0.0, r=0.0, d=0.0, *atr=NULL, t=0.0; 
	for (j = 0; j < ks; j++) d = d + hf; //cout << "d = " << d << "\t"; //dkospv = RaschAlphaTvKarVer(); 
	vtvk=cemv;
	for (vtv = vtvn; vtv<vtvk; vtv++) { //пробегаем по температуре
mu=rasPorpoRazVer(porver, vyfv, 1, vysv, isrp, vpkf); //0 - старые, 1 - новые значения
k=0; ras=mu[k]; k++; srra=mu[k]; k++; prgr=mu[k]; k++; legr=mu[k]; k++; po=mu[k]; 
j=0; srp=po[j]; k++; if (po) delete[]po; po=mu[k]; marp=po[j]; if (po) delete[]po; if (mu) delete[]mu; cout << "Sred razm por = " << ko << "\t";
r=d; t=ko; po=KorrZnachVozdPros(ko, r, porver, 0); j=0; ko=po[j]; j++; p=po[j]; if (po) delete[]po; cout << "Korrektirovannyi sredniy razmer por = " << ko << "\tDolya Ploschadi = " << p << "\t";
ko=t; ko=KorrZnachVozdProsVermik(ko, r, porver); hvna = ko; hvko = ko; k=getchar();
		if (vtv>vtvn) {
			q = jk - jn; ka = (tholv[vtv] - tgorv[vtv]) / y0ver; kb = tgorv[vtv]; hkver = hnver;
			for (k = 0; k < q; k++) { temrasv[k] = kb + ka*hkver; hkver = hkver + dhk; }
		}
		hvove = hvna; hver = hvove*(hf - porver) / porver;
		while (hvove <= hvko) {
			j = jn; hkver = hnver*dhk; //пробегаем по размерам пор
			while ((j < jk) && (hkver < y0ver)) {
				cout << "hk = " << hkver << endl; sctxv = j; //пробегаем по координате
				atr = zadrktVer(j, ks, d, vtv, hver, hvove, 1, 0, etev[vtv], 0, 0, atr);
				kxver[j - jn] = hkver; hkver = hkver + dhk; j++; delete[]atr;
			}
			vyvodfile(lxver, jk - jn, 0, hvove, sfov);
			for (j = 0; j<jk - jn - 1; j++) {
				p = kxver[j + 1] - kxver[j];
				r = txver[j] - txver[j + 1]; if (fabs(p)>e) r = fabs(r / p); else r = 0.0;
				p = fabs(qxver[j + 1] + qxver[j]) / 2e0; lxver[j] = p / r;
			}
			vyvodfile(ooxver, jk - jn, 2, hvove, sfov);
			vyvodfile(txver, jk - jn, 2, hvove, sfov);
			vyvodfile(lxver, jk - jn - 1, 2, hvove, sfov);
			hvove = hvove + hvh;
		}
	}
}
//---------
double **PoiskZavVelTem(int v, double **mu, int rmu, int n, double *efte, double h)
{
	int k=1, f=rmu, j=0, q=j, w=j, b=10; double *poin=new double[k], cf=0.0, hf=1e0, nf=cf, tf=cf, e=1e-1, cfp=cf, cft=cf;
	double *temvht=NULL, *temvct=NULL, *tepvt=NULL, *po=NULL, *temhq=NULL, *ts=NULL, ***muu=NULL;
	double *temcq=NULL, *temvs=NULL, *tepv=NULL, *cemt=NULL, *ktpq=NULL, *ktp=NULL, **muv=NULL;
	int vyfrve=0, vysove=0, vymivmf=1, vyukve=1, c=0, nnvyfv=2, nnvysv=3, nnvmivmf=2;
	int nnvyuv=2, kvf=0, jvsv=0, qvmi=0, qvuv=0, d=0;
	double *nvyfv=new double[nnvyfv], *nvysv=new double[nnvysv], *nvmivmf=new double[nnvmivmf];
	double *nvyuv=new double[nnvyuv]; temhq=new double[n]; temcq=new double[n]; temvs=new double[n]; 
	tepv=new double[n]; cemt=new double[n]; ktpq=new double[n]; ktp=new double[n]; muu=new double**[b]; 
	if ((!nvyfv) || (!nvysv) || (!nvmivmf) || (!nvyuv) || (!temhq) || (!temcq) || (!temvs) || (!tepv) || (!cemt) || (!ktpq) || (!ktp) || (!muu)) 
	{ cout << snmv << endl; k=getchar(); exit(1); } 
	for (k=0; k<n; k++) { temhq[k]=0.0; temcq[k]=0.0; temvs[k]=0.0; tepv[k]=0.0; cemt[k]=0.0; ktpq[k]=0.0; }
	k=0; nvyfv[k]=0; k++; nvyfv[k]=1; //фракции вермикулита
	k=0; nvysv[k]=0; k++; nvysv[k]=1; k++; nvysv[k]=2; //состояния вермикулита
	k=0; nvmivmf[k]=1; k++; nvmivmf[k]=2; //стационарные методы измерений - 2019 и 2020
	k=0; nvyuv[k]=1; k++; nvyuv[k]=2; k=0; j=-1; b=0; //укладка вермикулита
	for (kvf=0; kvf<nnvyfv; kvf++) {
		vyfrve=nvyfv[kvf];
		for (jvsv=0; jvsv<nnvysv; jvsv++) {
			vysove=nvysv[jvsv];
			if ((!vyfrve) || (vyfrve==2)) {
				for (qvmi=0; qvmi<nnvmivmf; qvmi++) {
					vymivmf=nvmivmf[qvmi]; 
					muv=new double*[f]; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
					muv=napMasEKTP(vyfrve, vysove, vymivmf, efte, k, h, vyukve, n, k, muv, f, dmkov); muu[b]=muv; b++; if (vysove>0) break; } }
			else if (vyfrve==1) {
				for (qvuv=0; qvuv<nnvyuv; qvuv++) {
					vyukve=nvyuv[qvuv];
					muv=new double*[f]; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
					muv=napMasEKTP(vyfrve, vysove, vymivmf, efte, k, h, vyukve, n, k, muv, f, dmkov); muu[b]=muv; b++; } } } } 
	for (j=0; j<b; j++) { muv=muu[j];
	d=0; temvct=muv[d]; d++; temvht=muv[d]; d++; tepvt=muv[d]; d++; ktp=muv[d]; d++; ts=muv[d]; d++; po=muv[d]; k=0; nf=po[k]; 
					cf=e; while (cf<nf) { cf=cf+hf; k++; } q=k; 
					d=1; k=0; cfp=tepvt[k]; for (k=d; k<q; k++) { cft=tepvt[k]; if ((cft<=cfp) && (d>0)) { d=-1; break; } cfp=cft; } //for (k=0; k<q; k++) cout << "c = " << c << "\tqo = " << tepvt[k] << "\tts = " << ts[k] << "\t"; //vyvodfile(ts, q, k, nf, "C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\tmp.txt");
						if (d>0) { for (w=0; w<q; w++) { 
							tf=ts[w]; cf=tepvt[w]*tf;
							for (k=0; k<n; k++) {
							if (fabs(tf-efte[k])<=hf) {
								temvs[k]=temvs[k]+tf;
								temhq[k]=temhq[k]+tf*temvht[k]; temcq[k]=temcq[k]+tf*temvct[k]; ktpq[k]=ktpq[k]+tf*ktp[k];
								tepv[k]=tepv[k]+cf; cemt[k]=cemt[k]+hf; //cout << "tf = " << tf << "\tptp = " << tepvt[k] << "\t";
								break; } } } } c++; }
	b=c; for (k=0; k<b; k++) { muv=muu[k]; for (j=0; j<f; j++) { po=muv[j]; delete[]po; } if (muv) delete[]muv; } //for (k=0; k<n; k++) cout << "cemt = " << cemt[k] << "\t";
		for (k=0; k<n; k++) {
			cf=temvs[k]; if (cf>e) {
				temhq[k]=temhq[k]/cf; temcq[k]=temcq[k]/cf; ktpq[k]=ktpq[k]/cf; tepv[k]=tepv[k]/cf; }
			cf=cemt[k]; if (cf>e) temvs[k]=temvs[k]/cf; else temvs[k]=0.0; }
	cf=0.0; for (k=0; k<n; k++) cf=cf+hf;
	k=0; poin[k]=cf;
	mu[k]=temcq; k++; mu[k]=temhq; k++; mu[k]=tepv; k++; 
	mu[k]=ktpq; k++; mu[k]=temvs; k++; mu[k]=poin; 
	if (nvyfv) delete[]nvyfv; if (nvysv) delete[]nvysv;
	if (nvmivmf) delete[]nvmivmf; if (nvyuv) delete[]nvyuv;
	if (muu) delete[]muu; if (cemt) delete[]cemt;
	return mu; 
}
double **opredTempHolGor(double *ektp, double *ete, int n, int l, double h0, int v, double **mu, int rmu, int ni, double *efte, int dmkoef, double tn) //моделирование процесса теплообмена в образце
{ //n - длина массива ektp, l - длина массивов temvs, qob, ktp, temvh, temvc, ni - длина efte, qob - плотность теплового потока, которую может создать лабораторная установка 
	int cemf=rmu, k=0, j=1, w=0, rm=k, nk=n; 
	double g=0.0, p=0.0, *qon=NULL, *tena=NULL, nit=1e7, hf=1e0, *po=new double[j], del=hf;
	double ep=1e-2, d=1e-3, Thna=tn, Thnac=0.0, *koeq=new double[dmkoef], ts=0.0, laef=0.0;
	double tgor=0.0, thol=0.0, tsred=0.0, etem=0.0, dt=hf, **muv=new double*[cemf], nf=0.0;
	double *qob=new double[n], *temvs=new double[n], *temvc=new double[n], *temvh=new double[n], *ktp=new double[n];
	if ((!koeq) || (!po) || (!muv) || (!qob) || (!temvc) || (!temvh) || (!temvs)) { cout << snmv << endl; k=getchar(); exit(1); }
	muv=PoiskZavVelTem(k, muv, cemf, ni, efte, h0);
	k=2; qon=muv[k]; k=4; tena=muv[k]; k++; po=muv[k]; 
	k=0; nf=po[k]; w=k; g=ep; while (g<nf) { g=g+hf; w++; } //cout << "w = " << w << endl; for (k=0; k<w; k++) cout << "q ( " << k << " ) = " << qon[k] << "\tt = " << tena[k] << endl;
	for (k=0; k<dmkoef; k++) koeq[k]=0.0; koeq=koefPribSha(qon, tena, w, koeq, snmv);
		for (k=0; k<ni; k++) {
			ts=ete[k]; g=0.0; p=g;
			for (j=0; j<dmkoef; j++) { g=g+pow(ts, p)*koeq[j]; p=p+hf; } 
			qob[k]=g; } 
		if (koeq) delete[]koeq; for (k=0; k<cemf; k++) { po=muv[k]; delete[]po; } delete[]muv;
		g=0.0; for (k=0; k<n; k++) {
			if (qob[k]>ep) {
				laef=ektp[k]; p=0.0; Thnac=Thna+g*dt; del=hf; etem=ete[k]; ktp[k]=laef;
				while ((del>ep) && (p<nit)) {
					thol=Thnac+p*d; //Tg - массив температур горячих стенок, Th - массив температур холодных стенок
					tgor=thol+qob[k]*h0/laef;
					del=(2e0*etem-(thol+tgor));
					p=p+hf; } 
				g=g+hf; }
			else { thol=0.0; tgor=0.0; qob[k]=0.0; ktp[k]=0.0; } tsred=(tgor+thol)/2e0; 
			temvc[k]=thol; temvh[k]=tgor; temvs[k]=tsred; }  
	g=0.0; for (k=0; k<n; k++) g=g+hf; //cout << "n = " << g << endl;
	k=1; po=new double[k]; muv=new double*[cemf]; if ((!po) || (!muv)) { cout << snmv << endl; k=getchar(); exit(1); } k=0; po[k]=g;
	muv[k]=temvc; k++; muv[k]=temvh; k++; muv[k]=qob; k++; muv[k]=ktp; k++; muv[k]=temvs; k++; muv[k]=po; //for (k=0; k<cemf; k++) mu[k]=muv[k]; cout << "n = " << n << endl; //for (k=0; k<n; k++) cout << "qob = " << qob[k] << "\tektp = " << ektp[k] << "\tte = " << ete[k] << "tc = " << temvc[k] << "\ttg = " << temvh[k] << "\tts = " << temvs[k] << endl; 
	k=0; j=1; mu=vydelPol(k, n, muv, mu, rmu, j); 
	if (temvc) delete[]temvc; if (temvh) delete[]temvh; //if (qob) delete[]qob; 
	if (ktp) delete[]ktp; if (temvs) delete[]temvs; if (po) delete[]po; if (muv) delete[]muv; 
	return mu;
}
double **napMasEKTPVerNac(double wmg, double wsi, double wal, int vfv, int vuv, int vmi, int vsv, int fl)
{
	int k=1, nk=0, j=0, jk=0, f=cemdu, qn=0, q=0, c=3; 
	double hf=1e0, *po=NULL, s=0.0, *etevv=NULL, *isv=NULL, ys0=y0ver;
	double **mu=new double*[f], **muv=new double*[f]; if ((!mu) || (!muv)) { cout << snmv << endl; k=getchar(); exit(1); } 
	double nf=0.0, t=0.0, g=0.0, e=1e-6, tn=0.0;
	muv=napMasEKTP(vfv, vsv, vmi, etev, k, ys0, vuv, cemv, k, muv, f, dmkov);
	k=0; j=1; mu=vydelPol(k, cemv, muv, mu, f, j);
	etevv=new double[cemv]; isv=new double[cemv];
	if ((!etevv) || (!isv)) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k=0; k<cemv; k++) { etevv[k]=etev[k]; isv[k]=stchsrver[k]; }
	k=0; tholv=mu[k]; k++; tgorv=mu[k]; k++; qobv=mu[k]; k++; 
	ektpv=mu[k]; k++; etev=mu[k]; k++; po=mu[k]; 
	if (muv) delete[]muv; if (mu) delete[]mu; 
	jk=cemv; k=0; nf=po[k]; if (po) delete[]po; 
	t=e; while (t<nf) { t=t+hf; k++; } cemv=k; //cout << "n = " << cemv << "\t"; //for (k=0; k<cemv; k++) cout << "qob = " << qobv[k] << "\tektp = " << ektpv[k] << "\tte = " << etev[k] << "\ttc = " << tholv[k] << "\ttg = " << tgorv[k] << "\tts = " << etev[k] << endl;
	k=0; tn = etev[k]; tnav = tn; 
	if (stchsrver) delete[]stchsrver; kttkv = new double[cemv]; stchsrver = new double[cemv]; 
	if ((!stchsrver) || (!kttkv)) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < cemv; k++) { stchsrver[k] = 0.0; kttkv[k] = 0.0; }
	for (j=0; j<cemv; j++) { f=1; nk=0; for (k=0; k<jk; k++) 
		if (fabs(etev[j]-etevv[k])<hf) { f=-1; nk=k; break; } 
		if (f<0) { stchsrver[j]=isv[nk]; } 
		else { stchsrver[j]=0.0; } } 
	for (k = 0; k<cemv; k++) { g = epsisredver(etev[k], tkuscv, kuscv, dmkooscv, dkoscvt, dkoscvm, dkoscvl); stchsrver[k] = g; } 
	if (etevv) delete[]etevv; if (isv) delete[]isv; //for (k=0; k<cemv; k++) cout << "tem = " << etev[k] << "\tst_ch = " << stchsrver[k] << "\tmko = " << mkov[k] << "\t"; cout << endl;
	kttkv = opredKTPTverKarkVerm(etev, ektpv, porver, wsi, wal, wmg, vyfv, dmkov, tnoscv, dtoscv, dmkooscv, kuscv, tkuscv, stchsrver, cemv, vysv, vpkf, isrp, fl); 
	for (k=0; k<cemv; k++) cout << "tem = " << etev[k] << "\tkttkv = " << kttkv[k] << "\tektp = " << ektpv[k] << endl;
	k=1; po=new double[k]; k=0; po[k]=nf; 
	mu=new double*[c]; 
	k=0; mu[k]=kttkv; k++; mu[k]=etev; k++; mu[k]=po;
	if (fl>0) { if (tholv) delete[]tholv; if (tgorv) delete[]tgorv; if (qobv) delete[]qobv; }
	return mu;
}
double **napMasEKTP(int vyfrve, int vysove, int vymeizvemafr, double *te, int v, double h, int vyukve, int nt, int wa, double **mu, int rmu, int dmk)
{
	int k=0, nvm=dmk, j=0, vy=0, n84=0, n207=0, cedu=rmu, q=0, cedumi=cemdum, nn=0;
	double F=13.85*1e-4, nf84=0.0, nf207=0.0, *po=NULL, *tc84=NULL, s=0.0, t=0.0, e=1e-1; 
	double *th84=NULL, hf=1e0, *tp84=NULL, *th207=NULL, *tc207=NULL, *tp207=NULL;
	double *koeft=NULL, *koefh=NULL, *koefc=NULL, *koefq=NULL, *koefs=NULL, **muv=NULL, r=0.0; //cout << "vf = " << vyfrve << "\tvs = " << vysove << "\tvmi = " << vymeizvemafr << "\tvu = " << vyukve << endl; //for (k=0; k<nt; k++) cout << "te ( " << k << " ) = " << te[k] << "\t";
	double *ktp=NULL, *temvs=NULL, *temvh=NULL, *temvc=NULL, *tepv=NULL, *vm=NULL, *ma=NULL; 
	//-------
	if ((!vyfrve) || (vyfrve==2)) { //фракция 2-0,7 мм или фракция 1,6-0,35 мм
		th207=arrTemHigh207(); 
		k=cedumi; muv=new double*[k]; if (!muv) { cout << snmv << endl; k = getchar(); exit(1); }
		muv=arrTemCold207(muv); k=0; po=muv[k]; k++; tc207=muv[k]; 
		tp207=arrTepPot207();
		k=0; t=0.0; nf207=po[k]; while (t<nf207) { t=t+hf; k++; } n207=k; nn=k; //cout << "n207 = " << n207 << "\t";
		temvs=new double[n207]; temvh=new double[n207]; temvc=new double[n207]; 
		tepv=new double[n207]; ktp=new double[n207]; 
		if ((!temvs) || (!temvh) || (!temvc) || (!tepv) || (!ktp)) { cout << snmv << endl; k=getchar(); exit(1); } //for (k=0; k<n207; k++) cout << "th = " << th207[k] << "\ttc = " << tc207[k] << "\tQ = " << tp207[k] << endl;
		for (k=0; k<n207; k++) {
		temvh[k]=th207[k]; temvc[k]=tc207[k]; temvs[k]=(temvh[k]+temvc[k])/2e0; 
		tepv[k]=tp207[k]/F; ktp[k]=fabs(temvh[k]-temvc[k])/h; ktp[k]=tepv[k]/ktp[k]; }
		if (tc207) delete[]tc207; if (th207) delete[]th207; if (tp207) delete[]tp207; if (po) delete[]po; 
		if (!vysove) {
		if (!vymeizvemafr) { //установка Netzsch
		if (muv) delete[]muv; muv=new double*[cedu]; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
		muv=oprkoefKTPiskhchao(vymeizvemafr, k, te, h, n207, muv, cedu, nt, dmk); 
		k=0; koefc=muv[k]; k++; koefh=muv[k]; k++; koefq=muv[k]; k++;
		koeft=muv[k]; k++; koefs=muv[k]; k++; po=muv[k]; }
		else if (vymeizvemafr==1) { //данные 2020 года - ГОСТ 12170 - стационарный метод
		if (ktp) delete[]ktp; ktp=arrKTP_2020(); if (temvh) delete[]temvh; temvh=arrTem1_2020(); 
		k=cedumi; muv=new double*[k]; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
		muv=arrTem2_2020(muv);
		k=0; po=muv[k]; nf207=po[k]; k++; if (temvc) delete[]temvc; temvc=muv[k]; 
		if (temvs) delete[]temvs; temvs=arrTem3_2020(); 
		k=0; s=e; while (s<nf207) { s=s+hf; k++; } n207=k; k=0; q=3; //cout << "nk 2020 = " << n207 << endl; for (k=0; k<n207; k++) cout << "ts ( " << k << " ) = " << temvs[k] << "\tktp = " << ktp[k] << endl;
		double *ts=new double[n207], *tst=new double[n207]; if ((!ts) || (!tst)) { cout << snmv << endl; k=getchar(); exit(1); }
		for (k=0; k<n207; k++) { tst[k]=(temvh[k]+temvc[k])/2e0; ts[k]=tst[k]; temvs[k]=ts[k]; }
		ts=danIskh207(ts, tst, k, q, q); ktp=danIskh207(ktp, tst, k, q, q); temvh=danIskh207(temvh, tst, k, q, q); 
		temvc=danIskh207(temvc, tst, k, q, q); temvs=danIskh207(temvs, tst, k, q, q); if (tst) delete[]tst; //for (k=0; k<q; k++) cout << "ts ( " << k << " ) = " << ts[k] << "\tktp = " << ktp[k] << endl;
		koeft=new double[dmk]; koefh=new double[dmk];
		koefc=new double[dmk]; koefs=new double[dmk]; koefq=new double[dmk];
		if ((!koeft) || (!koefh) || (!koefc) || (!koefs) || (!koefq)) { cout << snmv << endl; k=getchar(); exit(1); }
		for (k=0; k<dmk; k++) { koefh[k]=0.0; koefc[k]=0.0; koefq[k]=0.0; koeft[k]=0.0; koefs[k]=0.0; } 
		if (ktp) koeft=koefPribSha(ktp, ts, q, koeft, snmv); 
		if (temvh) koefh=koefPribSha(temvh, ts, q, koefh, snmv);
		if (temvc) koefc=koefPribSha(temvc, ts, q, koefc, snmv);
		if (temvs) koefs=koefPribSha(temvs, ts, q, koefs, snmv);
		tepv=new double[q]; if (!tepv) { cout << snmv << endl; k=getchar(); exit(1); }
		for (k=0; k<q; k++) {
			s=0.0; r=0.0; for (j=0; j<dmk; j++) { s=s+koeft[j]*pow(ts[k], r); r=r+hf; }
			tepv[k]=s*fabs(temvh[k]-temvc[k])/h; }
			koefq=koefPribSha(tepv, ts, q, koefq, snmv); }
		else if (vymeizvemafr==2) { //данные 2019 года - ГОСТ 12170
		koefq=danPoTemTepl2071(temvs, tepv, nvm); koefh=danPoTemTepl2071(temvs, temvh, nvm);
		koefc=danPoTemTepl2071(temvs, temvc, nvm); koeft=danPoTemTepl2071(temvs, ktp, nvm);
		koefs=danPoTemTepl2071(temvs, temvs, nvm); } }
				if (vysove==1) { //после повторных измерений
			koefc=danPoTemTepl2072(temvs, temvc, nvm); koefh=danPoTemTepl2072(temvs, temvh, nvm); 
			koefq=danPoTemTepl2072(temvs, tepv, nvm); koeft=danPoTemTepl2072(temvs, ktp, nvm); 
			koefs=danPoTemTepl2072(temvs, temvs, nvm); }
				else if (vysove == 2) { //после обжига при 1000 °С
			koefc=danPoTemTepl2073(temvs, temvc, nvm); koefh=danPoTemTepl2073(temvs, temvh, nvm); 
			koefq=danPoTemTepl2073(temvs, tepv, nvm); koeft=danPoTemTepl2073(temvs, ktp, nvm); 
			koefs=danPoTemTepl2073(temvs, temvs, nvm); }
				else if (vysove == 3) 
				{ //после повторного обжига при 1000 °С
			koefc=danPoTemTepl2074(temvs, temvc, nvm); koefh=danPoTemTepl2074(temvs, temvh, nvm); 
			koefq=danPoTemTepl2074(temvs, tepv, nvm); koeft=danPoTemTepl2074(temvs, ktp, nvm); 
			koefs=danPoTemTepl2074(temvs, temvs, nvm); } }
	else if (vyfrve==1) { //фракция 8-4 мм
		k=cedumi; muv=new double*[k]; muv=arrTemCold84(muv); k=0; po=muv[k]; k++; tc84=muv[k]; 
		th84=arrTemHigh84(); tp84=arrTepPot84();
		k=0; t=0.0; nf84=po[k]; while (t<nf84) { t=t+hf; k++; } n84=k; nn=k; //for (k=0; k<n84; k++) cout << "th = " << th84[k] << "\ttc = " << tc84[k] << "\ttp = " << tp84[k] << endl;
		temvs=new double[n84]; temvh=new double[n84]; temvc=new double[n84]; tepv=new double[n84]; ktp=new double[n84];
		if ((!temvs) || (!temvh) || (!temvc) || (!tepv) || (!ktp)) { cout << snmv << endl; k=getchar(); exit(1); }
		for (k=0; k<n84; k++) { 
		temvs[k]=(th84[k]+tc84[k])/2e0; 
		temvh[k]=th84[k]; temvc[k]=tc84[k]; tepv[k]=tp84[k]; 
		tepv[k]=tepv[k]/F; ktp[k]=fabs(th84[k]-tc84[k])/h; 
		ktp[k]=tepv[k]/ktp[k]; }
		if (po) delete[]po; if (th84) delete[]th84; if (tp84) delete[]tp84; if (tc84) delete[]tc84; 
		if (vyukve==1) { //плоско-параллельная засыпка
			if (!vysove) { //исходный
				koefc=danPoTemTepl840(temvs, temvc, nvm); koefh=danPoTemTepl840(temvs, temvh, nvm); 
				koefq=danPoTemTepl840(temvs, tepv, nvm); koeft=danPoTemTepl840(temvs, ktp, nvm); 
				koefs=danPoTemTepl840(temvs, temvs, nvm); }
			else if (vysove==1) { //после повторных измерений
				koefc=danPoTemTepl842(temvs, temvc, nvm); koefh=danPoTemTepl842(temvs, temvh, nvm); 
				koefq=danPoTemTepl842(temvs, tepv, nvm); koeft=danPoTemTepl842(temvs, ktp, nvm); 
				koefs=danPoTemTepl842(temvs, temvs, nvm); }
			else if (vysove==2) { //после обжига
				koefc=danPoTemTepl845(temvs, temvc, nvm); koefh=danPoTemTepl845(temvs, temvh, nvm); 
				koefq=danPoTemTepl845(temvs, tepv, nvm); koeft=danPoTemTepl845(temvs, ktp, nvm); 
				koefs=danPoTemTepl845(temvs, temvs, nvm); } }
		else if (vyukve==2) { //вертикальная засыпка
			if (!vysove) { //исходный
				koefc=danPoTemTepl841(temvs, temvc, nvm); koefh=danPoTemTepl841(temvs, temvh, nvm); 
				koefq=danPoTemTepl841(temvs, tepv, nvm); koeft=danPoTemTepl841(temvs, ktp, nvm); 
				koefs=danPoTemTepl841(temvs, temvs, nvm); }
			else if (vysove == 1) { //после повторных измерений
				koefc = danPoTemTepl844(temvs, temvc, nvm); koefh = danPoTemTepl844(temvs, temvh, nvm); 
				koefq = danPoTemTepl844(temvs, tepv, nvm); koeft = danPoTemTepl844(temvs, ktp, nvm); 
				koefs = danPoTemTepl844(temvs, temvs, nvm); }
			else if (vysove == 2) { //после обжига
				koefc = danPoTemTepl843(temvs, temvc, nvm); koefh = danPoTemTepl843(temvs, temvh, nvm); 
				koefq = danPoTemTepl843(temvs, tepv, nvm); koeft = danPoTemTepl843(temvs, ktp, nvm); 
				koefs = danPoTemTepl843(temvs, temvs, nvm); } } }
	if (temvs) delete[]temvs; if (temvh) delete[]temvh; 
	if (temvc) delete[]temvc; if (tepv) delete[]tepv; 
	if (ktp) delete[]ktp; if (muv) delete[]muv; 
	t=0.0; for (k=0; k<cedu; k++) t=t+hf; 
	k=1; muv=new double*[cedu]; po=new double[k]; if ((!muv) || (!po)) { cout << snmv << endl; k=getchar(); exit(1); }
	k=0; po[k]=t; muv[k]=koefc; k++; muv[k]=koefh; k++; 
	muv[k]=koefq; k++; muv[k]=koeft; k++; muv[k]=koefs; k++; muv[k]=po; //vy=cedu-1; for (k=0; k<vy; k++) { ma=muv[k]; for (j=0; j<dmk; j++) cout << "k = " << k << "\tkoef ( " << j << " ) = " << ma[j] << "\t"; cout << endl; }
	q=cedu-1; for (vy=0; vy<q; vy++) { 
	vm=muv[vy]; ma=new double[nt]; if (!ma) { cout << snmv << endl; k=getchar(); exit(1); }
	for (k=0; k<nt; k++) { t=te[k]; s=0.0; r=s;
	for (j=0; j<nvm; j++) { s=s+vm[j]*pow(t, r); r=r+hf; } 
	ma[k]=s; //if (!vy) cout << "te_s = " << t << "\tt_c = " << s << "\t"; if (vy==1) cout << "te_s = " << t << "\tt_h = " << s << "\t"; if (vy==2) cout << "te_s = " << t << "\tq = " << s << "\t"; if (vy==3) cout << "te_s = " << t << "\tktp = " << s << "\t";if (vy==4) cout << "te_s = " << t << "\tts = " << s << "\t";
	} //cout << endl;
	mu[vy]=ma; if (vm) delete[]vm; } po=muv[q]; if (po) delete[]po; if (muv) delete[]muv;
	k=1; po=new double[k]; if (!po) { cout << snmv << endl; k=getchar(); exit(1); }
	s=0.0; for (k=0; k<nt; k++) s=s+hf; k=0; po[k]=s; mu[q]=po; //cout << "nt = " << nt << "\t";
	muv=new double*[cedu]; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
	k=0; j=-1; muv=vydelPol(k, nt, mu, muv, rmu, j); for (k=0; k<cedu; k++) { po=muv[k]; mu[k]=po; } //k=5; po=mu[k]; k=0; t=po[k]; //cout << "s = " << s << "\tt = " << t << "\t"; 
	if (muv) delete[]muv; 
	return mu;
}
double **oprkoefKTPiskhchao(int vmiv, int v, double *efte, double h, int n207, double **mu, int rmu, int cem, int dmk) //vmiv - выбор метода измерений
{
	int nk=0, k=0, j=0, q=0, qn=0, nn=n207, cedumi=cemdum, w=rmu-1, f=rmu; 
	double hf=1e0, nf207=0.0, *po=NULL, *koeft=NULL, *koefh=NULL, s=0.0, t=0.0;
	double *koefc=NULL, *koefs=NULL, *koefq=NULL, **muv=NULL, **muvv=NULL, r=0.0, e=1e-1;
	double *tsv=NULL, *tgv=NULL, *thv=NULL, *qov=NULL, *ktpv=NULL, *mt=NULL, *ts=NULL; 
	if (!vmiv) {  //0 - установка Netzsch - нестационарный метод
		k=cedumi; muv=new double*[k]; if (!muv) { cout << snmv << endl; k = getchar(); exit(1); }
		muv=arrTem_Netzsch(muv);
		k=0; po=muv[k]; k++; mt=muv[k];
		k=0; s=e; nf207=po[k]; while (s<nf207) { s=s+hf; k++; } nk=k; //cout << "nk = " << nk << "\t"; //nk=8 - длина массива ktpv
		ktpv=arrKTP_Netzsch(); 
		if (muv) delete[]muv; if (po) delete[]po;
		k=rmu; muv=new double*[k]; k=0; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
		muv=opredTempHolGor(ktpv, mt, nk, nn, h, k, muv, rmu, cem, efte, dmk, tkv); //cem - длина массива efte //for (k=0; k<f; k++) { po=muv[k]; delete[]po; }
		if (ktpv) delete[]ktpv; if (mt) delete[]mt;
		k=0; thv=muv[k]; k++; tgv=muv[k]; k++; qov=muv[k]; k++; ktpv=muv[k]; k++; 
		tsv=muv[k]; k++; po=muv[k]; if (muv) delete[]muv;
		q=0; nf207=po[q]; s=0.0; while (s<nf207) { s=s+hf; q++; } qn=q; 
		koefh=new double[dmk]; koefc=new double[dmk]; koefq=new double[dmk]; koeft=new double[dmk]; koefs=new double[dmk];
		if ((!koefh) || (!koefc) || (!koefq) || (!koeft) || (!koefs)) { cout << snmv << endl; k=getchar(); exit(1); }
		for (k=0; k<dmk; k++) { koefh[k]=0.0; koefc[k]=0.0; koefq[k]=0.0; koeft[k]=0.0; koefs[k]=0.0; }
		koefc=koefPribSha(thv, tsv, qn, koefc, snmv); koefh=koefPribSha(tgv, tsv, qn, koefh, snmv); 
		koefq=koefPribSha(qov, tsv, qn, koefq, snmv); koeft=koefPribSha(ktpv, tsv, qn, koeft, snmv); 
		koefs=koefPribSha(tsv, tsv, qn, koefs, snmv); 
		if (thv) delete[]thv; if (tgv) delete[]tgv; if (qov) delete[]qov; 
		if (ktpv) delete[]ktpv; if (tsv) delete[]tsv; if (po) delete[]po; 
	}
	k=1; po=new double[k]; if (!po) { cout << snmv << endl; k=getchar(); exit(1); }
	s=0.0; for (k=0; k<rmu; k++) s=s+hf; k=0; po[k]=s;
	k=0; mu[k]=koefc; k++; mu[k]=koefh; k++; mu[k]=koefq; k++; mu[k]=koeft; k++; mu[k]=koefs; k++; mu[k]=po;
	return mu; 
}
double *danIskh207(double *ma, double *x, int v, int n, int np)
{
	if (ma) 
	{
		int k=0, p=0, q=0; 
		double s=0.0, t=0.0, hf=1e0, *m=new double[n];
		if (!m) { cout << snmv << endl; k=getchar(); exit(1); }
		for (k=0; k<np; k++) { s=s+ma[k]*x[k]; t=t+x[k]; } s=s/t;
		k=0; m[k]=s; k++;
		p=3; m[k]=ma[p]; k++;
		p=4; q=5; m[k]=(ma[p]*x[p]+ma[q]*x[q])/(x[p]+x[q]); //600, 800 и 1000 °C
		if (ma) delete[]ma; return m;
	}
	else return ma;
}
//---------
double **arrTem_Netzsch(double **mu)
{
	int k = 0, le = 8; double *tem = new double[le], hf = 1e0, nf = 0.0, *po = new double[1];
	if ((!tem) || (!po)) { cout << snmv << endl; k = getchar(); exit(1); }
	tem[k]=27.967;  k++; tem[k]=93.833;  k++; tem[k]=192.5; k++; 
	tem[k]=341.667; k++; tem[k]=491.467; k++; tem[k]=641.2; k++; 
	tem[k]=790.933; k++; tem[k]=991.133; k++;
	nf=0.0; for (k=0; k<le; k++) { tem[k]=tem[k]+tev0; nf=nf+hf; } k=0; po[k]=nf;
	k=0; mu[k]=po; k++; mu[k]=tem;
	return mu;
}
double *arrKTP_Netzsch()
{
	int k=0, le=8; double *ktp=new double[le]; 
	if (!ktp) { cout << snmv << endl; k = getchar(); exit(1); }
	ktp[k] = 7.0*1e-2;   k++; ktp[k] = 8.0*1e-2;   k++; ktp[k] = 103.0*1e-3; k++; 
	ktp[k] = 15.0*1e-2;  k++; ktp[k] = 207.0*1e-3; k++; ktp[k] = 283.0*1e-3; k++; 
	ktp[k] = 373.0*1e-3; k++; ktp[k] = 477.0*1e-3; 
	return ktp;
}
double *arrKTP_2020()
{
	int k=0, n=6; double *a = new double[n]; if (!a) { cout << snmv << endl; k = getchar(); exit(1); }
	a[k] = 0.175566644058715; k++; a[k] = 0.176801537812368; k++; a[k] = 0.179324717653617; k++;
	a[k] = 0.211768953068592; k++; a[k] = 0.237194543621728; k++; a[k] = 0.237231989760775;
	return a;
}
double *arrTem1_2020()
{
	int k=0, n=6; double *a=new double[n]; if (!a) { cout << snmv << endl; k = getchar(); exit(1); }
	a[k] = 585.0; k++; a[k] = 600.0;  k++; a[k] = 585.0; k++; 
	a[k] = 800.0; k++; a[k] = 1000.0; k++; a[k] = 1000.0;
	for (k = 0; k < n; k++) a[k] = a[k] + tev0;
	return a;
}
double **arrTem2_2020(double **mu)
{
	int k=1, n=6; double *a=new double[n], nf=0.0, hf=1e0, *po=new double[k];
	if ((!a) || (!po)) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < n; k++) nf = nf + hf; k=0; po[k]=nf;
	a[k] = 119.75; k++; a[k] = 138.0; k++; a[k] = 129.5; k++; 
	a[k] = 200.0;  k++; a[k] = 273.0; k++; a[k] = 261.0;
	for (k=0; k<n; k++) a[k]=a[k]+tev0;
	k=0; mu[k]=po; k++; mu[k]=a;
	return mu;
}
double *arrTem3_2020()
{
	int k=0, n=6; double *a=new double[n]; if (!a) { cout << snmv << endl; k = getchar(); exit(1); }
	a[k] = 377.0; k++; a[k] = 396.0; k++; a[k] = 383.5; k++; 
	a[k] = 548.0; k++; a[k] = 703.0; k++; a[k] = 697.25;
	for (k = 0; k < n; k++) a[k] = a[k] + tev0;
	return a;
}
//---------
double *danPoTemTepl840(double *temvs, double *temvh, int n) //Засыпка плоско-параллельная, исходный
{
	double tho1=0.0, tho2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n]; 
	int k=0, p=0, q=0;
	if (!isdan) { cout << snmv << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=14; q=16; tem1=temvs[p]+temvs[q];
	p=15; q=17; tem2=temvs[p]+temvs[q];
	p=14; q=16; tho1=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem1;
	p=15; q=17; tho2=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem2;
	tem1=tem1/2e0; tem2=tem2/2e0;
	double k1=0.0, k2=0.0; 
	if (fabs(tem2-tem1)>0.0) k1=(tho2-tho1)/(tem2-tem1); else k1=0.0; 
	k2=tho2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl841(double *temvs, double *temvh, int n) //Засыпка вертикальная, исходный
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan = new double[n]; 
	int k=0, p=0, q=0, r=0;
	if (!isdan) { cout << snmv << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=4; q=6; r=10; tem1=temvs[p]+temvs[q]+temvs[r];
	p=5; q=7; r=11; tem2=temvs[p]+temvs[q]+temvs[r];
	p=4; q=6; r=10; if (fabs(tem1)>0.0) tc1=(temvh[p]*temvs[p]+temvh[q]*temvs[q]+temvh[r]*temvs[r])/tem1;
	p=5; q=7; r=11; if (fabs(tem2)>0.0) tc2=(temvh[p]*temvs[p]+temvh[q]*temvs[q]+temvh[r]*temvs[r])/tem2;
	tem1=tem1/3e0; tem2=tem2/3e0;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tc2-tc1)/(tem2-tem1); else k1=0.0; k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl842(double *temvs, double *temvh, int n) //Засыпка плоско-параллельная, повторы
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n]; 
	int k=0, p=0, q=0, r=0;
	if (!isdan) { cout << snmv << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=18; q=20; r=22; tem1=temvs[p]+temvs[q]+temvs[r];
	p=19; q=21; r=23; tem2=temvs[p]+temvs[q]+temvs[r];
	p=18; q=20; r=22; if (fabs(tem1)>0.0) tc1=(temvh[p]*temvs[p]+temvh[q]*temvs[q]+temvh[r]*temvs[r])/tem1; else tc1=0.0;
	p=19; q=21; r=23; if (fabs(tem2)>0.0) tc2=(temvh[p]*temvs[p]+temvh[q]*temvs[q]+temvh[r]*temvs[r])/tem2; else tc2=0.0;
	tem1=tem1/3e0; tem2=tem2/3e0;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tc2-tc1)/(tem2-tem1); else k1=0.0; k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl843(double *temvs, double *temvh, int n) //Засыпка вертикальная, после обжига при 1000 °С
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n]; int k=0, p=0, q=0;
	if (!isdan) { cout << snmv << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=0; q=2; tem1=temvs[p]+temvs[q];
	p=1; q=3; tem2=temvs[p]+temvs[q];
	p=0; q=2; if (fabs(tem1)>0.0) tc1=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem1;
	p=1; q=3; if (fabs(tem2)>0.0) tc2=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem2;
	tem1=tem1/2e0; tem2=tem2/2e0;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tc2-tc1)/(tem2-tem1); else k1=0.0; k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl844(double *temvs, double *temvh, int n) //Засыпка вертикальная, повторы
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n], k1=0.0, k2=0.0; 
	int k=0, p=0, q=0;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=8; q=12; tem1=temvs[p]+temvs[q];
	p=9; q=13; tem2=temvs[p]+temvs[q];
	p=8; q=12; if (fabs(tem1)>0.0) tc1=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem1;
	p=9; q=13; if (fabs(tem2)>0.0) tc2=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem2;
	tem1=tem1/2e0; tem2=tem2/2e0;
	if (fabs(tem2-tem1)>0.0) k1=(tc2-tc1)/(tem2-tem1); else k1=0.0; k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl845(double *temvs, double *temvh, int n) //Засыпка плоско-параллельная, после обжига при 1000 °С
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n], hf=1e0; int k=0;
	if (!isdan) { cout << snmv << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	k=24; tem1=temvs[k];
	k=25; tem2=temvs[k];
	k=24; if (fabs(tem1)>0.0) tc1=(temvh[k]*temvs[k])/tem1;
	k=25; if (fabs(tem2)>0.0) tc2=(temvh[k]*temvs[k])/tem2;
	tem1=tem1/hf; tem2=tem2/hf;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tc2-tc1)/(tem2-tem1); else k1=0.0; k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
//---------
double *danPoTemTepl2071(double *temvs, double *tepv, int n) //Засыпка исходная, фракция 2-0,7 мм
{
	double tepv1=0.0, tepv2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n], hf=1e0; int k=0;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k] = 0.0;
	k=0; tem1=temvs[k];
	k=1; tem2=temvs[k];
	k=0; if (fabs(tem1)>0.0) tepv1=(tepv[k]*temvs[k])/tem1;
	k=1; if (fabs(tem2)>0.0) tepv2=(tepv[k]*temvs[k])/tem2;
	tem1=tem1/hf; tem2=tem2/hf;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tepv2-tepv1)/(tem2-tem1); else k1=0.0; k2=tepv2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1;
	return isdan;
}
double *danPoTemTepl2072(double *temvs, double *tepv, int n) //Фракция 2-0,7 мм (повторные измерения)
{
	double tepv1=0.0, tepv2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n]; 
	int k=0, p=0, q=0;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=2; q=4; tem1=temvs[p]+temvs[q];
	p=3; q=5; tem2=temvs[p]+temvs[q];
	p=2; q=4; if (fabs(tem1)>0.0) tepv1=(tepv[p]*temvs[p]+tepv[q]*temvs[q])/tem1;
	p=3; q=5; if (fabs(tem2)>0.0) tepv2=(tepv[p]*temvs[p]+tepv[q]*temvs[q])/tem2;
	tem1=tem1/2e0; tem2=tem2/2e0;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tepv2-tepv1)/(tem2-tem1); 
	else k1=0.0; k2=tepv2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k] = k1; return isdan;
}
double *danPoTemTepl2073(double *temvs, double *tepv, int n) //Фракция 2-0,7 мм, после обжига при 1000 °С
{
	double tepv1=0.0, tepv2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n]; 
	int k=0, p=0, q=0;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=6; q=8; tem1=temvs[p]+temvs[q];
	p=7; q=9; tem2=temvs[p]+temvs[q];
	p=6; q=8; tepv1=(tepv[p]*temvs[p]+tepv[q]*temvs[q])/tem1;
	p=7; q=9; tepv2=(tepv[p]*temvs[p]+tepv[q]*temvs[q])/tem2;
	tem1=tem1/2e0; tem2=tem2/2e0;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tepv2-tepv1)/(tem2-tem1); 
	else k1=0.0; k2=tepv2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl2074(double *temvs, double *tepv, int n) //Фракция 2-0,7 мм, после повторного обжига при 1000 °С
{
	double tepv1=0.0, tepv2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n], hf=1e0; 
	int k=0, p=0;
	if (!isdan) { cout << snmv << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=2; tem1=temvs[p];
	p=3; tem2=temvs[p];
	p=2; tepv1=(tepv[p]*temvs[p])/tem1;
	p=3; tepv2=(tepv[p]*temvs[p])/tem2;
	tem1=tem1/hf; tem2=tem2/hf;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tepv2-tepv1)/(tem2-tem1); 
	else k1=0.0; k2=tepv2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
//---------
double **arrTemCold84(double **mu)
{
	int k=0, n=26; double *temcol=new double[n], *po=new double[1], nf=0.0, hf=1e0;
	if ((!temcol) || (!po)) { cout << snmv << endl; k = getchar(); exit(1); }
	temcol[k] = 168.0;  k++; temcol[k] = 369.0;  k++; 
	temcol[k] = 168.0;  k++; temcol[k] = 369.0;  k++;
	temcol[k] = 148.6;  k++; temcol[k] = 356.5;  k++; 
	temcol[k] = 184.0;  k++; temcol[k] = 396.0;  k++;
	temcol[k] = 148.75; k++; temcol[k] = 350.0;  k++; 
	temcol[k] = 166.0;  k++; temcol[k] = 375.0;  k++;
	temcol[k] = 171.0;  k++; temcol[k] = 383.5;  k++; 
	temcol[k] = 106.75; k++; temcol[k] = 242.0;  k++;
	temcol[k] = 123.0;  k++; temcol[k] = 294.0;  k++; 
	temcol[k] = 111.0;  k++; temcol[k] = 240.0;  k++;
	temcol[k] = 109.0;  k++; temcol[k] = 232.75; k++; 
	temcol[k] = 127.0;  k++; temcol[k] = 291.0;  k++;
	temcol[k] = 221.0;  k++; temcol[k] = 443.75;
	for (k=0; k<n; k++) temcol[k]=temcol[k] + tev0;
	for (k=0; k<n; k++) nf=nf+hf; k=0; po[k]=nf;
	k=0; mu[k]=po; k++; mu[k]=temcol;
	return mu;
}
double *arrTemHigh84()
{
	int k=0, n=26; double *temh=new double[n];
	if (!temh) { cout << snmv << endl; k=getchar(); exit(1); }
	temh[k] = 545.0; k++; temh[k] = 927.0;  k++; 
	temh[k] = 530.0; k++; temh[k] = 925.0;  k++;
	temh[k] = 560.0; k++; temh[k] = 950.0;  k++; 
	temh[k] = 560.0; k++; temh[k] = 900.0;  k++;
	temh[k] = 558.0; k++; temh[k] = 950.0;  k++; 
	temh[k] = 540.0; k++; temh[k] = 920.0;  k++;
	temh[k] = 540.0; k++; temh[k] = 920.0;  k++; 
	temh[k] = 600.0; k++; temh[k] = 1000.0; k++;
	temh[k] = 580.0; k++; temh[k] = 1000.0; k++; 
	temh[k] = 587.0; k++; temh[k] = 1000.0; k++;
	temh[k] = 590.0; k++; temh[k] = 1000.0; k++; 
	temh[k] = 580.0; k++; temh[k] = 1000.0; k++;
	temh[k] = 483.0; k++; temh[k] = 850.0;
	for (k=0; k<n; k++) temh[k]=temh[k]+tev0;
	return temh;
}
double *arrTepPot84()
{
	int k=0, n=26; double *tep=new double[n];
	if (!tep) { cout << snmv << endl; k=getchar(); exit(1); }
	tep[k] = 3.9805;  k++; tep[k] = 9.5532;   k++; 
	tep[k] = 3.6447;  k++; tep[k] = 9.4779;   k++;
	tep[k] = 3.55732; k++; tep[k] = 11.4997;  k++; 
	tep[k] = 4.4624;  k++; tep[k] = 11.6021;  k++;
	tep[k] = 4.7766;  k++; tep[k] = 11.5016;  k++; 
	tep[k] = 3.4023;  k++; tep[k] = 10.2068;  k++;
	tep[k] = 3.92812; k++; tep[k] = 11.17333; k++; 
	tep[k] = 3.5144;  k++; tep[k] = 8.6593;   k++;
	tep[k] = 2.977;   k++; tep[k] = 8.0448;   k++; 
	tep[k] = 3.352;   k++; tep[k] = 9.218;    k++;
	tep[k] = 3.0313;  k++; tep[k] = 7.7946;   k++; 
	tep[k] = 3.0671;  k++; tep[k] = 6.1342;   k++;
	tep[k] = 1.73466; k++; tep[k] = 4.32967;
	return tep;
}
double **arrTemCold207(double **mu)
{
	int k=1, n=10; double *temcol=new double[n], *po=new double[k], nf=0.0, hf=1e0;
	if ((!temcol) || (!po)) { cout << snmv << endl; k=getchar(); exit(1); } k=0;
	temcol[k] = 109.0;  k++; temcol[k] = 235.0; k++; 
	temcol[k] = 101.0;  k++; temcol[k] = 199.0; k++;
	temcol[k] = 108.75; k++; temcol[k] = 238.0; k++; 
	temcol[k] = 124.0;  k++; temcol[k] = 266.0; k++;
	temcol[k] = 111.0;  k++; temcol[k] = 262.0;
	for (k=0; k<n; k++) temcol[k]=temcol[k]+tev0;
	for (k=0; k<n; k++) nf=nf+hf; k=0; po[k]=nf;
	k=0; mu[k]=po; k++; mu[k]=temcol;
	return mu;
}
double *arrTemHigh207()
{
	int k=0, n=10; double *temh=new double[n];
	if (!temh) { cout << snmv << endl; k=getchar(); exit(1); }
	temh[k] = 585.0; k++; temh[k] = 1000.0; k++; 
	temh[k] = 603.0; k++; temh[k] = 1000.0; k++;
	temh[k] = 603.0; k++; temh[k] = 1000.0; k++; 
	temh[k] = 571.5; k++; temh[k] = 1000.75; k++;
	temh[k] = 583.0; k++; temh[k] = 1000.0;
	for (k=0; k<n; k++) temh[k]=temh[k]+tev0;
	return temh;
}
double *arrTepPot207()
{
	int k=0, n=10; double *tepot=new double[n];
	if (!tepot) { cout << snmv << endl; k = getchar(); exit(1); }
	tepot[k] = 3.9596;  k++; tepot[k] = 8.6377; k++; 
	tepot[k] = 2.3003;  k++; tepot[k] = 5.3674; k++;
	tepot[k] = 3.56149; k++; tepot[k] = 7.123;  k++; 
	tepot[k] = 2.12992; k++; tepot[k] = 7.6956; k++;
	tepot[k] = 2.3003;  k++; tepot[k] = 6.9009;
	return tepot;
}
double **vydelPol(int v, int n, double **mu, double **muv, int f, int fl)
{ 
	int q=0, qn=n, k=0, m=0, w=f-1, j=1, nk=n, *mui=NULL, x=0, qk=0;
	double *po=NULL, nf=0.0, t=0.0, hf=1e0, *vm=NULL, e=1e-6, r=0.0, s=0.0;
	k=w; po=mu[k]; k=0; nf=po[k]; 
	t=e; while (t<nf) { k++; t=t+hf; } nk=k; mui=new int[nk]; if (!mui) { cout << snmv << endl; k=getchar(); exit(1); }
	q=0; for (k=0; k<nk; k++) {
		x=1; for (m=0; m<w; m++) { 
			po=mu[m]; if ((po[k]<e) && (x>0)) { x=-1; break; } } 
		if (x>0) { mui[k]=k; q++; } else mui[k]=-1; } qn=q;
		if (fl>0) { m=1; po=mu[m]; for (k=0; k<nk; k++) if (po[k]>templa) { mui[k]=-1; qn--; } } //if (fl>0) { cout << "qn = " << qn << "\tnf = " << nk << "\tw = " << w << "\tn = " << n << "\t"; } 
	for (m=0; m<w; m++) {
	vm=new double[qn]; if (!vm) { cout << snmv << endl; k=getchar(); exit(1); } //for (k=0; k<qn; k++) vm[k]=0.0;
	q=0; po=mu[m]; for (k=0; k<nk; k++) {
		x=mui[k]; if (x>=0) { vm[q]=po[x]; q++; } } muv[m]=vm; } qn=q; qk=q; if (mui) delete[]mui; 
	nf=0.0; for (k=0; k<qn; k++) nf=nf+hf; 
	k=1; po=new double[k]; if (!po) { cout << snmv << endl; k=getchar(); exit(1); } k=0; po[k]=nf; muv[w]=po; //if (fl>0) { q=w; po=muv[q]; q=0; t=po[q]; cout << "r = " << t << "\t"; s=e; q=0; while (s<t) { s=s+hf; q++; } qk=q; q=0; po=muv[q]; for (q=0; q<qk; q++) cout << "temc = " << po[q] << "\t"; cout << endl; q=1; po=muv[q]; for (q=0; q<qk; q++) cout << "temh = " << po[q] << "\t"; cout << endl; q=2; po=muv[q]; for (q=0; q<qk; q++) cout << "qo = " << po[q] << "\t"; cout << endl; q=3; po=muv[q]; for (q=0; q<qk; q++) cout << "ktp = " << po[q] << "\t"; cout << endl; q=4; po=muv[q]; for (q=0; q<qk; q++) cout << "tems = " << po[q] << "\t"; cout << endl; q=w; po=mu[q]; q=0; t=po[q]; cout << "r = " << t << endl; s=e; q=0; while (s<t) { s=s+hf; q++; } qk=q; q=0; po=mu[q]; for (q=0; q<qk; q++) cout << "temc = " << po[q] << "\t"; cout << endl; q=1; po=mu[q]; for (q=0; q<qk; q++) cout << "temh = " << po[q] << "\t"; cout << endl; q=2; po=mu[q]; for (q=0; q<qk; q++) cout << "qo = " << po[q] << "\t"; cout << endl; q=3; po=mu[q]; for (q=0; q<qk; q++) cout << "ktp = " << po[q] << "\t"; cout << endl; q=4; po=mu[q]; for (q=0; q<qk; q++) cout << "tems = " << po[q] << "\t"; cout << endl;  cout << mu << "\t" << muv << "\t"; }
	for (k=0; k<f; k++) { vm=mu[k]; //if ((fl>0) && (vm) && (k<w)) { for (q=0; q<nk; q++) cout << "k = "<< k << "\tq = " << q << "\tp = " << vm[q] << "\t"; cout << endl; } 
	if ((fl<0) && (vm)) delete[]vm; } //if (fl>0) { cout << "qn = " << qn << "\tn = " << nf << "\t"; } 
	return muv;
}
double **gotSredKTPTK()
{
int f=3, k=5, j=1; 
double **mu=new double*[f], *te=new double[k], *ktptk=new double[k], nf=0.0, hf=1e0, *po=new double[j], dt=1e2;
if ((!ktptk) || (!te) || (!mu) || (!po)) { cout << snmv << endl; k=getchar(); exit(1); }
j=0; ktptk[j]=1.0; j++; ktptk[j]=1.32; j++; ktptk[j]=1.66; j++; ktptk[j]=2.02; j++; ktptk[j]=2.39;
j=0; te[j]=573.15; for (j=1; j<k; j++) te[j]=te[j-1]+dt;
nf=0.0; for (j=0; j<k; j++) nf=nf+hf; j=0; po[j]=nf;
k=j; mu[k]=ktptk; k++; mu[k]=te; k++; mu[k]=po;
return mu;
}