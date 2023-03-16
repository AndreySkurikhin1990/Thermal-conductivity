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
#define vyfv 0 //выбор фракции: 0 - фракция 2-0,7 мм, 1 - фракция 8-4 мм, 2 - фракция 1,6-0,35 мм
#define vyuv 1 //выбор укладки: 1 - плоскопараллельная, 2 - вертикальная
#define vysv 1 //выбор состояния: 0 - исходное, 1 - после повторных измерений, 2 - после прокаливания при 1000 град С
using namespace std;
const int dsv = 50, dmkvover = 28, dmkooscv = 14, N = 13, ks = 10, dmkov = 4, qtn = 8, nxtv = 30, nnxtv = 0;
const int vtvn = 0, isrp=0, vpkf=0, vpmf=0; //0 - старые, 1 - новые
const double detev = 1e2, tev0 = 273.15, nnxtfv = 0.0, ssiv84 = 467.0*1e-3, salv84 = 129.0*1e-3, smgv84 = 282.0*1e-3;
const double ssiv207 = 458.0*1e-3, salv207 = 14.0*1e-2, smgv207 = 29.0*1e-2;
const double tnoscv = 3e2, tnacv = 2e2; //dkoscv - максимальное отклонение от экспериментальных данных в литературных источниках
const double dtoscv = 1e2, tocrasver = 1e-4, por207 = 66.35*1e-2, poris84 = 55.75*1e-2, porin84=81.53*1e-2, templa = 134.0*1e1 + tev0;
const double epsi = 1e-15, y0ver = 3e1*1e-3, pksvv = 0.0, dkoalv = 0.444905, poro84=86.61*1e-2; //dkoalv - определение КП из КО
const double poro16035=84.36*1e-2, por16035=83.97*1e-2;
const char nfv = 'A';
struct derevo {
	int otre; //1 - отражение или 0 - пропускание
	int ste; //номер стенки
	int vis; //видимость: 1 - виден, 0 - нет
	int lev; //номер уровня
	struct derevo *back; //указатель назад
	struct derevo *next; //указатель вперед
};
class KoePog { public: double alp, tem; KoePog *nex; };
//----------------------------------------
double *alphaVer = NULL, *Tver = NULL, *Aver = NULL, *Rver = NULL, hk = 0.0, *Raver = NULL, *Taver = NULL, *Aaver = NULL, *Rtver = NULL, *Ttver = NULL, *Atver = NULL;
double *mgover = NULL, *siover = NULL, *alover = NULL; 
int vtv = 0, sctxv = 0, cel = 4, dkoscvl = 6, cemv = 11, vtvk=cemv; /*выбор температуры вермикулита*/
char *snv = NULL, *ssv = NULL, *skptv = NULL, *svsv = NULL, *snmv = NULL, *sfnov = NULL;
char *sfatv = NULL, *sfov = NULL, *svfdv = NULL, *svfdvu = NULL;
double *qobv = NULL, *etev = NULL, *txver = NULL, *qxver = NULL, *lxver = NULL, *kxver = NULL;
double *dtxver = NULL, *stchsrver = NULL, *tgorv = NULL, *tholv = NULL;
double *ktrv = NULL, *dkoscvm = NULL, *dkoscvt = NULL, hver = 0.0, hvove = 0.0, hkver = 0.0, qobver = 0.0, *ktpvover;
double *vtever, porver = 0.0, *ecktpver, *tsredver, dkospv = 1.3002913762; //dkospv - дополнительный коэффициент для КО с учетом пористой структуры;
double *kttkv = NULL, cov = 0.0, bov = 0.0, *mkov = NULL, *ektpv = NULL, *kektpv = NULL, *mtsv = NULL, tmav, tmiv, dpctv = 1e0;
double *Tpctv, tmiver, tmaver, tnav = tnacv + tev0, *tkuscv, *kuscv, cnv = 1e-2, ckv = 1e2, bnv = 1e-2;
double bkv = 1e2, *temrasv = NULL, tnrv, *ooxver = NULL, *dpctvm, nxtfv, ssiv = 368.0*1e-3, salv = 132.0*1e-3, smgv = 217.0*1e-3;
//----------------------------------------
void zadrktVerNac();
double *podschchieleSha(int, int, int, double *, double *);
double *izstNVer(int, int, int, int);
double *kopoVer(double, double, int, double, int, int, double, double, int, double *);
double opredLuchSostVer(double *, double *, int, int, int);
double *koefPribSha(double *, double *, int, double *);
double *zadrktVer(int, int, double, int, double, double, int, int, double, int, int, double *);
double **chaRTAVer(int, double *, double *, double *, double *, double *, double *, int);
double **izmRTAVer(double *, int, int, double *, double *, double *, double *, double *, double *, int);
double *opredKoefOtr(double *, double *, int, double, int, int, double *, double *, int);
void napstrver();
double *oprkoefKTPiskhchao(int, int, double *, double *, double, double *, double *, double *, int);
double opredKTPTKTochSha(double *, double *, double, int);
double *opredKTPTverKarkVerm(double *, double *, double, double, double, double, int, int, double, double, int, double *, double *, double *, int, int, int);
double *reshnewtrafs(double *, double *, double *, double *, int, double *, double *, int, double *, double *, int, double *, double *, int, double, double *, double *, int, double, double, double, double *, int, double, double, double);
void initarrver(int, double, double, double);
void NapMasVozdSha(double*, double *, int);
void vyvodfile(double *, int, int, double, char *);
void novNapMas(double);
void zapisvfile(double *, int, char *);
double BolTochRasAlpha(int, int, double, double, double, double *, double *, double *, int, int, char *, double *, double *, int, double *, double *, int, double *, double *, double *, double *, double *, double *);
void osvpamver();
double **RaschRTAVer(int, double, double, double, int, int, double, int, double, int, int, int);
double epsisredver(double, double *, double *, int, double *, double *, int);
double *KoefPoglRosselNac(double *, int, int, double, double, double, double *, double *, int, double, double, double *, double *, int, int, int);
double **opredTempLPStenSha(double *, double *, double *, double *, double, int, double *, double *, double *, double *, int, int, double, double, double, int, char *);
double **RaschVneshIzluchSha(double *, double *, double *, double *, double *, int, double, double, char *);
double **RasLuchPloTepPot(int, double *, double *, double *, double *, double *, double *, int, double *, double *, double *, char *);
double **RaschSobLuchPlotTepPot(int, double *, double *, double *, double *, double *, double *, double, double, double *, int, char *);
double *opredTempStenShaFragm(double *, int, double *, double *, double *, double *, int, int, double, double, double, double, double);
double LuchKTPChudnovsky(double *, double, int, double);
double *rasPorpoRazVer(double, int, int, int, int, int);
double RaschAlphaTvKarVer();
double *EffectTols(double *, double *, double *, double, double, int);
double RasFracXeffVer60(int);
double RasFracXeffVer60_100(int);
double RasFracXeffVer100_150(int);
double RasFracXeffVer150_200(int);
double *koefoslab(double, double, double, double *, int, double *);
double ReflSredVer(double);
double reflver(double);
double oprProcSoderpoPoris(double *, double *, double, int);
double SeryeStenkiRasIzl(double, double, double, double *, double *, double *, double *, double *, int *);
double F0_lamT(double);
void opredtemphcVer(double *, double *, double *, double *, double *, int, int, double);
double *oprRasTemNachVer(int, int, double *, double *, int, double *, double, double, double, double, double, int);
double *oprEffDoliTepPerenVer(double, double, double);
double *KorrZnachVozdProsVer(double, double, double, int);
double *FuncRaschIzl(double, double, double, double, int, double *, double, double *, double, double, double, double, int, double *, double *, int, int, char *, double, int, double *, double, double *, double *, int, double *, double *, int, double *, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, int, double, double, int, double *, double *, double, double, double, char *, int, double);
double KorrZnachVozdProsVermik(double, double, double);
double *vydelPol(double *, double *, double *, double *, double *, int, int);
//----------------
double *napMasEKTPVer(int, int, int, double *, int, double, int, int, int);
double *arrTem_VVF2(int); //массив температур - экспериментальные данные на Netzsch - хаотичная засыпка фракции 2-0,7 мм (исходный)
double *arrKTP_VVF2(); //массив КТП вермикулита
double *arrKTP_VVF1();
double *arrTem1VVF1();
double *arrTem2VVF1(int);
double *arrTem3VVF1();
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
double *arrTemCold84(int);
double *arrTemHigh84();
double *arrTepPot84();
double *arrTemCold207(int);
double *arrTemHigh207();
double *arrTepPot207();
double *danIskh207(double *, int, int);
double *PoiskZavVelTemVer(int, double*);
void napMasEKTPVerNac(double, double, double);
double **RaschRTAitom(int, double, double, double, int, int, double, int, double, int, int, int);
double bbfn(double);
//----------------
void napstrver()
{
	if ((!sfatv) || (!sfov) || (!snv) || (!ssv) || (!skptv) || (!svsv) || (!snmv) || (!sfnov)) { cout << "No_memory!" << endl; getchar(); exit(1); }
	int k = 0; nxtfv = 0.0; while (k < nxtv) { nxtfv = nxtfv + 1e0; k++; } //cout << "nxtfv = " << nxtfv << endl;
	k=0;
	snv[k] = 'D'; k++; snv[k] = ':';  k++; snv[k] = '\\'; k++; snv[k] = '\\'; k++; snv[k] = '_'; k++;
	snv[k] = 'А'; k++; snv[k] = 'с';  k++; snv[k] = 'п';  k++; snv[k] = 'и';  k++; snv[k] = 'р'; k++;
	snv[k] = 'а'; k++; snv[k] = 'н';  k++; snv[k] = 'т';  k++; snv[k] = 'у';  k++; snv[k] = 'р'; k++;
	snv[k] = 'а'; k++; snv[k] = '\\'; k++; snv[k] = '\\'; k++; snv[k] = 't';  k++; snv[k] = 'm'; k++;
	snv[k] = 'p'; k++; snv[k] = '\\'; k++; snv[k] = '\\'; k++; snv[k] = '\0';
	k = 0;
	ssv[k] = 'C';  k++; ssv[k] = ':'; k++; ssv[k] = '\\'; k++; ssv[k] = '\\'; k++; ssv[k] = 'U';  k++;
	ssv[k] = 's';  k++; ssv[k] = 'e'; k++; ssv[k] = 'r';  k++; ssv[k] = 's';  k++; ssv[k] = '\\'; k++;
	ssv[k] = '\\'; k++; ssv[k] = 'А'; k++; ssv[k] = 'н';  k++; ssv[k] = 'д';  k++; ssv[k] = 'р';  k++;
	ssv[k] = 'е';  k++; ssv[k] = 'й'; k++; ssv[k] = '\\'; k++; ssv[k] = '\\'; k++; ssv[k] = 'D';  k++;
	ssv[k] = 'o';  k++; ssv[k] = 'c'; k++; ssv[k] = 'u';  k++; ssv[k] = 'm';  k++; ssv[k] = 'e';  k++;
	ssv[k] = 'n';  k++; ssv[k] = 't'; k++; ssv[k] = 's';  k++; ssv[k] = '\\'; k++; ssv[k] = '\\'; k++;
	ssv[k] = '_';  k++; ssv[k] = 'А'; k++; ssv[k] = 'с';  k++; ssv[k] = 'п';  k++; ssv[k] = 'и';  k++;
	ssv[k] = 'р';  k++; ssv[k] = 'а'; k++; ssv[k] = 'н';  k++; ssv[k] = 'т';  k++; ssv[k] = 'у';  k++;
	ssv[k] = 'р';  k++; ssv[k] = 'а'; k++; ssv[k] = '\\'; k++; ssv[k] = '\\'; k++; ssv[k] = 't';  k++;
	ssv[k] = 'm';  k++; ssv[k] = 'p'; k++; ssv[k] = '\\'; k++; ssv[k] = '\\'; k++; ssv[k] = '\0';
	k = 0;
	skptv[k] = 'K'; k++; skptv[k] = 'o'; k++; skptv[k] = 'e'; k++; skptv[k] = 'f'; k++; skptv[k] = 'f'; k++;
	skptv[k] = 'i'; k++; skptv[k] = 'c'; k++; skptv[k] = 'i'; k++; skptv[k] = 'e'; k++; skptv[k] = 'n'; k++;
	skptv[k] = 't'; k++; skptv[k] = '_'; k++; skptv[k] = 'p'; k++; skptv[k] = 'o'; k++; skptv[k] = 'g'; k++;
	skptv[k] = 'l'; k++; skptv[k] = 'o'; k++; skptv[k] = 's'; k++; skptv[k] = 'c'; k++; skptv[k] = 'h'; k++;
	skptv[k] = 'e'; k++; skptv[k] = 'n'; k++; skptv[k] = 'i'; k++; skptv[k] = 'y'; k++; skptv[k] = 'a'; k++;
	skptv[k] = '_'; k++; skptv[k] = 'v'; k++; skptv[k] = 'e'; k++; skptv[k] = 'r'; k++; skptv[k] = '_'; k++;
	skptv[k] = 'T'; k++; skptv[k] = '.'; k++; skptv[k] = 't'; k++; skptv[k] = 'x'; k++; skptv[k] = 't'; k++;
	skptv[k] = '\0';
	k = 0;
	svsv[k] = 'D'; k++; svsv[k] = 'o'; k++; svsv[k] = 'l'; k++; svsv[k] = 'i'; k++; svsv[k] = '_'; k++;
	svsv[k] = 'p'; k++; svsv[k] = 'r'; k++; svsv[k] = 'o'; k++; svsv[k] = 'p'; k++; svsv[k] = '_'; k++;
	svsv[k] = 'V'; k++; svsv[k] = 'e'; k++; svsv[k] = 'r'; k++; svsv[k] = 'm'; k++; svsv[k] = 'i'; k++;
	svsv[k] = 'k'; k++; svsv[k] = '-'; k++; svsv[k] = nfv; k++; svsv[k] = '.'; k++; svsv[k] = 't'; k++;
	svsv[k] = 'x'; k++; svsv[k] = 't'; k++; svsv[k] = '\0';
	k = 0;
	snmv[k] = 'N'; k++; snmv[k] = 'o'; k++; snmv[k] = '_'; k++; snmv[k] = 'm'; k++; snmv[k] = 'e'; k++;
	snmv[k] = 'm'; k++; snmv[k] = 'o'; k++; snmv[k] = 'r'; k++; snmv[k] = 'y'; k++; snmv[k] = '!'; k++; snmv[k] = '\0';
	k = 0;
	sfnov[k] = 'F'; k++; sfnov[k] = 'i'; k++; sfnov[k] = 'l'; k++; sfnov[k] = 'e'; k++; sfnov[k] = '_'; k++;
	sfnov[k] = 'i'; k++; sfnov[k] = 's'; k++; sfnov[k] = '_'; k++; sfnov[k] = 'n'; k++; sfnov[k] = 'o'; k++;
	sfnov[k] = 't'; k++; sfnov[k] = '_'; k++; sfnov[k] = 'o'; k++; sfnov[k] = 'p'; k++; sfnov[k] = 'e'; k++;
	sfnov[k] = 'n'; k++; sfnov[k] = '!'; k++; sfnov[k] = '\0';
	k = 0;
	svfdvu[k] = 'V'; k++; svfdvu[k] = 'y'; k++; svfdvu[k] = 'v'; k++; svfdvu[k] = 'o'; k++; svfdvu[k] = 'd'; k++;
	svfdvu[k] = 'v'; k++; svfdvu[k] = 'F'; k++; svfdvu[k] = 'i'; k++; svfdvu[k] = 'l'; k++; svfdvu[k] = 'e'; k++;
	svfdvu[k] = '.'; k++; svfdvu[k] = 't'; k++; svfdvu[k] = 'x'; k++; svfdvu[k] = 't'; k++; svfdvu[k] = '\0';
	for (k = 0; k < (2 * dsv); k++) { sfatv[k] = '\0'; sfov[k] = '\0'; }
	strcpy(sfatv, ssv); strcat(sfatv, skptv); sfatv[strlen(sfatv) + 1] = '\0'; 
	strcpy(sfov, ssv); strcat(sfov, svsv); sfov[strlen(sfov) + 1] = '\0'; 
	strcpy(svfdv, ssv); strcat(svfdv, svfdvu); svfdv[strlen(svfdv) + 1] = '\0'; 
}
void initarrver(int koel, double wmg, double wsi, double wal)
{
	if ((!vyfv) || (vyfv==2)) { ssiv = ssiv207; salv = salv207; smgv = smgv207; }
	else if (vyfv==1) { ssiv = ssiv84; salv = salv84; smgv = smgv84; }
	snv = new char[dsv]; ssv = new char[dsv]; skptv = new char[dsv]; svsv = new char[dsv];
	snmv = new char[dsv]; sfnov = new char[dsv]; svfdvu = new char[dsv];
	int j; if ((!snv) || (!ssv) || (!skptv) || (!svsv) || (!snmv) || (!sfnov) || (!svfdvu)) 
	{ cout << "No  memory!" << endl; j = getchar(); exit(1); }
	for (j = 0; j < dsv; j++) { snv[j] = '\0'; ssv[j] = '\0'; skptv[j] = '\0'; 
	svsv[j] = '\0'; snmv[j] = '\0'; sfnov[j] = '\0'; svfdvu[j] = '\0'; }
	j = 2 * dsv; sfatv = new char[j]; sfov = new char[j]; svfdv = new char[j];
	if ((!sfatv) || (!sfov) || (!svfdv)) { cout << "No memory!" << endl; j = getchar(); exit(1); }
	for (j = 0; j < (2 * dsv); j++) { sfatv[j] = '\0'; sfov[j] = '\0'; svfdv[j] = '\0'; }
	napstrver();
	dkoscvm = new double[dkoscvl]; dkoscvt = new double[dkoscvl]; int k=0;
	if ((!dkoscvm) || (!dkoscvt)) { cout << snmv << endl; k = getchar(); exit(1); }
	double tnd = 6e2, dtd = 2e2, tm; k=0; dkoscvt[k] = tnd; for (k = 1; k < dkoscvl; k++) dkoscvt[k] = dkoscvt[k - 1] + dtd;
	k = 0; dkoscvm[k] = 4.68; k++; dkoscvm[k] = 4.69; k++; dkoscvm[k] = 5.65; k++; dkoscvm[k] = 13.17; k++; dkoscvm[k] = 20.2; k++; dkoscvm[k] = 27.81;
	for (k = 0; k < dkoscvl; k++) {
		tm = dkoscvm[k] / 1e2; dkoscvm[k] = 1e0 - tm; //cout << "tem = " << dkoscvt[k] << "\tdkosc = " << dkoscvm[k] << endl; 
	} 
	qobv = new double[cemv]; etev = new double[cemv]; mkov = new double[cemv]; ektpv = new double[cemv];
	tgorv = new double[cemv]; tholv = new double[cemv]; tsredver = new double[cemv]; stchsrver = new double[cemv];
	Tver = new double[ks]; Aver = new double[ks]; Rver = new double[ks]; Raver = new double[ks];
	Taver = new double[ks]; Aaver = new double[ks]; Rtver = new double[ks]; Ttver = new double[ks];
	Atver = new double[ks]; alphaVer = new double[ks]; Tpctv = new double[ks];
	if ((!Tver) || (!Aver) || (!Rver) || (!Raver) || (!Taver) || (!Aaver) || (!Ttver) || (!Rtver) || (!Atver) || (!alphaVer) || (!Tpctv)) { cout << snmv << endl; j = getchar(); exit(1); }
	if ((!qobv) || (!etev) || (!mkov) || (!ektpv) || (!tgorv) || (!tholv) || (!stchsrver) || (!tsredver)) { cout << snmv << endl; j = getchar(); exit(1); }
	for (j = 0; j < cemv; j++) {
		qobv[j] = 0.0; etev[j] = 0.0; mkov[j] = 0.0;
		ektpv[j] = 0.0; stchsrver[j] = 0.0; tgorv[j] = 0.0; tholv[j] = 0.0; tsredver[j] = 0.0;
	}
	for (j = 0; j < ks; j++) {
		Aver[j] = 0.0; Rver[j] = 0.0; Raver[j] = 0.0; Taver[j] = 0.0;
		Aaver[j] = 0.0; Rtver[j] = 0.0; Ttver[j] = 0.0; Atver[j] = 0.0; Tver[j] = 0.0;
	}
	txver = new double[koel]; qxver = new double[koel]; lxver = new double[koel]; kxver = new double[koel]; dtxver = new double[koel]; ooxver = new double[koel];
	if ((!txver) || (!qxver) || (!lxver) || (!kxver) || (!dtxver) || (!ooxver)) { cout << snmv << endl; j = getchar(); exit(1); }
	for (j = 0; j < koel; j++) { txver[j] = 0.0; qxver[j] = 0.0; lxver[j] = 0.0; kxver[j] = 0.0; dtxver[j] = 0.0; ooxver[j] = 0.0; }
	tkuscv = new double[dmkooscv]; kuscv = new double[dmkooscv];
	if ((!tkuscv) || (!kuscv)) { cout << snmv << endl; k = getchar(); exit(1); }
	k=0; tkuscv[k] = tnoscv; for (k = 1; k < dmkooscv; k++) tkuscv[k] = tkuscv[k - 1] + dtoscv;
	kuscv = koefoslab(wmg, wsi, wal, tkuscv, dmkooscv, kuscv); //for (k=0; k<dmkooscv; k++) cout << "te = " << tkuscv[k] << "\tku = " << kuscv[k] << endl;
	ktpvover = new double[dmkvover]; vtever = new double[dmkvover];
	for (j = 0; j < dmkvover; j++) { ktpvover[j] = 0.0; vtever[j] = 0.0; }
	novNapMas(tnav); napMasEKTPVerNac(wmg, wsi, wal); for (k=0; k<cemv; k++) cout << "te = " << etev[k] << "\tktp = " << ektpv[k] << endl; 
}
void novNapMas(double tnach)
{
	NapMasVozdSha(ktpvover, vtever, dmkvover);
	int j, k;
	if ((!vysv) || (vysv==1)) {  //исходный или после повторных измерений
	if (!vyfv) { if (!vpmf) porver = por207; else if (vpmf==1) porver=por16035; }
	else if (vyfv==1) { if (!vpkf) porver = poris84; else if (vpkf==1) porver=porin84; }
	else if (vyfv==2) porver = poro16035; }
	if (vysv==2) { //после обжига
	if ((!vyfv) || (vyfv==2)) porver = poro16035;
	else if (vyfv==1) porver = poro84; }
	etev[0] = tnach; 
	for (j = 1; j < cemv; j++) etev[j] = etev[j - 1] + detev; //cout << "por = " << porver << endl;
}
void osvpamver()
{
	delete[]kektpv; delete[]kttkv; delete[]mkov; delete[]qobv; delete[]snv; delete[]ssv;
	delete[]skptv; delete[]svsv; delete[]snmv; delete[]sfnov; delete[]svfdvu; delete[]sfatv; delete[]sfov;
	delete[]svfdv; delete[]etev; delete[]ektpv; delete[]txver; delete[]qxver; delete[]lxver; delete[]kxver; delete[]ooxver;
	delete[]dtxver; delete[]ktpvover; delete[]vtever; delete[]Tpctv; delete[]alphaVer; delete[]tgorv; delete[]tholv;
	delete[]Tver; delete[]Aver; delete[]Rver; delete[]Raver; delete[]Taver; delete[]Aaver; delete[]dkoscvm; delete[]dkoscvt;
	delete[]Rtver; delete[]Ttver; delete[]Atver; delete[]stchsrver; delete[]tsredver;
}
void zadrktVerNac()
{
	int j, jk = nxtv, jn = nnxtv, k, q, f = 6; 
	initarrver(jk, smgv, ssiv - pksvv, salv + pksvv); double hf = 1e0, nf = 0.0, ka, kb;
	double dhk = y0ver / fabs(nxtfv - nnxtfv), hnver = nnxtfv*dhk, hvko = (13e1)*(1e-6);
	double hvh = 1e-6, hvna, p, r, d = 0.0, *atr = NULL, t, ko, *po, **mu; 
	for (j = 0; j < ks; j++) d = d + hf; //cout << "d = " << d << "\t";
	dkospv = RaschAlphaTvKarVer(); 
	vtvk=cemv;
	for (vtv = vtvn; vtv<vtvk; vtv++) { //пробегаем по температуре
		po = rasPorpoRazVer(porver, 0, 1, vysv, isrp, vpkf); ko = po[0]; delete[]po; cout << "Sred razm por = " << ko << "\t";
		r = d; t = ko; po = KorrZnachVozdProsVer(ko, r, porver, 0); ko=po[0]; cout << "Korr sred razm por = " << ko << "\t";
		r = d; ko = t; po = KorrZnachVozdProsVer(ko, r, porver, 1); p=po[vtv]; cout << "Dol Plo = " << p << "\t";
		r = d; ko = t; ko = KorrZnachVozdProsVermik(ko, r, porver); hvna = ko; hvko = ko;
		if (vtv>vtvn) {
			q = jk - jn; ka = (tholv[vtv] - tgorv[vtv]) / y0ver; kb = tgorv[vtv]; hkver = hnver;
			for (k = 0; k < q; k++) { temrasv[k] = kb + ka*hkver; hkver = hkver + dhk; }
		}
		hvove = hvna; hver = hvove*(1e0 - porver) / porver;
		while (hvove <= hvko) {
			j = jn; hkver = hnver*dhk; //пробегаем по размерам пор
			while ((j < jk) && (hkver < y0ver)) {
				cout << "hk = " << hkver << endl; sctxv = j; //пробегаем по координате
				mu = RaschRTAVer(ks, hver, 0.0, 0.0, 1, vtv, hvove, 0, etev[vtv], 0, 0, 1);
				Raver = mu[0]; Taver = mu[1]; Aaver = mu[2]; Rver = mu[3]; Tver = mu[4]; Aver = mu[5]; Rtver = mu[6]; Ttver = mu[7]; Atver = mu[8]; delete[]mu;
				atr = zadrktVer(j, ks, d, vtv, hver, hvove, 1, 0, etev[vtv], 0, 0, atr);
				kxver[j - jn] = hkver; hkver = hkver + dhk; j++; delete[]atr;
			}
			vyvodfile(lxver, jk - jn, 0, hvove, sfov);
			for (j = 0; j<jk - jn - 1; j++) {
				p = kxver[j + 1] - kxver[j];
				r = txver[j] - txver[j + 1]; if (fabs(p)>0.0) r = fabs(r / p); else r = 0.0;
				p = fabs(qxver[j + 1] + qxver[j]) / 2e0; lxver[j] = p / r;
			}
			vyvodfile(ooxver, jk - jn, 2, hvove, sfov);
			vyvodfile(txver, jk - jn, 2, hvove, sfov);
			vyvodfile(lxver, jk - jn - 1, 2, hvove, sfov);
			hvove = hvove + hvh;
		}
	}
}
double opredLuchSostVer(double *prot, double *Tsr, int le, int zf, int vyte)
{
	int k = 0, j = 0, q = 0, m = 0, l = 3, mma = 3000, mm, mt = 0, vy = 9, idr=1;
	double **mauk, *pt, s = 0.0, r = 0.0, p = 0.0, rez = 0.0, rap = 0.0;
	double *ao = new double[2 * N], *po = NULL, d = 0.0, hb = 5e-1, hb0 = hb, hf = 1e0;
	for (j = 0; j < ks; j++) d = d + hf;
	double de = hf, ta = 0.0, tb = 0.0, tc = 0.0, fa = 0.0, fb = 0.0, fc = 0.0;
	double ba = bnv, bb = bkv, bc = (ba + bb) / 2e0, xk = 0.0, dxk = 1e-3;
	double ca = cnv, cb = ckv, cc = (ca + cb) / 2e0, hc = 1e-1, hc0 = hc, tt = 0.0;
	double mk = 1e2, mkm = mk, c00, b0 = bc, c0 = cc, kdc = 5e0, kdb = 5e0, tminver, tmaxver;
	if (!ao) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < (2 * N); k++) ao[k] = 0.0;
	qobver = opredKTPTKTochSha(qobv, etev, etev[vtv], cemv);
	cout << "qo = " << qobver << "\n";
	if ((vtv == vtvn) && (sctxv <= 2) && (!temrasv)) {
		j = nxtv - nnxtv; temrasv = new double[j];
		if (!temrasv) { cout << snmv << endl; k = getchar(); exit(1); }
		ta = tholv[vtv] - tgorv[vtv]; ta = ta / y0ver; tb = tgorv[vtv];
		for (k = 0; k < j; k++) { temrasv[k] = tb + ta*xk; xk = xk + dxk; }
	}
	if ((vtv == vtvn) && (sctxv <= 2) && (!ktrv)) {
		ktrv = KoefPoglRosselNac(etev, idr, cemv, smgv, ssiv - pksvv, salv + pksvv, dkoscvt, dkoscvm, dkoscvl, dkospv, dkoalv, tkuscv, kuscv, dmkooscv, idr, 1);
		for (k = 0; k < cemv; k++) cout << "ktr = " << ktrv[k] << "\ttem = " << etev[k] << endl;
	} 
	l = vtv; tminver = tholv[l]; tmaxver = tgorv[l];
	ta = tminver; tb = tmaxver; tc = temrasv[sctxv];
	mauk = RaschRTAVer(ks, hver, 0.0, 0.0, 1, vtv, hvove, 0, tc, 0, 0, 1);
	k = 0; Raver = mauk[k]; k++; Taver = mauk[k]; k++; Aaver = mauk[k]; k++; Rver = mauk[k]; k++; Tver = mauk[k]; k++;
	Aver = mauk[k]; k++; Rtver = mauk[k]; k++; Ttver = mauk[k]; k++; Atver = mauk[k]; delete[]mauk;
	ca = cnv; cb = ckv; m = 0; q = 0;
	while ((m<mma) && (fabs(cb - ca)>1e-3) && (!q)) {
		cc = (ca + cb) / 2e0;
		dpctv=opredKTPTKTochSha(dpctvm, etev, tc, cemv);
		pt = FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, cc, vy, qobv, etev, cemv, 0, snmv, dpctv, vtv, ektpv, y0ver, ktrv, kttkv, ks, ktpvover, vtever, dmkvover, ecktpver, hver, hvove, Rver, Tver, Aver, Rtver, Ttver, Atver, txver, qxver, dtxver, sctxv, tocrasver, porver, 2 * N, stchsrver, Tpctv, b0, c0, mk, svfdv, 0, hkver);
		k = 3; fc = pt[k]; delete[]pt;
		if (fabs(fc) < 1e0) { c0 = cc; break; }
		pt = FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, ca, vy, qobv, etev, cemv, 0, snmv, dpctv, vtv, ektpv, y0ver, ktrv, kttkv, ks, ktpvover, vtever, dmkvover, ecktpver, hver, hvove, Rver, Tver, Aver, Rtver, Ttver, Atver, txver, qxver, dtxver, sctxv, tocrasver, porver, 2 * N, stchsrver, Tpctv, b0, c0, mk, svfdv, 0, hkver);
		k = 3; fa = pt[k]; delete[]pt;
		if (fabs(fa) < 1e0) { c0 = ca; break; }
		pt = FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, cb, vy, qobv, etev, cemv, 0, snmv, dpctv, vtv, ektpv, y0ver, ktrv, kttkv, ks, ktpvover, vtever, dmkvover, ecktpver, hver, hvove, Rver, Tver, Aver, Rtver, Ttver, Atver, txver, qxver, dtxver, sctxv, tocrasver, porver, 2 * N, stchsrver, Tpctv, b0, c0, mk, svfdv, 0, hkver);
		k = 3; fb = pt[k]; delete[]pt;
		if (fabs(fb) < 1e0) { c0 = cb; break; }
		if ((fa*fc) < 0.0) { cb = cc; q = 0; }
		if ((fb*fc) < 0.0) { ca = cc; q = 0; }
		if (((fb*fc)<0.0) && ((fa*fc)<0.0)) { cb = cc; q = 0; }
		if (((fb*fc)>0.0) && ((fa*fc)>0.0)) q = 1;
		m++;
	}
	c0 = cc; c00 = c0; cov = c00; cout << "c0 = " << c0 << "\ttb = " << tb << "\ttc = " << tc; 
	m = 0; ba = bnv; bb = bkv; l = vtv; ta = tholv[l]; tb = tgorv[l]; tc = temrasv[sctxv];
	while ((m<mma) && ((bb - ba)>1e-1))
	{
		bc = (ba + bb) / 2e0; b0 = bc; c0 = c00; cc = c0; dpctv=opredKTPTKTochSha(dpctvm, etev, tc, cemv);
		pt = FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, c0, vy, qobv, etev, cemv, 0, snmv, dpctv, vtv, ektpv, y0ver, ktrv, kttkv, ks, ktpvover, vtever, dmkvover, ecktpver, hver, hvove, Rver, Tver, Aver, Rtver, Ttver, Atver, txver, qxver, dtxver, sctxv, tocrasver, porver, 2 * N, stchsrver, Tpctv, b0, c0, mk, svfdv, 0, hkver);
		k = 9; s = pt[k]; k = 7; mk = pt[k]; k = 5; b0 = pt[k]; k = 4; tt = pt[k]; k = 10; r = pt[k]; delete[]pt;
		if ((mk<mkm) && (mk>0.0)) { mkm = mk; bov = bc; b0 = bc; lxver[sctxv] = r; ooxver[sctxv] = s; txver[sctxv] = tt; qxver[sctxv] = qobver; } //cout << "\tbc = " << bc << "\tmk = " << mk << "\tmkm = " << mkm << endl;
		m++; if (bc < 1e0) hb = hb0 / kdb; else hb = hb0; bb = bb - hb;
	} cout << "\tbop = " << bov << "\tmop = " << mkm << "\tcop = " << cov << endl;
	if (fabs(tnrv - tt)>1e-1)
		tnrv = tt; else {
			cnv = cov - hc0; if (cnv <= 0.0) cnv = 1e-2;
			ckv = cov + 2e0*hc0;
			bnv = bov - hb0; if (bnv <= 0.0) bnv = 1e-2;
			bkv = bov + 2e0*hb0;
			tnrv = tt;
		}
	delete[]ao; return rez;
}
double *opredKoefOtr(double *tem, double *ko, int kost, double ep, int vyve, int cemve, double *mkove, double *eteve, int vyma)
{
	int k = 0, j, q, f;
	q = 1; for (k = 0; k < cemve; k++) if (mkove[k] < ep) q = 0;
	if (!q) for (k = 0; k < cemve; k++) {
		if (vyve) mkove[k] = ReflSredVer(eteve[k]); }
	for (k = 0; k < kost; k++) ko[k] = opredKTPTKTochSha(mkove, eteve, tem[k], cemve);
	return ko;
}
double *zadrktVer(int zf, int kost, double d, int vyte, double htk, double hvo, int prod, int vy, double tc, int c, int u, double *rta)
{
	int j; double **mu = RaschRTAVer(kost, htk, 0.0, 0.0, 1, vyte, hvo, 0, tc, 1, 0, 0), *tt = mu[0], *te = new double[kost];
	if (!te) { cout << snmv << endl; j = getchar(); exit(1); }
	else for (j = 0; j < kost; j++) te[j] = tt[j];
	if (!u) {
		int k = 0, kst, q = 0, m = 4, b;
		double *prs = new double[m*kost*kost], *pr = new double[m], Er;
		if ((!prs) || (!pr)) { cout << snmv << endl; k = getchar(); exit(1); }
		for (k = 0; k < (m*kost*kost); k++) prs[k] = 0.0;
		q = 0; for (kst = 1; kst <= kost; kst++)
		{
			for (b = 0; b < m; b++) pr[b] = 0.0;
			for (k = 1; k <= kost; k++) {
				pr = izstNVer(k, kst, m, kost);
				for (b = 0; b < m; b++) prs[q + b] = pr[b]; q = q + m;
			}
		}
		if (prod == 1) Er = opredLuchSostVer(prs, te, q, zf, vyte);
		delete[]pr; delete[]te; delete[]tt; delete[]mu; return prs;
	}
	else { delete[]tt; delete[]mu; return te; }
}
double **izmRTAVer(double *tere, int kost, int izm, double *Ra, double *Ta, double *Aa, double *Rb, double *Tb, double *Ab, int v) //izm = 0 - нет изменений, izm - учитываются изменения //поиск изменения степени черноты или безразмерного коэффициента поглощения
{
	double **mu, ko, dkoscv; int k, rt = dmkooscv;
	for (k = 0; k < kost; k++) {
		ko = opredKTPTKTochSha(kuscv, tkuscv, tere[k], rt);
		dkoscv = opredKTPTKTochSha(dkoscvm, dkoscvt, tere[k], dkoscvl);
		if ((ko<0.0) || (ko>1e0) || (dkoscv<0.0) || (dkoscv>1e0) || (!izm)) { dkoscv = 1e0; ko = 1e0; }
		Aa[k] = Aa[k] * ko*dkoscv; Ta[k] = 1e0 - Aa[k] - Ra[k]; Tb[k] = Ta[k]; Ab[k] = Aa[k];
	}
	mu = chaRTAVer(kost, Ra, Ta, Aa, Rb, Tb, Ab, v);
	return mu;
}
double **chaRTAVer(int kost, double *Ra, double *Ta, double *Aa, double *Rb, double *Tb, double *Ab, int v)
{
	int f = 9, k = kost; double *tmp = new double[k], **mu = new double*[f]; 
	if ((!tmp) || (!mu)) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < kost; k++) { if (Ta[k] * Ra[k] >= 1e0) v = 3; if (v == 3) { Ab[k] = Aa[k]; Rb[k] = Ra[k]; Tb[k] = Ta[k]; } }
	if (v < 3) {
		for (k = 0; k < kost; k++) {
			Ab[k] = (1e0 - Ta[k] + Ra[k] * Ta[k] - Ra[k]); Ab[k] = Ab[k] / (1e0 - Ta[k] * Ra[k]);
			tmp[k] = pow((1e0 - Ra[k])*Ta[k], 2e0)*Ra[k] / (1e0 - pow(Ra[k] * Ta[k], 2e0)) + Ra[k]; Rb[k] = tmp[k];
			Tb[k] = pow(1e0 - Ra[k], 2e0)*Ta[k] / (1e0 - pow(Ra[k] * Ta[k], 2e0));
		}
	}
	delete[]tmp; k = 0; mu[k] = Ra; k++; mu[k] = Ta; k++; mu[k] = Aa; k++; mu[k] = Rb; k++; mu[k] = Tb; k++; mu[k] = Ab; k++;
	return mu;
}
double *kopoVer(double ta, double tb, int vyb, double tl, int kost, int vyte, double htk, double hvo, int w, double *al)
{
	int k, q = 100, j, r, f;
	double t0 = 2e1, te = tev0 + t0, dte = 1e0, *p=NULL, e = 1e-3, t, a1, a2, t1, t2;
	char *s = new char[q]; KoePog *kp = new KoePog, *ne=NULL, *roo=NULL, *pre=NULL;
	if ((!kp) || (!s)) { cout << snmv << endl; k = getchar(); exit(1); } for (j = 0; j < q; j++) s[j] = '\0';
	ifstream fin; fin.open(sfatv); if (!fin.is_open()) { cout << sfnov << endl; k = getchar(); exit(1); }
	roo = kp; k = 0; while (!fin.eof()) {
		fin.getline(s, q, '\n'); ne = new KoePog; if (!ne) { cout << snmv << endl; j = getchar(); exit(1); }
		kp->alp = atof(s)*dkoalv*dkospv; kp->tem = te; kp->nex = ne; pre = kp; kp = ne; k++; te = te + dte; kp->nex = NULL;
	}
	delete[]ne; kp = NULL; pre->nex = kp; fin.close(); delete[]s; r = k;
	if ((!vyb) || (vyb == 2)) {
		f = 2; double *xi = new double[kost], *teks = new double[kost], *koe = new double[f], knat;
		if ((!teks) || (!koe) || (!xi)) { cout << snmv << endl; k = getchar(); exit(1); }
		teks = oprRasTemNachVer(cemv, dmkov, teks, koe, kost, xi, htk, hvo, hkver, ta, tb, w);
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
		j = 1; p = new double[j]; if (p) p[0] = t; else { cout << snmv << endl; j = getchar(); exit(1); }
	}
	kp = roo; while (kp) { ne = kp->nex; delete kp; kp = ne; } //удаление списка
	return p;
}
double *izstNVer(int izst, int kst, int l, int ocs)
{
	int k; double *o;
	if (abs(izst - kst) <= N) o = podschchieleSha(izst, kst, ocs, Rver, Tver);
	else {
		o = new double[l]; if (!o) { cout << "No memory!"; k = getchar(); exit(1); }
		for (k = 0; k < l; k++) o[k] = 0.0;
	}
	return o;
}
double BolTochRasAlpha(int vyte, int kost, double hvp, double htk, double tc, double *ktpvo, double *vte, double *ete, int ce, int dmk, char *snm, double *Tra, double *Ref, int nao, double *ktptk, double *Ab, int vyv, double *Refa, double *Traa, double *Aba, double *Reft, double *Trat, double *Abt)
{
	int j, l = nao, k; double *prot = NULL, hlr0 = 1e0, hrl0 = 0.0, alks, ra = fabs(hlr0 - hrl0), tol = 1e0;
	double *ao = new double[l], *Ts = NULL, *Tss = NULL, *Tsr = NULL, *Tna = NULL, **mu, *po, *sislr = NULL, *sisrl = NULL;
	double *reiz, *hrl1 = NULL, *hlr1 = NULL, d = 0.0, slr = 0.0, srl = 0.0, hf = 1e0; for (j = 0; j < kost; j++) d = d + hf;
	if (!vyv) {
		prot = zadrktVer(1, kost, d, 0, htk, hvp, 0, 0, tc, 1, 0, prot); //vyv - выбоор вещества: 0 - вермикулит, 1 - шамот
		mu = RaschRTAVer(kost, htk, 0.0, 0.0, 1, 0, hvp, 0, tc, 0, 0, 1); }
	else { cout << "Net takogo veschestva!" << endl; k = getchar(); exit(1); }
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
double RasFracXeffVer60(int v)
{
	int l = 6, k; double *mkbr = new double[l], *mv = new double[l], *tol = new double[l];
	k = 0;	mkbr[k] = 250.239; k++; mkbr[k] = 249.740; k++; mkbr[k] = 249.223; k++;
	mkbr[k] = 250.395; k++; mkbr[k] = 250.336; k++; mkbr[k] = 249.55;
	k = 0;	mv[k] = 0.283; k++; mv[k] = 0.464; k++; mv[k] = 0.812; k++;
	mv[k] = 0.22; k++; mv[k] = 0.547; k++; mv[k] = 0.777;
	k = 0;	tol[k] = 0.73; k++; tol[k] = 0.72; k++; tol[k] = 0.72; k++;
	tol[k] = 0.72; k++; tol[k] = 0.71; k++; tol[k] = 0.7;
	double rokbr = 2.75, rov = 0.49, *xv = EffectTols(mkbr, mv, tol, rov, rokbr, l), xsr = 0.0;
	if (!v) xsr = (xv[0] + xv[3]) / 2e0;
	if (v == 1) xsr = (xv[1] + xv[4]) / 2e0;
	if (v == 2) xsr = (xv[2] + xv[5]) / 2e0;
	delete[]mkbr; delete[]mv; delete[]tol; delete[]xv; return xsr;
}
double RasFracXeffVer60_100(int v)
{
	int l = 9, k; double *mkbr = new double[l], *mv = new double[l], *tol = new double[l], xsr;
	k = 0;	mkbr[k] = 250.0;	k++; mkbr[k] = 249.629; k++; mkbr[k] = 249.294; k++; mkbr[k] = 249.706; k++;
	mkbr[k] = 249.51; k++; mkbr[k] = 249.307; k++; mkbr[k] = 250.328; k++; mkbr[k] = 249.604; k++; mkbr[k] = 249.206;
	k = 0;	mv[k] = 0.255;	k++; mv[k] = 0.539;	  k++; mv[k] = 0.809;		k++; mv[k] = 0.295;	  k++;
	mv[k] = 0.517;	k++; mv[k] = 0.756;	  k++; mv[k] = 0.36;		k++; mv[k] = 0.534;	  k++; mv[k] = 0.843;
	k = 0;	tol[k] = 0.72;	k++; tol[k] = 0.71;	  k++; tol[k] = 0.7;		k++; tol[k] = 0.7;	  k++;
	tol[k] = 0.73;	k++; tol[k] = 0.72;	  k++; tol[k] = 0.74;	    k++; tol[k] = 0.7;	  k++; tol[k] = 0.76;
	double rokbr = 2.75, rov = 0.52, *xv = EffectTols(mkbr, mv, tol, rov, rokbr, l);
	if (!v) xsr = (xv[0] + xv[3] + xv[6]) / 3e0;
	if (v == 1) xsr = (xv[1] + xv[4] + xv[7]) / 3e0;
	if (v == 2) xsr = (xv[2] + xv[5] + xv[8]) / 3e0;
	delete[]mkbr; delete[]mv; delete[]tol; delete[]xv; return xsr;
}
double RasFracXeffVer100_150(int v)
{
	int l = 9, k; double *mkbr = new double[l], *mv = new double[l], *tol = new double[l], xsr;
	k = 0;	mkbr[k] = 249.913; k++; mkbr[k] = 249.607; k++; mkbr[k] = 249.218; k++; mkbr[k] = 249.929; k++;
	mkbr[k] = 249.695; k++; mkbr[k] = 249.306; k++; mkbr[k] = 250.405; k++; mkbr[k] = 249.625; k++; mkbr[k] = 249.348;
	k = 0;	mv[k] = 0.315;	 k++; mv[k] = 0.473;     k++; mv[k] = 0.709;     k++; mv[k] = 0.293;     k++;
	mv[k] = 0.528;     k++; mv[k] = 0.83;      k++; mv[k] = 0.27;      k++; mv[k] = 0.493;     k++; mv[k] = 0.764;
	k = 0;	tol[k] = 0.74;     k++; tol[k] = 0.74;     k++; tol[k] = 0.72;     k++; tol[k] = 0.72;     k++;
	tol[k] = 0.71;     k++; tol[k] = 0.7;      k++; tol[k] = 0.78;     k++; tol[k] = 0.73;     k++; tol[k] = 0.76;
	double rokbr = 2.75, rov = 0.53, *xv = EffectTols(mkbr, mv, tol, rov, rokbr, l);
	if (!v) xsr = (xv[0] + xv[3] + xv[6]) / 3e0;
	if (v == 1) xsr = (xv[1] + xv[4] + xv[7]) / 3e0;
	if (v == 2) xsr = (xv[2] + xv[5] + xv[8]) / 3e0;
	delete[]mkbr; delete[]mv; delete[]tol; delete[]xv; return xsr;
}
double RasFracXeffVer150_200(int v)
{
	int l = 6, k; double *mkbr = new double[l], *mv = new double[l], *tol = new double[l], xsr;
	k = 0;	mkbr[k] = 250.882; k++; mkbr[k] = 249.590; k++; mkbr[k] = 249.213; k++;
	mkbr[k] = 250.299; k++; mkbr[k] = 249.441; k++; mkbr[k] = 249.365;
	k = 0;	mv[k] = 0.320;	 k++; mv[k] = 0.533;	   k++; mv[k] = 0.849;	 k++;
	mv[k] = 0.223;	 k++; mv[k] = 0.502;	   k++; mv[k] = 0.797;
	k = 0;	tol[k] = 0.76;	 k++; tol[k] = 0.72;	   k++; tol[k] = 0.69;	 k++;
	tol[k] = 0.73;	 k++; tol[k] = 0.73;	   k++; tol[k] = 0.73;
	double rokbr = 2.75, rov = 0.56, *xv = EffectTols(mkbr, mv, tol, rov, rokbr, l);
	if (!v) xsr = (xv[0] + xv[3]) / 2e0;
	if (v == 1) xsr = (xv[1] + xv[4]) / 2e0;
	if (v == 2) xsr = (xv[2] + xv[5]) / 2e0;
	delete[]mkbr; delete[]mv; delete[]tol; delete[]xv; return xsr;
}
double *EffectTols(double *mkbr, double *mv, double *tol, double rov, double rokbr, int n)
{
	int k; double *vkbr = new double[n], *vv = new double[n], *xvo = new double[n], t=1e3;
	for (k = 0; k < n; k++) {
		vkbr[k] = mkbr[k] / (t*rokbr);
		vv[k] = mv[k] / (t*rov);
		xvo[k] = (vv[k] / (vv[k] + vkbr[k]))*tol[k] * t;
	}
	delete[]vkbr; delete[]vv; return xvo;
}
double RaschAlphaTvKarVer()
{
	int l = 3, k, fr = 0, nfr = 4, d = 10, q, cfk = l*nfr, cei;
	long j, jk = 100 * 100 * 100, pj, h;
	double x = 0.0, y = 0.0, p = 0.0, xt = 0.0, yp = 0.0, hf = 1e0, tc = 22.0 + tev0, altc;
	double *srra = rasPorpoRazVer(porver, vyfv, 3, vysv, isrp, vpkf), *porv = new double[nfr];
	double *pn = new double[d], *alx = new double[d], *rcf = new double[cfk];
	double *rapo = rasPorpoRazVer(porver, vyfv, 1, vysv, isrp, vpkf), **unau;
	double srp = rapo[0], ce, cet, *pp = NULL, lamtem; delete[]rapo;
	double *legr = rasPorpoRazVer(porver, vyfv, 4, vysv, isrp, vpkf), dvpn = 1448.0*(1e-6) / tc / 4e0; //меньше 1 %
	rapo = rasPorpoRazVer(porver, vyfv, 2, vysv, isrp, vpkf); ce = rapo[0]; cet = ce; delete[]rapo;
	j = 0; while (cet > 0.0) { cet = cet - hf; j++; } cei = j;
	rapo = rasPorpoRazVer(porver, vyfv, 0, vysv, isrp, vpkf); //cout << "cem_srp = " << cei << "\tsrp = " << srp << "\t"; cout << endl; //for (j=0; j<cei; j++) if (j<10) cout << "j = " << j << "\trpr = " << rapo[j] << "\tlegr = " << legr[j] << "\t"; cout << endl;
	if (porv) {
		x = 0.0; yp = 0.0; xt = 0.0;
		for (j = 0; j < cei; j++) {
			if (yp <= 3e1) x = x + srra[j] * rapo[j];
			if (yp <= 8e1) xt = xt + srra[j] * rapo[j]; yp = yp + hf;
		}
		porv[0] = x*porver / srp; porv[1] = xt*porver / srp; porv[2] = porver; porv[3] = porver;
	}
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
	delete[]rapo; delete[]srra; delete[]rcf; delete[]pn; delete[]alx; delete[]alsf; delete[]legr; delete[]porv;
	x = x / altc; printf("sr_okp = %0.10lf\n", x); //ослабление КП за счет пористой структуры вермикулита
	return x;
}
double **RaschRTAVer(int kost, double htm, double kta, double ktb, int izm, int vyte, double hvo, int v, double ti, int c, int w, int u)
{
	int k = kost, j, r; double *te = NULL, e = tocrasver, **mauk=NULL;
	if (!v) {
		alphaVer = kopoVer(kta, ktb, 2, 0.0, kost, vyte, htm, hvo, w, alphaVer);
		te = kopoVer(kta, ktb, 0, 0.0, kost, vyte, htm, hvo, w, alphaVer);
		if (!c) {
			k = kost;
			if (Raver) delete[]Raver; Raver = new double[k];
			if (Rver) delete[]Rver; Rver = new double[k];
			if (Rtver) delete[]Rtver; Rtver = new double[k];
			k = kost; Raver = opredKoefOtr(te, Raver, k, tocrasver, 1, cemv, mkov, etev, 0);
			for (j = 0; j < k; j++) { Rtver[j] = Raver[j]; Rver[j] = Raver[j]; }
			if (Taver) delete[]Taver; Taver = new double[k];
			if (Ttver) delete[]Ttver; Ttver = new double[k];
			if (Tver) delete[]Tver; Tver = new double[k];
			if (Aaver) delete[]Aaver; Aaver = new double[k];
			if (Atver) delete[]Atver; Atver = new double[k];
			if (Aver) delete[]Aver; Aver = new double[k];
			k = kost; for (j = 0; j < k; j++) {
				Taver[j] = -alphaVer[j] * htm; Taver[j] = exp(Taver[j]); Aaver[j] = 1e0 - Taver[j] - Raver[j];
				if (Aaver[j] < 0.0) { Aaver[j] = alphaVer[j] * htm; Taver[j] = 1e0 - Raver[j] - Aaver[j]; }
				Atver[j] = Aaver[j]; Ttver[j] = Taver[j]; Aver[j] = Aaver[j]; Tver[j] = Taver[j];
			}
			k = kost; mauk = izmRTAVer(te, k, 1, Rver, Tver, Aver, Raver, Taver, Aaver, 0);
			r = 6; mauk[r] = Rtver; r++; mauk[r] = Ttver; r++; mauk[r] = Atver;
		}
		else if (c == 1) { r = 1; mauk = new double*[r]; 
		if (!mauk) { cout << snmv << endl; j = getchar(); exit(1); } mauk[0] = te; }
	}
	else if (v == 1) {
		te = kopoVer(0.0, 0.0, 1, ti, kost, c, htm, hvo, w, alphaVer); r = 1;
		mauk = new double*[r]; if (!mauk) { cout << snmv << endl; j = getchar(); exit(1); } mauk[0] = te;
	}
	if (!u) mauk[0] = te; return mauk;
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
double *oprEffDoliTepPerenVer(double ko, double d, double por)
{
	int k, j, f, kost = ks, ksu = 2 * ks; double hvozd = ko, htch, hf = 1e0, e = tocrasver;
	ecktpver = new double[kost]; if (!ecktpver) { cout << snmv << endl; k = getchar(); exit(1); }
	double *tepo = new double[ksu], *dol = new double[cemv], r, p, sa, sb, sc, t, fa, fb, fc;
	if ((!tepo) || (!dol)) { cout << snmv << endl; k = getchar(); exit(1); }
	for (j = 0; j < ksu; j++) tepo[j] = 0.0;
	hvozd = ko; htch = hvozd*(hf - por) / por;
	k = 1; for (j = 0; j < cemv; j++) if (qobv[j] < e) k = 0;
	if (!k) opredtemphcVer(ektpv, etev, tgorv, tholv, qobv, cemv, dmkov, y0ver); //cout << "hvoz = " << hvozd << "\tht = " << htch << endl;
	for (j = 0; j < cemv; j++) { //определяем ЭКТП многослойной стенки
		tepo = opredTempStenShaFragm(tepo, 2 * kost, ktpvover, vtever, etev, kttkv, cemv, dmkvover, htch, hvozd, qobv[j], etev[j], -hf); //в середине слоя
		r = 0.0; for (k = 0; k < ksu; k++) for (f = 0; f<ksu; f++) { p = tepo[k] - tepo[f]; if (p>r) r = p; }
		p = hvozd*(d - hf) + htch*d; r = r / p; t = qobv[j] / r; ecktpver[j] = t;
	}
	delete[]tepo;
	f = 1000; for (j = 0; j<cemv; j++)
	{
		sa = 0.0; sb = 1e0; k = 0;
		do {
			sc = (sa + sb) / 2e0;
			fa = kttkv[j] * sa + ecktpver[j] * (hf - sa) - ektpv[j]; //эффективные КТП многослойной стенки и  перемычки должны сравняться
			fb = kttkv[j] * sb + ecktpver[j] * (hf - sb) - ektpv[j]; //чтобы найти относительные доли площадей сечения переноса общей ПТП
			fc = kttkv[j] * sc + ecktpver[j] * (hf - sc) - ektpv[j];
			if ((fc*fb>0.0) && (fa*fc<0.0)) sb = sc; if ((fc*fa>0.0) && (fb*fc<0.0)) sa = sc;
			r = fabs(sa - sb); k++;
		} while ((r>e) && (k<f));
		dol[j] = sc;
	}
	return dol;
}
double *KorrZnachVozdProsVer(double hps, double ksf, double por, int vy)
{
	int j = 0, k = 1000, q=0; double pa = 1e-2, pb = 1e0, *po, pc, ra = fabs(pa - pb), e = tocrasver, hf = 1e0;
	double fa, fb, fc, ta, tb, tc, tca, tcb, tcc, ka = hps*pa / por, kb = pb*hps / por, kc, *pkc=new double[1]; //cout << "hps = " << hps << "\tksuf = " << ksf << "\tpor = " << por << endl;
	dpctvm=new double[cemv];
	for (q=0; q<cemv; q++) { j=0; pa = 1e-2; pb = 1e0; ra = fabs(pa - pb); ka = hps*pa / por; kb = pb*hps / por;
	while ((ra>e) && (j < k)) { //подтягиваем пористость к значению, которое задали изначально, во время подстройки ЭКТП
		pc = (pa + pb) / 2e0;
		kc = hps*pc / por;
		po = oprEffDoliTepPerenVer(kc, ksf, pc); tc = po[q]; delete[]po; //при 373 К
		tcc = kc*(hf - pc) / pc;
		fc = (hf - tc)*kc / (kc + tcc) - por;
		ka = hps*pa / por;
		po = oprEffDoliTepPerenVer(ka, ksf, pa); ta = po[q]; delete[]po; //определяем долю площади сечения перемычки
		tca = ka*(hf - pa) / pa;
		fa = (hf - ta)*ka / (ka + tca) - por;
		kb = hps*pb / por;
		po = oprEffDoliTepPerenVer(kb, ksf, pb); tb = po[q]; delete[]po; //через перемычку тепло распространяется чистой теплопроводностью
		tcb = kb*(hf - pb) / pb;
		fb = (hf - tb)*kb / (kb + tcb) - por;
		if ((fc*fb > 0.0) && (fa*fc<0.0)) pb = pc; if ((fc*fa>0.0) && (fb*fc < 0.0)) pa = pc;
		j++; ra = fabs(pa - pb);
	}
	dpctvm[q] = tc; cout << "Dol Plo CTP = " << tc << endl; }
	for (j=0; j<cemv; j++) cout << "EC KTP TK ( " << j << " ) = " << ecktpver[j] << endl;
	if (!vy) { pkc[0]=kc; return pkc; } //скорректированное значение толщины воздушной прослойки (размер поры), когда ввели перемычку
	else if (vy == 1) return dpctvm;
} //доля площади, через которую происходит перенос тепла чистой теплопроводностью
double *oprRasTemNachVer(int ce, int cee, double *teks, double *koe, int kost, double *xi, double htch, double hvozd, double hkover, double a, double b, int w)
{
	double e = tocrasver, hkx = 0.0, ht = 0.0; int k, j;
	if (!w) {
		k = 1; for (j = 0; j < ce; j++) if ((tgorv[j] < e) || (tholv[j] < e) || (qobv[j] < e)) { k = 0; break; }
		if (!k) opredtemphcVer(ektpv, etev, tgorv, tholv, qobv, ce, cee, y0ver);
		hkx = hkover; koe[0] = (tholv[vtv] - tgorv[vtv]) / y0ver; koe[1] = tgorv[vtv]; tmav = tgorv[vtv]; tmiv = tholv[vtv];
	}
	else if (w == 1) { hkx = y0ver / 2e0; koe[0] = a; koe[1] = b; tmav = b; tmiv = b - a*y0ver; }
	xi[0] = hkx + htch / 2e0; ht = hvozd + htch;
	for (k = 1; k < kost; k++) xi[k] = xi[k - 1] + ht; //массив середин каждой из стенок по толщине
	for (k = 0; k<kost; k++) teks[k] = koe[0] * xi[k] + koe[1]; //for (k=0; k<kost; k++) cout << "k = " << k << "\txk = " << xi[k] << "\ttex = " << teks[k] << endl; //линеаризация поля температур
	qobver = opredKTPTKTochSha(qobv, etev, (teks[0] + teks[kost - 1]) / 2e0, ce); return teks;
}
double KorrZnachVozdProsVermik(double hps, double ksf, double por)
{
	int j = 0, k = 1000, q=0; double pa = 1e-3, pb = 1e0, pc, ra = fabs(pa - pb); dpctv=opredKTPTKTochSha(dpctvm, etev, etev[vtv], cemv);
	double fa, fb, fc, tca, tcb, tcc, ka = hps*pa / por, kb = pb*hps / por, kc, e = tocrasver; //cout << "hps = " << hps << "\tksuf = " << ksf << "\tpor = " << por << endl;
	while ((ra>e) && (j < k)) { //подтягиваем пористость к значению, которое задали изначально, во время подстройки ЭКТП
		pc = (pa + pb) / 2e0;
		kc = hps*pc / por;
		tcc = kc*(1e0 - pc) / pc;
		fc = (1e0 - dpctv)*kc / (kc + tcc) - por;
		ka = hps*pa / por;
		tca = ka*(1e0 - pa) / pa;
		fa = (1e0 - dpctv)*ka / (ka + tca) - por;
		kb = hps*pb / por;
		tcb = kb*(1e0 - pb) / pb;
		fb = (1e0 - dpctv)*kb / (kb + tcb) - por;
		if ((fc*fb > 0.0) && (fa*fc<0.0)) pb = pc;
		if ((fc*fa>0.0) && (fb*fc < 0.0)) pa = pc;
		j++; ra = fabs(pa - pb);
	} //скорректированное значение толщины воздушной прослойки (размер поры), когда ввели перемычку
	return kc;
}
//---------
double *PoiskZavVelTemVer(int v, double *poin)
{
	int k, n = qtn; double *Te = new double[n], delT = 1e2, Tnacha = tnav, y0 = y0ver, cf = 0.0, hf = 1e0, nf = 0.0;
	if (!Te) { cout << snmv << endl; k = getchar(); exit(1); }
	Te[0] = Tnacha; for (k = 1; k < n; k++) Te[k] = Te[k - 1] + delT; //double *temvh=new double[n], *temvc=new double[n]; if ((!temvh) || (!temvc)) { cout << snmv << endl; k=getchar(); exit(1); } //for (k=0; k<n; k++) { temvc[k]=0.0; temvh[k]=0.0; } //double *ktpv=new double[n]; if (!ktpv) { cout << snmv << endl; k=getchar(); exit(1); } for (k=0; k<n; k++) { ktpv[k]=0.0; }
	double *temvs = new double[n], *tepv = new double[n]; if ((!tepv) || (!temvs)) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < n; k++) { temvs[k] = 0.0; tepv[k] = 0.0; }
	double *temvht, *temvct, *tepvt; //double *ktpvt;
	int vyfrve = 0, vysove = 0, vymivmf = 1, vyukve = 1, c = 0, nnvyfv = 2, nnvysv = 3, nnvmivmf = 2, nnvyuv = 2, kvf, jvsv, qvmi, qvuv;
	double *nvyfv = new double[nnvyfv], *nvysv = new double[nnvysv], *nvmivmf = new double[nnvmivmf], *nvyuv = new double[nnvyuv];
	if ((!nvyfv) || (!nvysv) || (!nvmivmf) || (!nvyuv)) { cout << snmv << endl; k = getchar(); exit(1); }
	nvyfv[0] = 0; nvyfv[1] = 1; nvysv[0] = 0; nvysv[1] = 1; nvysv[2] = 2; nvmivmf[0] = 1; nvmivmf[1] = 2; nvyuv[0] = 1; nvyuv[1] = 2;
	for (kvf = 0; kvf < nnvyfv; kvf++) {
		vyfrve = nvyfv[kvf];
		for (jvsv = 0; jvsv < nnvysv; jvsv++) {
			vysove = nvysv[jvsv];
			if (!vyfrve) {
				for (qvmi = 0; qvmi < nnvmivmf; qvmi++) {
					vymivmf = nvmivmf[qvmi];
					temvct = napMasEKTPVer(vyfrve, vysove, vymivmf, Te, 0, y0, vyukve, n, 1);
					temvht = napMasEKTPVer(vyfrve, vysove, vymivmf, Te, 1, y0, vyukve, n, 1); //ktpvt=napMasEKTPVer(vyfrve,vysove,vymivmf,Te,3,y0,vyukve,n,1);
					tepvt = napMasEKTPVer(vyfrve, vysove, vymivmf, Te, 2, y0, vyukve, n, 1);
					for (k = 0; k < n; k++) { cf = (temvct[k] + temvht[k]) / 2e0; temvs[k] = temvs[k] + cf; cf = tepvt[k]; tepv[k] = tepv[k] + cf; } //for (k=0; k<n; k++) { ktpv[k]=ktpv[k]+ktpvt[k]; temvc[k]=temvc[k]+temvct[k]; temvh[k]=temvht[k]+temvh[k];}
					delete[]temvht; delete[]temvct; delete[]tepvt; //delete []ktpvt; 
					c++;
				}
			}
			else if (vyfrve == 1) {
				for (qvuv = 0; qvuv < nnvyuv; qvuv++) {
					vyukve = nvyuv[qvuv];
					temvct = napMasEKTPVer(vyfrve, vysove, vymivmf, Te, 0, y0, vyukve, n, 1);
					temvht = napMasEKTPVer(vyfrve, vysove, vymivmf, Te, 1, y0, vyukve, n, 1); //ktpvt=napMasEKTPVer(vyfrve,vysove,vymivmf,Te,3,y0,vyukve,n,1); 
					tepvt = napMasEKTPVer(vyfrve, vysove, vymivmf, Te, 2, y0, vyukve, n, 1);
					for (k = 0; k < n; k++) { cf = (temvct[k] + temvht[k]) / 2e0; temvs[k] = temvs[k] + cf; cf = tepvt[k]; tepv[k] = tepv[k] + cf; } //for (k=0; k<n; k++) { ktpv[k]=ktpv[k]+ktpvt[k]; temvc[k]=temvc[k]+temvct[k]; temvh[k]=temvht[k]+temvh[k]; }
					delete[]temvht; delete[]temvct; delete[]tepvt; //delete []ktpvt; 
					c++;
				}
			}
		}
	}
	delete[]nvyfv; delete[]nvysv; delete[]nvmivmf; delete[]nvyuv;
	cf = 0.0; for (k = 0; k < c; k++) cf = cf + hf; //if (cf>0.0) { for (k=0; k<c; k++) { temvc[k]=temvc[k]/cf; temvh[k]=temvh[k]/cf; temvs[k]=(temvc[k]+temvh[k])/2e0; } } //for (k=0; k<c; k++) ktpv[k]=ktpv[k]/cf; //for (k=0; k<n; k++) { cout << "T_hol ( " << k << " ) = " << temvc[k] << endl; } for (k=0; k<n; k++) { cout << "T_gor ( " << k << " ) = " << temvh[k] << endl; } for (k=0; k<n; k++) { cout << "T_sr ( " << k << " ) = " << temvs[k] << endl; } for (k=0; k<n; k++) { cout << "KTP ( " << k << " ) = " << ktpv[k] << endl; } for (k=0; k<n; k++) { cout << "PTP ( " << k << " ) = " << tepv[k] << endl; } 
	if (cf>0.0) { for (k = 0; k < n; k++) { tepv[k] = tepv[k] / cf; temvs[k] = temvs[k] / cf; } } delete[]Te; //delete []temvc; delete []temvh; 
	cf = 0.0; for (k = 0; k < n; k++) cf = cf + hf;
	if (!v) { for (k = 0; k < n; k++) poin[k] = tepv[k]; }
	else if (v == 1) { for (k = 0; k < n; k++) poin[k] = temvs[k]; }
	else if (v == 2) { poin[0] = cf; }
	delete[]tepv; delete[]temvs; return poin;
}
void opredtemphcVer(double *efktpv, double *eftev, double *tgv, double *thv, double *qon, int dlma, int n, double h) //n=3 - длина массива коэффициентов приближающего многочлена, tgv - температура горячей стенки, thv - температура холодной стенки, dlma - длина массива ЭКТП
{
	int k = 0, j; double *koeq = new double[n], *kho = new double[n], *kgo = new double[n];
	double *tsv = new double[n], t = tev0, ts, g, p, hf = 1e0, r, s;
	if ((!koeq) || (!kho) || (!kgo) || (!tsv)) { cout << snmv << endl; j = getchar(); exit(1); }
	for (k = 0; k < n; k++) tsv[k] = (tgv[k] + thv[k]) / 2e0;
	koeq = koefPribSha(qon, tsv, n, koeq);
	kgo = koefPribSha(tgv, tsv, n, kgo);
	kho = koefPribSha(thv, tsv, n, kho);
	for (k = 0; k < dlma; k++) {
		ts = eftev[k];
		g = 0.0; p = 0.0; for (j = 0; j < n; j++) { g = g + pow(ts, p)*koeq[j]; p = p + hf; } if (g<0.0) g = 0.0; if (qobv) qobv[k] = g;
		p = efktpv[k]; if (fabs(p)>0.0) { g = fabs(g*h / p / 2e0); tgorv[k] = eftev[k] + g; tholv[k] = eftev[k] - g; }
		else {
			r = 0.0; s = 0.0; for (j = 0; j < n; j++) { r = r + pow(ts, s)*kgo[j]; s = s + hf; } 
			if (r < 0.0) r = 0.0; if (tgorv) tgorv[k] = r;
			r = 0.0; s = 0.0; for (j = 0; j < n; j++) { r = r + pow(ts, s)*kho[j]; s = s + hf; } 
			if (r < 0.0) r = 0.0; if (tholv) tholv[k] = r;
		}
	}
	delete[]tsv; delete[]koeq; delete[]kho; delete[]kgo;
}
double *opredTempHolGorVer(double *ektp, double *ete, int n, int l, double h0, int v) //моделирование процесса теплообмена в вермикулитовом образце
{
	int k, j, w = 0; double del = 1e0, g, p, *qon, *tena, nit = 1e6, hf = 1e0, *po = new double[1], nf;
	double ep = 1e-2, d = 1e-3, Thna = tnoscv, Thnac = 0.0, *koeq = new double[l], ts, laef, dt = 1e0;
	if ((!koeq) || (!po)) { cout << snmv << endl; k = getchar(); exit(1); }
	po = PoiskZavVelTemVer(2, po); nf = po[0]; w = 0; while (nf>0.0) { nf = nf - hf; w++; }
	qon = new double[w]; tena = new double[w];
	if ((!tena) || (!qon)) { cout << snmv << endl; k = getchar(); exit(1); }
	qon = PoiskZavVelTemVer(0, qon); tena = PoiskZavVelTemVer(1, tena); //for (k = 0; k < w; k++) cout << "q = " << qon[k] << "\tt = " << tena[k] << endl;
	for (k = 0; k < l; k++) koeq[k] = 0.0;
	koeq = koefPribSha(qon, tena, w, koeq);
	double *tgc = new double[n], *thc = new double[n], *qobc = new double[n], *tsc = new  double[n];
	if ((!tgc) || (!thc) || (!qobc) || (!tsc)) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < n; k++) { tgc[k] = 0.0; thc[k] = 0.0; qobc[k] = 0.0; tsc[k] = 0.0; }
	if (n <= cemv) {
		for (k = 0; k < n; k++) {
			ts = ete[k]; g = 0.0; p = 0.0;
			for (j = 0; j<l; j++) { g = g + pow(ts, p)*koeq[j]; p = p + hf; } qobv[k] = g;
		} //for (k=0; k<n; k++) cout << "qob = " << qobv[k] << "\tektp = " << ektp[k] << "\tte = " << ete[k] << endl; 
		delete[]koeq; delete[]qon; delete[]tena; delete[]po;
		g = 0.0; for (k = 0; k<n; k++) {
			if (qobv[k]>0.0) {
				laef = ektp[k]; p = 0.0; Thnac = Thna + g*dt; del = 1e0;
				while ((del>ep) && (p < nit)) {
					tholv[k] = Thnac + p*d; //Tg - массив температур горячих стенок, Th - массив температур холодных стенок
					tgorv[k] = tholv[k] + qobv[k] * h0 / laef;
					del = (2e0*ete[k] - (tholv[k] + tgorv[k]));
					p = p + hf;
				} g = g + hf;
			}
			else { tholv[k] = 0.0; tgorv[k] = 0.0; }
		}
		for (k = 0; k < n; k++) tsredver[k] = (tgorv[k] + tholv[k]) / 2e0;
		for (k = 0; k < n; k++) { tgc[k] = tgorv[k]; thc[k] = tholv[k]; qobc[k] = qobv[k]; tsc[k] = tsredver[k]; }
	} //for (k=0; k<n; k++) cout << "t_h = " << thc[k] << "\tt_g = " << tgc[k] << "\tq_o = " << qobc[k] << "\tts = " << tsc[k] << "\tktp = " << ektp[k] << endl;
	if (!v) { delete[]thc; delete[]qobc; delete[]tsc; return tgc; }
	else if (v == 1) { delete[]tgc; delete[]qobc; delete[]tsc; return thc; }
	else if (v == 2) { delete[]thc; delete[]tgc; delete[]tsc; return qobc; }
	else if (v == 3) { delete[]thc; delete[]tgc; delete[]qobc; return tsc; }
}
void napMasEKTPVerNac(double wmg, double wsi, double wal)
{
	kektpv = napMasEKTPVer(vyfv, vysv, vmivmf, etev, 3, y0ver, vyuv, cemv, 0);
	int k, nk, j; double hf = 1e0;
	double *tholvs = napMasEKTPVer(vyfv, vysv, vmivmf, etev, 0, y0ver, vyuv, cemv, 1);
	double *tgorvs = napMasEKTPVer(vyfv, vysv, vmivmf, etev, 1, y0ver, vyuv, cemv, 1);
	double *qobvs = napMasEKTPVer(vyfv, vysv, vmivmf, etev, 2, y0ver, vyuv, cemv, 1);
	double *ektpvs = napMasEKTPVer(vyfv, vysv, vmivmf, etev, 3, y0ver, vyuv, cemv, 1);
	double *tsredvers = napMasEKTPVer(vyfv, vysv, vmivmf, etev, 4, y0ver, vyuv, cemv, 1);
	delete[]tholv; delete[]tgorv; delete[]qobv; delete[]tsredver; delete[]ektpv;
	tholv = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, 0, cemv);
	tgorv = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, 1, cemv);
	qobv = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, 2, cemv);
	tsredver = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, 3, cemv);
	ektpv = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, 4, cemv);
	double nf = 0.0, *p = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, 5, cemv), tn, g;
	delete[]tholvs; delete[]tgorvs; delete[]qobvs; delete[]ektpvs; delete[]tsredvers;
	for (k = 0; k<cemv; k++) tsredver[k] = (tholv[k] + tgorv[k]) / 2e0;
	nf = p[0]; delete[]p; k = 0; while (nf>0.0) { nf = nf - hf; k++; } nk = k; cemv = nk; tn = tsredver[0]; tnav = tn; //for (k=0; k<cemv; k++) cout << "t_hol = " << tholv[k] << "\tt_gor = " << tgorv[k] << "\tq_ob = " << qobv[k] << "\tts = " << tsredver[k] << "\tktp = " << ektpv[k] << endl;
	if (etev) delete[]etev; //cout << "cemv = " << cemv << endl;
	if (mkov) delete[]mkov; if (stchsrver) delete[]stchsrver;
	etev = new double[cemv]; kttkv = new double[cemv]; mkov = new double[cemv]; stchsrver = new double[cemv];
	if ((!stchsrver) || (!etev) || (!kttkv) || (!mkov)) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < cemv; k++) { stchsrver[k] = 0.0; kttkv[k] = 0.0; mkov[k] = 0.0; etev[k] = 0.0; }
	etev[0] = tn; for (k = 1; k<cemv; k++) etev[k] = etev[k - 1] + detev;
	for (k = 0; k<cemv; k++) { g = epsisredver(etev[k], tkuscv, kuscv, dmkooscv, dkoscvt, dkoscvm, dkoscvl); stchsrver[k] = g; } 
	for (k=0; k<cemv; k++) cout << "tem = " << etev[k] << "\tst_ch = " << stchsrver[k] << "\t"; cout << endl;
	kttkv = opredKTPTverKarkVerm(etev, ektpv, porver, wsi, wal, wmg, vyfv, dmkov, tnoscv, dtoscv, dmkooscv, kuscv, tkuscv, stchsrver, cemv, vysv, vpkf);
	for (k=0; k<cemv; k++) cout << "tem = " << etev[k] << "\tkttkv = " << kttkv[k] << endl;
}
double *napMasEKTPVer(int vyfrve, int vysove, int vymeizvemafr, double *te, int v, double h, int vyukve, int nt, int wa)
{
	int k, nvm = dmkov, j; double F = 13.85e-4, nf84, *po = arrTemCold84(0), *tc84 = arrTemCold84(1);
	double *th84 = arrTemHigh84(), hf = 1e0, *tp84 = arrTepPot84(); nf84 = po[0];
	int n84 = 0; while (nf84>0.0) { nf84 = nf84 - hf; n84++; }
	double nf207, *th207 = arrTemHigh207(), *tc207 = arrTemCold207(1), *tp207 = arrTepPot207(); po = arrTemCold207(0); nf207 = po[0];
	int n207 = 0; while (nf207>0.0) { nf207 = nf207 - hf; n207++; } //nvm - длина массива коэффициентов
	double *ts207 = new double[n207], *ts84 = new double[n84];
	if ((!ts207) || (!ts84)) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < n84; k++) ts84[k] = (th84[k] + tc84[k]) / 2e0; for (k = 0; k < n207; k++) ts207[k] = (th207[k] + tc207[k]) / 2e0;
	double *koeft = NULL, *koefh = NULL, *koefc = NULL, *koefq = NULL, *koefs = NULL;
	double *ktp = NULL, *temvs = NULL, *temvh = NULL, *temvc = NULL, *tepv = NULL, *vm = NULL;
	//-------
	if ((!vyfrve) || (vyfrve==2)) { //фракция 2-0,7 мм или фракция 1,6-0,35 мм
		if (tc84) delete[]tc84; if (th84) delete[]th84; if (tp84) delete[]tp84;
		temvs = new double[n207]; temvh = new double[n207]; temvc = new double[n207]; tepv = new double[n207]; ktp = new double[n207];
		if ((!temvs) || (!temvh) || (!temvc) || (!tepv) || (!ktp)) { cout << snmv << endl; k = getchar(); exit(1); }
		for (k = 0; k < n207; k++) { temvs[k] = ts207[k]; temvh[k] = th207[k]; temvc[k] = tc207[k]; tepv[k] = tp207[k]; }
		for (k = 0; k < n207; k++) { ktp[k] = fabs(temvh[k] - temvc[k]) / h; tepv[k] = tepv[k] / F; ktp[k] = tepv[k] / ktp[k]; }
		if (!vysove) { //исходный
			koeft = oprkoefKTPiskhchao(vymeizvemafr, 0, temvs, tepv, h, temvh, temvc, ktp, nvm); koefh = oprkoefKTPiskhchao(vymeizvemafr, 1, temvs, tepv, h, temvh, temvc, ktp, nvm); koefc = oprkoefKTPiskhchao(vymeizvemafr, 2, temvs, tepv, h, temvh, temvc, ktp, nvm); koefs = oprkoefKTPiskhchao(vymeizvemafr, 3, temvs, tepv, h, temvh, temvc, ktp, nvm); koefq = oprkoefKTPiskhchao(vymeizvemafr, 4, temvs, tepv, h, temvh, temvc, ktp, nvm);
		}
		else if (vysove == 1) { //после повторных измерений
			koefq = danPoTemTepl2072(temvs, tepv, nvm); koeft = danPoTemTepl2072(temvs, ktp, nvm); koefh = danPoTemTepl2072(temvs, temvh, nvm); koefc = danPoTemTepl2072(temvs, temvc, nvm); koefs = danPoTemTepl2072(temvs, temvs, nvm);
		}
		else if (vysove == 2) { //после обжига при 1000 °С
			koefq = danPoTemTepl2073(temvs, tepv, nvm); koefh = danPoTemTepl2073(temvs, temvh, nvm); koefc = danPoTemTepl2073(temvs, temvc, nvm); koeft = danPoTemTepl2073(temvs, ktp, nvm); koefs = danPoTemTepl2073(temvs, temvs, nvm);
		}
		else if (vysove == 3) { //после повторного обжига при 1000 °С
			koefq = danPoTemTepl2074(temvs, tepv, nvm); koefh = danPoTemTepl2074(temvs, temvh, nvm); koefc = danPoTemTepl2074(temvs, temvc, nvm); koeft = danPoTemTepl2074(temvs, ktp, nvm); koefs = danPoTemTepl2074(temvs, temvs, nvm);
		}
	}
	else if (vyfrve == 1) { //фракция 8-4 мм
		if (tc207) delete[]tc207; if (th207) delete[]th207; if (tp207) delete[]tp207;
		temvs = new double[n84]; temvh = new double[n84]; temvc = new double[n84]; tepv = new double[n84]; ktp = new double[n84];
		if ((!temvs) || (!temvh) || (!temvc) || (!tepv) || (!ktp)) { cout << snmv << endl; k = getchar(); exit(1); }
		for (k = 0; k < n84; k++) { temvs[k] = ts84[k]; temvh[k] = th84[k]; temvc[k] = tc84[k]; tepv[k] = tp84[k]; tepv[k] = tepv[k] / F; ktp[k] = fabs(th84[k] - tc84[k]) / h; ktp[k] = tepv[k] / ktp[k]; }
		if (vyukve == 1) { //плоско-параллельная засыпка
			if (!vysove) { //исходный
				koefh = danPoTemTepl840(temvs, temvh, nvm); koefc = danPoTemTepl840(temvs, temvc, nvm); koeft = danPoTemTepl840(temvs, ktp, nvm); koefq = danPoTemTepl840(temvs, tepv, nvm); koefs = danPoTemTepl840(temvs, temvs, nvm);
			}
			else if (vysove == 1) { //после повторных измерений
				koefh = danPoTemTepl842(temvs, temvh, nvm); koefc = danPoTemTepl842(temvs, temvc, nvm); koeft = danPoTemTepl842(temvs, ktp, nvm); koefq = danPoTemTepl842(temvs, tepv, nvm); koefs = danPoTemTepl842(temvs, temvs, nvm);
			}
			else if (vysove == 2) { //после обжига
				koefh = danPoTemTepl845(temvs, temvh, nvm); koefc = danPoTemTepl845(temvs, temvc, nvm); koeft = danPoTemTepl845(temvs, ktp, nvm); koefq = danPoTemTepl845(temvs, tepv, nvm); koefs = danPoTemTepl845(temvs, temvs, nvm);
			}
		}
		else if (vyukve == 2) { //вертикальная засыпка
			if (!vysove) { //исходный
				koefh = danPoTemTepl841(temvs, temvh, nvm); koefc = danPoTemTepl841(temvs, temvc, nvm); koefq = danPoTemTepl841(temvs, tepv, nvm); koeft = danPoTemTepl841(temvs, ktp, nvm); koefs = danPoTemTepl841(temvs, temvs, nvm);
			}
			else if (vysove == 1) { //после повторных измерений
				koefh = danPoTemTepl844(temvs, temvh, nvm); koefc = danPoTemTepl844(temvs, temvc, nvm); koeft = danPoTemTepl844(temvs, ktp, nvm); koefq = danPoTemTepl844(temvs, tepv, nvm); koefs = danPoTemTepl844(temvs, temvs, nvm);
			}
			else if (vysove == 2) { //после обжига
				koefh = danPoTemTepl843(temvs, temvh, nvm); koefc = danPoTemTepl843(temvs, temvc, nvm); koeft = danPoTemTepl843(temvs, ktp, nvm); koefq = danPoTemTepl843(temvs, tepv, nvm); koefs = danPoTemTepl844(temvs, temvs, nvm);
			}
		}
	}
	if (!v) { vm = koefc; delete[]koefh; delete[]koefq; delete[]koeft; delete[]koefs; }
	else if (v == 1) { vm = koefh; delete[]koefc; delete[]koefq; delete[]koeft; delete[]koefs; }
	else if (v == 2) { vm = koefq; delete[]koefh; delete[]koefc; delete[]koeft; delete[]koefs; }
	else if (v == 3) { vm = koeft; delete[]koefc; delete[]koefh; delete[]koefq; delete[]koefs; }
	else if (v == 4) { vm = koefs; delete[]koefc; delete[]koefh; delete[]koefq; delete[]koeft; }
	double s, t, r, *ma = new double[nt]; if (!ma) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < nt; k++) { t = te[k]; s = 0.0; r = 0.0; 
	for (j = 0; j<nvm; j++) { s = s + vm[j] * pow(t, r); r = r + hf; } ma[k] = s; }
	if (po) delete[]po; if (temvs) delete[]temvs; if (temvh) delete[]temvh;
	if (temvc) delete[]temvc; if (tepv) delete[]tepv; if (ktp) delete[]ktp;
	if (!wa) { delete[]ma; return vm; }
	else if (wa == 1) { delete[]vm; return ma; }
}
double *oprkoefKTPiskhchao(int vmiv, int v, double *temvs, double *tepv, double h, double *temvh, double *temvc, double *ktp, int n) //n=4 - длина массива коэффициентов, vmiv - выбор метода измерений
{
	int nk = n, k, j, q, qn, nn = n; double hf = 1e0, nf207, *po, *koeft = NULL, *koefh = NULL, *koefc = NULL, *koefs = NULL, *koefq = NULL;
	double *ktp1 = NULL, *t10 = NULL, *t20 = NULL, *t30 = NULL, *t1 = NULL, *t2 = NULL, *t3 = NULL, *ktpn = NULL, *t = NULL, *mt = NULL, *ts = NULL;
	double *tsv = NULL, *tgv = NULL, *thv = NULL, *qov = NULL, *kt = NULL, *tgvn = NULL, *thvn = NULL, *qovn = NULL, *tsvn = NULL, *ktn = NULL;
	if (!vmiv) {  //установка Netzsch
		po = arrTem_VVF2(0);
		nk = 0; nf207 = po[nk]; while (nf207 > 0.0) { nf207 = nf207 - hf; nk++; }
		mt = arrTem_VVF2(1);
		ktpn = arrKTP_VVF2();
		tgv = opredTempHolGorVer(ktpn, mt, nk, nn, h, 0);
		thv = opredTempHolGorVer(ktpn, mt, nk, nn, h, 1);
		qov = opredTempHolGorVer(ktpn, mt, nk, nn, h, 2);
		tsv = opredTempHolGorVer(ktpn, mt, nk, nn, h, 3);
		thvn = vydelPol(thv, tgv, qov, ktpn, tsv, 0, nk);
		tgvn = vydelPol(thv, tgv, qov, ktpn, tsv, 1, nk);
		qovn = vydelPol(thv, tgv, qov, ktpn, tsv, 2, nk);
		tsvn = vydelPol(thv, tgv, qov, ktpn, tsv, 3, nk);
		ktn = vydelPol(thv, tgv, qov, ktpn, tsv, 4, nk);
		po = vydelPol(thv, tgv, qov, ktpn, tsv, 5, nk);
		nf207 = po[0]; q = 0; while (nf207 > 0.0) { nf207 = nf207 - hf; q++; } qn = q;
		koefh = new double[nn]; koefc = new double[nn]; koefq = new double[nn]; koeft = new double[nn]; koefs = new double[nn];
		if ((!koefh) || (!koefc) || (!koefq) || (!koeft) || (!koefs)) { cout << snmv << endl; k = getchar(); exit(1); }
		for (k = 0; k<nn; k++) { koefh[k] = 0.0; koefc[k] = 0.0; koefq[k] = 0.0; koeft[k] = 0.0; koefs[k] = 0.0; }
		koeft = koefPribSha(ktn, tsvn, qn, koeft); koefh = koefPribSha(tgvn, tsvn, qn, koefh); koefc = koefPribSha(thvn, tsvn, qn, koefc);
		koefs = koefPribSha(tsvn, tsvn, qn, koefs); koefq = koefPribSha(qovn, tsvn, qn, koefq);
		delete[]tsvn; delete[]tgvn; delete[]thvn; delete[]ktn; 
		delete[]qovn; delete[]po; delete[]tgv; delete[]thv; 
		delete[]qov; delete[]tsv; delete[]mt; delete[]ktpn;
	}
	else if (vmiv == 1) { //данные 2020 года
		ktp1 = arrKTP_VVF1(); t10 = arrTem1VVF1(); t20 = arrTem2VVF1(1); po = arrTem2VVF1(0); t30 = arrTem3VVF1();
		nf207 = po[0]; nk = 0; while (nf207>0.0) { nf207 = nf207 - hf; nk++; }
		ktpn = danIskh207(ktp1, 0, nk); t1 = danIskh207(t10, 1, nk); t2 = danIskh207(t20, 2, nk); t3 = danIskh207(t30, 3, nk);
		po = NULL; po = danIskh207(po, 0, 0); nf207 = po[0]; nk = 0; while (nf207 > 0.0) { nf207 = nf207 - hf; nk++; }
		delete[]ktp1; delete[]t10; delete[]t20; delete[]t30; delete[]po;
		ts = new double[nk]; if (!ts) { cout << snmv << endl; k = getchar(); exit(1); } //nk - размер массива ts
		for (k = 0; k < nk; k++) ts[k] = (t1[k] + t2[k]) / 2e0;
		koeft = new double[n]; koefh = new double[n]; koefc = new double[n]; koefs = new double[n]; koefq = new double[n];
		if ((!koeft) || (!koefh) || (!koefc) || (!koefs) || (!koefq)) { cout << snmv << endl; k = getchar(); exit(1); }
		for (k = 0; k < n; k++) { koefh[k] = 0.0; koefc[k] = 0.0; koefq[k] = 0.0; koeft[k] = 0.0; koefs[k] = 0.0; }
		koeft = koefPribSha(ktpn, ts, nk, koeft); koefh = koefPribSha(t1, ts, nk, koefh); koefc = koefPribSha(t2, ts, nk, koefc); koefs = koefPribSha(t3, ts, nk, koefs);
		double s, r, *qv = new double[n]; if (!qv) { cout << snmv << endl; k = getchar(); exit(1); }
		for (k = 0; k < nk; k++) {
			s = 0.0; r = 0.0;
			for (j = 0; j < n; j++) { s = s + koeft[j] * pow(ts[j], r); r = r + hf; }
			qv[k] = s*fabs(t1[k] - t2[k]) / h;
		}
		koefq = koefPribSha(qv, ts, nk, koefq);
		delete[]qv; delete[]t1; delete[]t2; delete[]t3; delete[]ktpn;
	}
	else if (vmiv == 2) { //данные 2019 года
		koefq = danPoTemTepl2071(temvs, tepv, n);
		koefh = danPoTemTepl2071(temvs, temvh, n);
		koefc = danPoTemTepl2071(temvs, temvc, n);
		koeft = danPoTemTepl2071(temvs, ktp, n);
		koefs = danPoTemTepl2071(temvs, temvs, n);
	}
	if (!v) { delete[]koefh; delete[]koefc; delete[]koefs; delete[]koefq; t = koeft; }
	else if (v == 1) { delete[]koeft; delete[]koefc; delete[]koefs; delete[]koefq; t = koefh; }
	else if (v == 2) { delete[]koefh; delete[]koeft; delete[]koefs; delete[]koefq; t = koefc; }
	else if (v == 3) { delete[]koefh; delete[]koefc; delete[]koeft; delete[]koefq; t = koefs; }
	else if (v == 4) { delete[]koefh; delete[]koefc; delete[]koefs; delete[]koeft; t = koefq; }
	return t;
}
double *danIskh207(double *ma, int v, int n)
{
	if (ma) {
		int k, np = 3; double s = 0.0, r = 0.0, hf = 1e0, *m = new double[n];
		if (!m) { cout << snmv << endl; k = getchar(); exit(1); }
		for (k = 0; k < np; k++) { s = s + ma[k]; r = r + hf; } s = s / r;
		m[0] = s; m[1] = ma[3]; m[2] = (ma[4] + ma[5]) / 2e0;
		return m;
	}
	else { int k; double *po = new double[1]; 
	if (!po) { cout << snmv << endl; k = getchar(); exit(1); } po[0] = 3e0; return po; }
}
//---------
double *arrTem_VVF2(int v)
{
	int k = 0, le = 8; double *tem = new double[le], hf = 1e0, nf = 0.0, *po = new double[1];
	if ((!tem) || (!po)) { cout << snmv << endl; k = getchar(); exit(1); }
	tem[k] = 27.967;  k++; tem[k] = 93.833; k++; tem[k] = 192.5;   k++; tem[k] = 341.667; k++;
	tem[k] = 491.467; k++; tem[k] = 641.2;  k++; tem[k] = 790.933; k++; tem[k] = 991.133; k++;
	nf = 0.0; for (k = 0; k < le; k++) { tem[k] = tem[k] + tev0; nf = nf + hf; } po[0] = nf;
	if (!v) { delete[]tem; return po; }
	else if (v == 1) { delete[]po; return tem; }
}
double *arrKTP_VVF2()
{
	int k = 0, le = 8; double *ktp = new double[le]; if (!ktp) { cout << snmv << endl; k = getchar(); exit(1); }
	ktp[k] = 7e-2;   k++; ktp[k] = 8e-2;   k++; ktp[k] = 103e-3; k++; ktp[k] = 15e-2; k++;
	ktp[k] = 207e-3; k++; ktp[k] = 283e-3; k++; ktp[k] = 373e-3; k++; ktp[k] = 477e-3; return ktp;
}
double *arrKTP_VVF1()
{
	int k = 0, n = 6; double *a = new double[n]; if (!a) { cout << snmv << endl; k = getchar(); exit(1); }
	a[k] = 0.175566644058715; k++; a[k] = 0.176801537812368; k++; a[k] = 0.179324717653617; k++;
	a[k] = 0.211768953068592; k++; a[k] = 0.237194543621728; k++; a[k] = 0.237231989760775;
	return a;
}
double *arrTem1VVF1()
{
	int k = 0, n = 6; double *a = new double[n]; if (!a) { cout << snmv << endl; k = getchar(); exit(1); }
	a[k] = 585.0; k++; a[k] = 600.0; k++; a[k] = 585.0; k++; a[k] = 800.0; k++; a[k] = 1000.0; k++; a[k] = 1000.0;
	for (k = 0; k < n; k++) a[k] = a[k] + tev0;
	return a;
}
double *arrTem2VVF1(int v)
{
	int k = 0, n = 6; double *a = new double[n], nf = 0.0, hf = 1e0, *po = new double[1];
	if ((!a) || (!po)) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < n; k++) nf = nf + hf; po[0] = nf;
	k = 0; a[k] = 119.75; k++; a[k] = 138.0; k++; a[k] = 129.5; k++; a[k] = 200.0; k++; a[k] = 273.0; k++; a[k] = 261.0;
	for (k = 0; k < n; k++) a[k] = a[k] + tev0;
	if (!v) { delete[]a; return po; }
	else { delete[]po; return a; }
}
double *arrTem3VVF1()
{
	int k = 0, n = 6; double *a = new double[n]; if (!a) { cout << snmv << endl; k = getchar(); exit(1); }
	a[k] = 377.0; k++; a[k] = 396.0; k++; a[k] = 383.5; k++; a[k] = 548.0; k++; a[k] = 703.0; k++; a[k] = 697.25;
	for (k = 0; k < n; k++) a[k] = a[k] + tev0;
	return a;
}
//---------
double *danPoTemTepl840(double *temvs, double *temvh, int n) //Засыпка плоско-параллельная, исходный
{
	double tho1 = 0.0, tho2 = 0.0, tem1 = 0.0, tem2 = 0.0, *isdan = new double[n]; int k;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < n; k++) isdan[k] = 0.0;
	tem1 = temvs[14] + temvs[16];
	tem2 = temvs[15] + temvs[17];
	tho1 = (temvh[14] * temvs[14] + temvh[16] * temvs[16]) / tem1;
	tho2 = (temvh[15] * temvs[15] + temvh[17] * temvs[17]) / tem2;
	tem1 = tem1 / 2e0; tem2 = tem2 / 2e0;
	double k1; if (tem2 != tem1) k1 = (tho2 - tho1) / (tem2 - tem1); else k1 = 0.0; double k2 = tho2 - k1*tem2;
	isdan[0] = k2; isdan[1] = k1; return isdan;
}
double *danPoTemTepl841(double *temvs, double *temvh, int n) //Засыпка вертикальная, исходный
{
	double tc1 = 0.0, tc2 = 0.0, tem1 = 0.0, tem2 = 0.0, *isdan = new double[n]; int k;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k<n; k++) isdan[k] = 0.0;
	tem1 = temvs[4] + temvs[6] + temvs[10];
	tem2 = temvs[5] + temvs[7] + temvs[11];
	if (fabs(tem1)>0.0) tc1 = (temvh[4] * temvs[4] + temvh[6] * temvs[6] + temvh[10] * temvs[10]) / tem1;
	if (fabs(tem2) > 0.0) tc2 = (temvh[5] * temvs[5] + temvh[7] * temvs[7] + temvh[11] * temvs[11]) / tem2;
	tem1 = tem1 / 3e0; tem2 = tem2 / 3e0;
	double k1; if (tem2 != tem1) k1 = (tc2 - tc1) / (tem2 - tem1); else k1 = 0.0; double k2 = tc2 - k1*tem2;
	isdan[0] = k2; isdan[1] = k1; return isdan;
}
double *danPoTemTepl842(double *temvs, double *temvh, int n) //Засыпка плоско-параллельная, повторы
{
	double tc1 = 0.0, tc2 = 0.0, tem1 = 0.0, tem2 = 0.0, *isdan = new double[n]; int k;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k<n; k++) isdan[k] = 0.0;
	tem1 = temvs[18] + temvs[20] + temvs[22];
	tem2 = temvs[19] + temvs[21] + temvs[23];
	if (fabs(tem1)>0.0) tc1 = (temvh[18] * temvs[18] + temvh[20] * temvs[20] + temvh[22] * temvs[22]) / tem1; else tc1 = 0.0;
	if (fabs(tem2) > 0.0) tc2 = (temvh[19] * temvs[19] + temvh[21] * temvs[21] + temvh[23] * temvs[23]) / tem2; else tc2 = 0.0;
	tem1 = tem1 / 3e0; tem2 = tem2 / 3e0;
	double k1; if (tem2 != tem1) k1 = (tc2 - tc1) / (tem2 - tem1); else k1 = 0.0; double k2 = tc2 - k1*tem2;
	isdan[0] = k2; isdan[1] = k1; return isdan;
}
double *danPoTemTepl843(double *temvs, double *temvh, int n) //Засыпка вертикальная, после обжига при 1000 °С
{
	double tc1 = 0.0, tc2 = 0.0, tem1 = 0.0, tem2 = 0.0, *isdan = new double[n]; int k;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k<n; k++) isdan[k] = 0.0;
	tem1 = temvs[0] + temvs[2];
	tem2 = temvs[1] + temvs[3];
	if (fabs(tem1)>0.0) tc1 = (temvh[0] * temvs[0] + temvh[2] * temvs[2]) / tem1;
	if (fabs(tem2) > 0.0) tc2 = (temvh[1] * temvs[1] + temvh[3] * temvs[3]) / tem2;
	tem1 = tem1 / 2e0; tem2 = tem2 / 2e0;
	double k1; if (tem2 != tem1) k1 = (tc2 - tc1) / (tem2 - tem1); else k1 = 0.0; double k2 = tc2 - k1*tem2;
	isdan[0] = k2; isdan[1] = k1; return isdan;
}
double *danPoTemTepl844(double *temvs, double *temvh, int n) //Засыпка вертикальная, повторы
{
	double tc1 = 0.0, tc2 = 0.0, tem1 = 0.0, tem2 = 0.0, *isdan = new double[n]; int k;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k<n; k++) isdan[k] = 0.0;
	tem1 = temvs[8] + temvs[12];
	tem2 = temvs[9] + temvs[13];
	if (fabs(tem1)>0.0) tc1 = (temvh[8] * temvs[8] + temvh[12] * temvs[12]) / tem1;
	if (fabs(tem2) > 0.0) tc2 = (temvh[9] * temvs[9] + temvh[13] * temvs[13]) / tem2;
	tem1 = tem1 / 2e0; tem2 = tem2 / 2e0;
	double k1; if (tem2 != tem1) k1 = (tc2 - tc1) / (tem2 - tem1); else k1 = 0.0; double k2 = tc2 - k1*tem2;
	isdan[0] = k2; isdan[1] = k1; return isdan;
}
double *danPoTemTepl845(double *temvs, double *temvh, int n) //Засыпка плоско-параллельная, после обжига при 1000 °С
{
	double tc1 = 0.0, tc2 = 0.0, tem1 = 0.0, tem2 = 0.0, *isdan = new double[n]; int k;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k<n; k++) isdan[k] = 0.0;
	tem1 = temvs[24];
	tem2 = temvs[25];
	if (fabs(tem1)>0.0) tc1 = (temvh[24] * temvs[24]) / tem1;
	if (fabs(tem2) > 0.0) tc2 = (temvh[25] * temvs[25]) / tem2;
	tem1 = tem1 / 1e0; tem2 = tem2 / 1e0;
	double k1; if (tem2 != tem1) k1 = (tc2 - tc1) / (tem2 - tem1); else k1 = 0.0; double k2 = tc2 - k1*tem2;
	isdan[0] = k2; isdan[1] = k1; return isdan;
}
//---------
double *danPoTemTepl2071(double *temvs, double *tepv, int n) //Засыпка исходная, фракция 2-0,7 мм
{
	double tepv1 = 0.0, tepv2 = 0.0, tem1 = 0.0, tem2 = 0.0, *isdan = new double[n]; int k;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k<n; k++) isdan[k] = 0.0;
	tem1 = temvs[0];
	tem2 = temvs[1];
	if (fabs(tem1)>0.0) tepv1 = (tepv[0] * temvs[0]) / tem1;
	if (fabs(tem2) > 0.0) tepv2 = (tepv[1] * temvs[1]) / tem2;
	tem1 = tem1 / 1e0; tem2 = tem2 / 1e0;
	double k1; if (tem2 != tem1) k1 = (tepv2 - tepv1) / (tem2 - tem1); else k1 = 0.0; double k2 = tepv2 - k1*tem2;
	isdan[0] = k2; isdan[1] = k1; return isdan;
}
double *danPoTemTepl2072(double *temvs, double *tepv, int n) //Фракция 2-0,7 мм (повторные измерения)
{
	double tepv1 = 0.0, tepv2 = 0.0, tem1 = 0.0, tem2 = 0.0, *isdan = new double[n]; int k;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k<n; k++) isdan[k] = 0.0;
	tem1 = temvs[2] + temvs[4];
	tem2 = temvs[3] + temvs[5];
	if (fabs(tem1)>0.0) tepv1 = (tepv[2] * temvs[2] + tepv[4] * temvs[4]) / tem1;
	if (fabs(tem2) > 0.0) tepv2 = (tepv[3] * temvs[3] + tepv[5] * temvs[5]) / tem2;
	tem1 = tem1 / 2e0; tem2 = tem2 / 2e0;
	double k1; if (tem2 != tem1) k1 = (tepv2 - tepv1) / (tem2 - tem1); else k1 = 0.0; double k2 = tepv2 - k1*tem2;
	isdan[0] = k2; isdan[1] = k1; return isdan;
}
double *danPoTemTepl2073(double *temvs, double *tepv, int n) //Фракция 2-0,7 мм, после обжига при 1000 °С
{
	double tepv1 = 0.0, tepv2 = 0.0, tem1 = 0.0, tem2 = 0.0, *isdan = new double[n]; int k;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < n; k++) isdan[k] = 0.0;
	tem1 = temvs[6] + temvs[8];
	tem2 = temvs[7] + temvs[9];
	tepv1 = (tepv[6] * temvs[6] + tepv[8] * temvs[8]) / tem1;
	tepv2 = (tepv[7] * temvs[7] + tepv[9] * temvs[9]) / tem2;
	tem1 = tem1 / 2e0; tem2 = tem2 / 2e0;
	double k1; if (tem2 != tem1) k1 = (tepv2 - tepv1) / (tem2 - tem1); else k1 = 0.0; double k2 = tepv2 - k1*tem2;
	isdan[0] = k2; isdan[1] = k1; return isdan;
}
double *danPoTemTepl2074(double *temvs, double *tepv, int n) //Фракция 2-0,7 мм, после повторного обжига при 1000 °С
{
	double tepv1 = 0.0, tepv2 = 0.0, tem1 = 0.0, tem2 = 0.0, *isdan = new double[n]; int k;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < n; k++) isdan[k] = 0.0;
	tem1 = temvs[2];
	tem2 = temvs[3];
	tepv1 = (tepv[2] * temvs[2]) / tem1;
	tepv2 = (tepv[3] * temvs[3]) / tem2;
	tem1 = tem1 / 1e0; tem2 = tem2 / 1e0;
	double k1; if (tem2 != tem1) k1 = (tepv2 - tepv1) / (tem2 - tem1); else k1 = 0.0; double k2 = tepv2 - k1*tem2;
	isdan[0] = k2; isdan[1] = k1; return isdan;
}
//---------
double *arrTemCold84(int v)
{
	int k = 0, n = 26; double *temcol = new double[n], *po = new double[1], nf = 0.0, hf = 1e0;
	if ((!temcol) || (!po)) { cout << snmv << endl; k = getchar(); exit(1); }
	temcol[k] = 168.0;  k++; temcol[k] = 369.0;  k++; temcol[k] = 168.0;  k++; temcol[k] = 369.0; k++;
	temcol[k] = 148.6;  k++; temcol[k] = 356.5;  k++; temcol[k] = 184.0;  k++; temcol[k] = 396.0; k++;
	temcol[k] = 148.75; k++; temcol[k] = 350.0;  k++; temcol[k] = 166.0;  k++; temcol[k] = 375.0; k++;
	temcol[k] = 171.0;  k++; temcol[k] = 383.5;  k++; temcol[k] = 106.75; k++; temcol[k] = 242.0; k++;
	temcol[k] = 123.0;  k++; temcol[k] = 294.0;  k++; temcol[k] = 111.0;  k++; temcol[k] = 240.0; k++;
	temcol[k] = 109.0;  k++; temcol[k] = 232.75; k++; temcol[k] = 127.0;  k++; temcol[k] = 291.0; k++;
	temcol[k] = 221.0;  k++; temcol[k] = 443.75;
	for (k = 0; k < n; k++) temcol[k] = temcol[k] + tev0;
	for (k = 0; k < n; k++) nf = nf + hf; po[0] = nf;
	if (!v) { delete[]temcol; return po; }
	else if (v == 1) { delete[]po; return temcol; }
}
double *arrTemHigh84()
{
	int k = 0, n = 26; double *temh = new double[n];
	if (!temh) { cout << snmv << endl; k = getchar(); exit(1); }
	temh[k] = 545.0; k++; temh[k] = 927.0;  k++; temh[k] = 530.0; k++; temh[k] = 925.0;  k++;
	temh[k] = 560.0; k++; temh[k] = 950.0;  k++; temh[k] = 560.0; k++; temh[k] = 900.0;  k++;
	temh[k] = 558.0; k++; temh[k] = 950.0;  k++; temh[k] = 540.0; k++; temh[k] = 920.0;  k++;
	temh[k] = 540.0; k++; temh[k] = 920.0;  k++; temh[k] = 600.0; k++; temh[k] = 1000.0; k++;
	temh[k] = 580.0; k++; temh[k] = 1000.0; k++; temh[k] = 587.0; k++; temh[k] = 1000.0; k++;
	temh[k] = 590.0; k++; temh[k] = 1000.0; k++; temh[k] = 580.0; k++; temh[k] = 1000.0; k++;
	temh[k] = 483.0; k++; temh[k] = 850.0;
	for (k = 0; k < n; k++) temh[k] = temh[k] + tev0;
	return temh;
}
double *arrTepPot84()
{
	int k = 0, n = 26; double *tep = new double[n];
	if (!tep) { cout << snmv << endl; k = getchar(); exit(1); }
	tep[k] = 3.9805;  k++; tep[k] = 9.5532;   k++; tep[k] = 3.6447; k++; tep[k] = 9.4779;  k++;
	tep[k] = 3.55732; k++; tep[k] = 11.4997;  k++; tep[k] = 4.4624; k++; tep[k] = 11.6021; k++;
	tep[k] = 4.7766;  k++; tep[k] = 11.5016;  k++; tep[k] = 3.4023; k++; tep[k] = 10.2068; k++;
	tep[k] = 3.92812; k++; tep[k] = 11.17333; k++; tep[k] = 3.5144; k++; tep[k] = 8.6593;  k++;
	tep[k] = 2.977;   k++; tep[k] = 8.0448;   k++; tep[k] = 3.352;  k++; tep[k] = 9.218;   k++;
	tep[k] = 3.0313;  k++; tep[k] = 7.7946;   k++; tep[k] = 3.0671; k++; tep[k] = 6.1342;  k++;
	tep[k] = 1.73466; k++; tep[k] = 4.329670;
	return tep;
}
double *arrTemCold207(int v)
{
	int k = 0, n = 10; double *temcol = new double[n], *po = new double[1], nf = 0.0, hf = 1e0;
	if ((!temcol) || (!po)) { cout << snmv << endl; k = getchar(); exit(1); }
	temcol[k] = 109.0;  k++; temcol[k] = 235.0; k++; temcol[k] = 101.0; k++; temcol[k] = 199.0; k++;
	temcol[k] = 108.75; k++; temcol[k] = 238.0; k++; temcol[k] = 124.0; k++; temcol[k] = 266.0; k++;
	temcol[k] = 111.0;  k++; temcol[k] = 262.0;
	for (k = 0; k < n; k++) temcol[k] = temcol[k] + tev0;
	for (k = 0; k < n; k++) nf = nf + hf; po[0] = nf;
	if (!v) { delete[]temcol; return po; }
	else if (v == 1) { delete[]po; return temcol; }
}
double *arrTemHigh207()
{
	int k = 0, n = 10; double *temh = new double[n];
	if (!temh) { cout << snmv << endl; k = getchar(); exit(1); }
	temh[k] = 585.0; k++; temh[k] = 1000.0; k++; temh[k] = 603.0; k++; temh[k] = 1000.0;  k++;
	temh[k] = 603.0; k++; temh[k] = 1000.0; k++; temh[k] = 571.5; k++; temh[k] = 1000.75; k++;
	temh[k] = 583.0; k++; temh[k] = 1000.0;
	for (k = 0; k < n; k++) temh[k] = temh[k] + tev0;
	return temh;
}
double *arrTepPot207()
{
	int k = 0, n = 10; double *tepot = new double[n];
	if (!tepot) { cout << snmv << endl; k = getchar(); exit(1); }
	tepot[k] = 3.9596;    k++; tepot[k] = 8.6377; k++; tepot[k] = 2.3003;  k++; tepot[k] = 5.3674; k++;
	tepot[k] = 3.5614925; k++; tepot[k] = 7.123;  k++; tepot[k] = 2.12992; k++; tepot[k] = 7.6956; k++;
	tepot[k] = 2.3003;    k++; tepot[k] = 6.9009;
	return tepot;
}
double *vydelPol(double *temvcs, double *temvhs, double *qos, double *ektpvs, double *temvss, int v, int n)
{
	int q, qn, k;
	q = 0; for (k = 0; k<n; k++) if ((temvcs[k]>0.0) && (temvhs[k] > 0) && (qos[k] > 0.0) && (ektpvs[k] > 0.0) && (temvss[k] > 0.0) && (temvhs[k]<templa)) q++; qn = q;
	double hf = 1e0, nf = 0.0, *p = new double[1], *vm; for (k = 0; k<qn; k++) nf = nf + hf;
	vm = new double[qn]; if ((!vm) || (!p)) { cout << snmv << endl; k = getchar(); exit(1); } p[0] = nf;
	q = 0; for (k = 0; k<n; k++) {
		if ((temvcs[k]>0.0) && (temvhs[k]>0) && (qos[k]>0.0) && (ektpvs[k] > 0.0) && (temvss[k] > 0.0) && (temvhs[k] < templa)) {
			if (!v) vm[q] = temvcs[k];
			else if (v == 1) vm[q] = temvhs[k];
			else if (v == 2) vm[q] = qos[k];
			else if (v == 3) vm[q] = temvss[k];
			else if (v == 4) vm[q] = ektpvs[k];
			q++;
		}
	}
	if (v == 5) { delete[]vm; vm = p; }
	return vm;
}