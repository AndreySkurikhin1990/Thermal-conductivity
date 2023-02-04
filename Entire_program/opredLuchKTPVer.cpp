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
#define vmivmf 0 //âûáîð ìåòîäà èçìåðåíèé äëÿ ôðàêöèè 2-0,7 ìì: 0 - íåñòàöèîíàðíûé, 1 - ñòàöèîíàðíûé
#define vyfv 0 //âûáîð ôðàêöèè: 0 - ôðàêöèÿ 2-0,7 ìì, 1 - ôðàêöèÿ 8-4 ìì, 2 - ôðàêöèÿ 1,6-0,35 ìì
#define vyuv 1 //âûáîð óêëàäêè: 1 - ïëîñêîïàðàëëåëüíàÿ, 2 - âåðòèêàëüíàÿ
#define vysv 0 //âûáîð ñîñòîÿíèÿ: 0 - èñõîäíîå, 1 - ïîñëå ïîâòîðíûõ èçìåðåíèé, 2 - ïîñëå ïðîêàëèâàíèÿ ïðè 1000 ãðàä Ñ
using namespace std;
const int dsv = 50, dmkvover = 28, dmkooscv = 14, N = 13, ks = 10, dmkov = 4, qtn = 8, nxtv = 30, nnxtv = 0;
const int vtvn = 0, isrp=0, vpkf=0, vpmf=0, cemdu=6, cemdumi=2; //0 - ñòàðûå, 1 - íîâûå
const double detev = 1e2, tev0 = 273.15, nnxtfv = 0.0, ssiv84 = 467.0*1e-3, salv84 = 129.0*1e-3, smgv84 = 282.0*1e-3;
const double ssiv207 = 458.0*1e-3, salv207 = 14.0*1e-2, smgv207 = 29.0*1e-2;
const double tnoscv = 3e2, tnacv = 2e2; //dkoscv - ìàêñèìàëüíîå îòêëîíåíèå îò ýêñïåðèìåíòàëüíûõ äàííûõ â ëèòåðàòóðíûõ èñòî÷íèêàõ
const double dtoscv = 1e2, tocrasver = 1e-8, por207 = 66.35*1e-2, poris84 = 55.75*1e-2, porin84=81.53*1e-2, templa = 134.0*1e1 + tev0;
const double epsi = 1e-15, y0ver = 3e1*1e-3, pksvv = 0.0, dkoalv = 0.444905, poro84=86.61*1e-2; //dkoalv - îïðåäåëåíèå ÊÏ èç ÊÎ
const double poro16035=84.36*1e-2, por16035=83.97*1e-2;
const char nfv = 'A';
struct derevo {
	int otre; //1 - îòðàæåíèå èëè 0 - ïðîïóñêàíèå
	int ste; //íîìåð ñòåíêè
	int vis; //âèäèìîñòü: 1 - âèäåí, 0 - íåò
	int lev; //íîìåð óðîâíÿ
	struct derevo *back; //óêàçàòåëü íàçàä
	struct derevo *next; //óêàçàòåëü âïåðåä
};
class KoePog { public: double alp, tem; KoePog *nex; };
//----------------------------------------
double *alphaVer = NULL, *Tver = NULL, *Aver = NULL, *Rver = NULL, hk = 0.0, *Raver = NULL, *Taver = NULL, *Aaver = NULL, *Rtver = NULL, *Ttver = NULL, *Atver = NULL;
double *mgover = NULL, *siover = NULL, *alover = NULL; 
int vtv = 0, sctxv = 0, cel = 4, dkoscvl = 6, cemv = 11, vtvk=cemv; /*âûáîð òåìïåðàòóðû âåðìèêóëèòà*/
char *snv = NULL, *ssv = NULL, *skptv = NULL, *svsv = NULL, *snmv = NULL, *sfnov = NULL;
char *sfatv = NULL, *sfov = NULL, *svfdv = NULL, *svfdvu = NULL;
double *qobv = NULL, *etev = NULL, *txver = NULL, *qxver = NULL, *lxver = NULL, *kxver = NULL;
double *dtxver = NULL, *stchsrver = NULL, *tgorv = NULL, *tholv = NULL;
double *ktrv = NULL, *dkoscvm = NULL, *dkoscvt = NULL, hver = 0.0, hvove = 0.0, hkver = 0.0, qobver = 0.0, *ktpvover;
double *vtever, porver = 0.0, *ecktpver, *tsredver, dkospv = 1.3002913762; //dkospv - äîïîëíèòåëüíûé êîýôôèöèåíò äëÿ ÊÎ ñ ó÷åòîì ïîðèñòîé ñòðóêòóðû;
double *kttkv = NULL, cov = 0.0, bov = 0.0, *mkov = NULL, *ektpv = NULL, *kektpv = NULL, *mtsv = NULL, tmav, tmiv, dpctv = 1e0;
double *Tpctv, tmiver, tmaver, tnav = tnacv + tev0, *tkuscv, *kuscv, cnv = 1e-2, ckv = 1e2, bnv = 1e-2;
double bkv = 1e2, *temrasv = NULL, tnrv, *ooxver = NULL, *dpctvm, nxtfv, ssiv = 368.0*1e-3, salv = 132.0*1e-3, smgv = 217.0*1e-3;
//----------------------------------------
void zadrktVerNac();
double *podschchieleSha(int, int, int, double *, double *);
double *izstNVer(int, int, int, int);
double rasotprpovsSha(int *, int *, int, int, int, int, double *, double *);
double otrprovsSha(int, int, double *, double *);
double *kopoVer(double, double, int, double, int, int, double, double, int, double *);
double opredLuchSostVer(double *, double *, int, int, int);
double *koefPribSha(double *, double *, int, double *);
double *zadrktVer(int, int, double, int, double, double, int, int, double, int, int, double *);
double *zadrktSha(int, int, double, int, double, double, int, double, int, int, double *);
double *zadrktitom(int, int, double, int, double, double, int, int, double, int, int, double *);
double *zadrktkvi(int, int, double, int, double, double, int, int, double, int, int, double *);
double **chaRTAVer(int, double *, double *, double *, double *, double *, double *, int);
double **izmRTAVer(double *, int, int, double *, double *, double *, double *, double *, double *, int);
double *opredKoefOtr(double *, double *, int, double, int, int, double *, double *, int);
void napstrver();
double opredKTPTKTochSha(double *, double *, double, int);
double *opredKTPTverKarkVerm(double *, double *, double, double, double, double, int, int, double, double, int, double *, double *, double *, int, int, int, int);
double *reshnewtrafs(double *, double *, double *, double *, int, double *, double *, int, double *, double *, int, double *, double *, int, double, double *, double *, int, double, double, double, double *, int, double, double, double);
void initarrver(int, double, double, double);
void NapMasVozdSha(double*, double *, int);
void vyvodfile(double *, int, int, double, char *);
void novNapMas(double);
void zapisvfile(double *, int, char *);
double BolTochRasAlpha(int, int, double, double, double, double *, double *, double *, int, int, char *, double *, double *, int, double *, double *, int, double *, double *, double *, double *, double *, double *);
void osvpamver();
double **RaschRTAVer(int, double, double, double, int, int, double, int, double, int, int, int);
double **RaschRTAkvi(int, double, double, double, int, int, double, int, double, int, int, int);
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
double *EffectTols(double *, double *, double *, double, double, int);
double RasFracXeffVer60(int);
double RasFracXeffVer60_100(int);
double RasFracXeffVer100_150(int);
double RasFracXeffVer150_200(int);
double *koefoslab(double, double, double, double *, int, double *);
double ReflSredVer(double);
double ReflSredkvi(double, int);
double ReflSredSha(double);
double ReflSreditom(double, int);
double reflsha(double);
double reflver(double);
double oprProcSoderpoPoris(double *, double *, double, int);
double SeryeStenkiRasIzl(double, double, double, double *, double *, double *, double *, double *, int *);
double F0_lamT(double);
void opredtemphcVer(double *, double *, double *, double *, double *, int, int, double);
double *oprRasTemNachVer(int, int, double *, double *, int, double *, double, double, double, double, double, int);
double *oprEffDoliTepPerenVer(double, double, double);
double *KorrZnachVozdProsVer(double, double, double, int);
double *opredPolTempTvKarShaFragm(double *, int, double *, double *, int, double, double, double, double, double);
double RasIzlSerStenNac(double *, double *, double *, double *, double *, double *, double *, double, double, double, double, int, double, double, double, double, int, double *, double *, int, int, double *, double *, double *, char *, int, double *, double *, int);
double *FuncRaschIzl(double, double, double, double, int, double *, double, double *, double, double, double, double, int, double *, double *, int, int, char *, double, int, double *, double, double *, double *, int, double *, double *, int, double *, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, int, double, double, int, double *, double *, double, double, double, char *, int);
double KorrZnachVozdProsVermik(double, double, double);
double **vydelPol(double *, double *, double *, double *, double *, int, int, double **, int);
double **RaschRTASha(int, double, double, double, double, double, int, int, int, int);
double **opredTempHolGorVer(double *, double *, int, int, double, int, double **, int, double *, double *);
//----------------
double **oprkoefKTPiskhchao(int, int, double *, double, int, double **, int, int, int);
double **napMasEKTPVer(int, int, int, double *, int, double, int, int, int, double **, int);
double *arrTem_VVF2(int); //ìàññèâ òåìïåðàòóð - ýêñïåðèìåíòàëüíûå äàííûå íà Netzsch - õàîòè÷íàÿ çàñûïêà ôðàêöèè 2-0,7 ìì (èñõîäíûé)
double *arrKTP_VVF2(); //ìàññèâ ÊÒÏ âåðìèêóëèòà
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
double **PoiskZavVelTemVer(int, double **, int, double *, double *, int);
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
	snv[k] = 'À'; k++; snv[k] = 'ñ';  k++; snv[k] = 'ï';  k++; snv[k] = 'è';  k++; snv[k] = 'ð'; k++;
	snv[k] = 'à'; k++; snv[k] = 'í';  k++; snv[k] = 'ò';  k++; snv[k] = 'ó';  k++; snv[k] = 'ð'; k++;
	snv[k] = 'à'; k++; snv[k] = '\\'; k++; snv[k] = '\\'; k++; snv[k] = 't';  k++; snv[k] = 'm'; k++;
	snv[k] = 'p'; k++; snv[k] = '\\'; k++; snv[k] = '\\'; k++; snv[k] = '\0';
	k = 0;
	ssv[k] = 'C';  k++; ssv[k] = ':'; k++; ssv[k] = '\\'; k++; ssv[k] = '\\'; k++; ssv[k] = 'U';  k++;
	ssv[k] = 's';  k++; ssv[k] = 'e'; k++; ssv[k] = 'r';  k++; ssv[k] = 's';  k++; ssv[k] = '\\'; k++;
	ssv[k] = '\\'; k++; ssv[k] = 'À'; k++; ssv[k] = 'í';  k++; ssv[k] = 'ä';  k++; ssv[k] = 'ð';  k++;
	ssv[k] = 'å';  k++; ssv[k] = 'é'; k++; ssv[k] = '\\'; k++; ssv[k] = '\\'; k++; ssv[k] = 'D';  k++;
	ssv[k] = 'o';  k++; ssv[k] = 'c'; k++; ssv[k] = 'u';  k++; ssv[k] = 'm';  k++; ssv[k] = 'e';  k++;
	ssv[k] = 'n';  k++; ssv[k] = 't'; k++; ssv[k] = 's';  k++; ssv[k] = '\\'; k++; ssv[k] = '\\'; k++;
	ssv[k] = '_';  k++; ssv[k] = 'À'; k++; ssv[k] = 'ñ';  k++; ssv[k] = 'ï';  k++; ssv[k] = 'è';  k++;
	ssv[k] = 'ð';  k++; ssv[k] = 'à'; k++; ssv[k] = 'í';  k++; ssv[k] = 'ò';  k++; ssv[k] = 'ó';  k++;
	ssv[k] = 'ð';  k++; ssv[k] = 'à'; k++; ssv[k] = '\\'; k++; ssv[k] = '\\'; k++; ssv[k] = 't';  k++;
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
	double tnd = 6e2, dtd = 2e2, tm, hf=1e0; k=0; dkoscvt[k] = tnd; for (k = 1; k < dkoscvl; k++) dkoscvt[k] = dkoscvt[k - 1] + dtd;
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
	for (j=0; j<dmkvover; j++) { ktpvover[j]=0.0; vtever[j]=0.0; } if ((!ktpvover) || (!vtever)) { cout << snmv << endl; j = getchar(); exit(1); }
	novNapMas(tnav);
	for (k = 0; k<cemv; k++) stchsrver[k] = hf; //for (k = 0; k<cemv; k++) { tm = epsisredver(etev[k], tkuscv, kuscv, dmkooscv, dkoscvt, dkoscvm, dkoscvl); stchsrver[k] = tm; } //for (k=0; k<cemv; k++) cout << "tem = " << etev[k] << "\tst_ch = " << stchsrver[k] << "\t"; cout << endl;
	napMasEKTPVerNac(wmg, wsi, wal); for (k=0; k<cemv; k++) cout << "te = " << etev[k] << "\tktp = " << ektpv[k] << endl; 
}
void novNapMas(double tnach)
{
	NapMasVozdSha(ktpvover, vtever, dmkvover);
	int j, k;
	if ((!vysv) || (vysv==1)) {  //èñõîäíûé èëè ïîñëå ïîâòîðíûõ èçìåðåíèé
	if (!vyfv) { if (!vpmf) porver = por207; else if (vpmf==1) porver=por16035; }
	else if (vyfv==1) { if (!vpkf) porver = poris84; else if (vpkf==1) porver=porin84; }
	else if (vyfv==2) porver = poro16035; }
	if (vysv==2) { //ïîñëå îáæèãà
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
	double dhk = y0ver / fabs(nxtfv - nnxtfv), hnver = nnxtfv*dhk, ko=1e-6, hvko = (13e1)*ko, srp, marp;
	double hvh = ko, hvna, p, r, d = 0.0, *atr = NULL, t, *po, **mu, *ras, *srra, *legr, *prgr; 
	for (j = 0; j < ks; j++) d = d + hf; //cout << "d = " << d << "\t";
	dkospv = RaschAlphaTvKarVer(); 
	vtvk=cemv;
	for (vtv = vtvn; vtv<vtvk; vtv++) { //ïðîáåãàåì ïî òåìïåðàòóðå
		mu=rasPorpoRazVer(porver, 0, 1, vysv, isrp, vpkf);
		k=0; ras=mu[k]; k++; srra=mu[k]; k++; prgr=mu[k]; k++; legr=mu[k]; k++; 
		po=mu[k]; j=0; srp=po[j]; if (po) delete[]po; k++; po=mu[k]; marp=po[j]; if (po) delete[]po; if (mu) delete[]mu;
		cout << "Sred razm por = " << ko << "\t";
		r = d; t = ko; po = KorrZnachVozdProsVer(ko, r, porver, 0); ko=po[0]; cout << "Korr sred razm por = " << ko << "\t";
		r = d; ko = t; po = KorrZnachVozdProsVer(ko, r, porver, 1); p=po[vtv]; cout << "Dol Plo = " << p << "\t";
		r = d; ko = t; ko = KorrZnachVozdProsVermik(ko, r, porver); hvna = ko; hvko = ko;
		if (vtv>vtvn) {
			q = jk - jn; ka = (tholv[vtv] - tgorv[vtv]) / y0ver; kb = tgorv[vtv]; hkver = hnver;
			for (k = 0; k < q; k++) { temrasv[k] = kb + ka*hkver; hkver = hkver + dhk; }
		}
		hvove = hvna; hver = hvove*(hf - porver) / porver;
		while (hvove <= hvko) {
			j = jn; hkver = hnver*dhk; //ïðîáåãàåì ïî ðàçìåðàì ïîð
			while ((j < jk) && (hkver < y0ver)) {
				cout << "hk = " << hkver << endl; sctxv = j; //ïðîáåãàåì ïî êîîðäèíàòå
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
		pt = FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, cc, vy, qobv, etev, cemv, 0, snmv, dpctv, vtv, ektpv, y0ver, ktrv, kttkv, ks, ktpvover, vtever, dmkvover, ecktpver, hver, hvove, Rver, Tver, Aver, Rtver, Ttver, Atver, txver, qxver, dtxver, sctxv, tocrasver, porver, 2 * N, stchsrver, Tpctv, b0, c0, mk, svfdv, 0);
		k = 3; fc = pt[k]; delete[]pt;
		if (fabs(fc) < 1e0) { c0 = cc; break; }
		pt = FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, ca, vy, qobv, etev, cemv, 0, snmv, dpctv, vtv, ektpv, y0ver, ktrv, kttkv, ks, ktpvover, vtever, dmkvover, ecktpver, hver, hvove, Rver, Tver, Aver, Rtver, Ttver, Atver, txver, qxver, dtxver, sctxv, tocrasver, porver, 2 * N, stchsrver, Tpctv, b0, c0, mk, svfdv, 0);
		k = 3; fa = pt[k]; delete[]pt;
		if (fabs(fa) < 1e0) { c0 = ca; break; }
		pt = FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, cb, vy, qobv, etev, cemv, 0, snmv, dpctv, vtv, ektpv, y0ver, ktrv, kttkv, ks, ktpvover, vtever, dmkvover, ecktpver, hver, hvove, Rver, Tver, Aver, Rtver, Ttver, Atver, txver, qxver, dtxver, sctxv, tocrasver, porver, 2 * N, stchsrver, Tpctv, b0, c0, mk, svfdv, 0);
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
		pt = FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, c0, vy, qobv, etev, cemv, 0, snmv, dpctv, vtv, ektpv, y0ver, ktrv, kttkv, ks, ktpvover, vtever, dmkvover, ecktpver, hver, hvove, Rver, Tver, Aver, Rtver, Ttver, Atver, txver, qxver, dtxver, sctxv, tocrasver, porver, 2 * N, stchsrver, Tpctv, b0, c0, mk, svfdv, 0);
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
		if (!vyve) mkove[k] = ReflSredSha(eteve[k]);
		else if (vyve == 1) mkove[k] = ReflSredVer(eteve[k]);
		else if (vyve == 2) mkove[k] = ReflSreditom(eteve[k], vyma);
		else if (vyve == 3) mkove[k] = ReflSredkvi(eteve[k], vyma);
	}
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
double **izmRTAVer(double *tere, int kost, int izm, double *Ra, double *Ta, double *Aa, double *Rb, double *Tb, double *Ab, int v) //izm = 0 - íåò èçìåíåíèé, izm - ó÷èòûâàþòñÿ èçìåíåíèÿ //ïîèñê èçìåíåíèÿ ñòåïåíè ÷åðíîòû èëè áåçðàçìåðíîãî êîýôôèöèåíòà ïîãëîùåíèÿ
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
	kp = roo; while (kp) { ne = kp->nex; delete kp; kp = ne; } //óäàëåíèå ñïèñêà
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
		prot = zadrktVer(1, kost, d, 0, htk, hvp, 0, 0, tc, 1, 0, prot); //vyv - âûáîîð âåùåñòâà: 0 - âåðìèêóëèò, 1 - øàìîò
		mu = RaschRTAVer(kost, htk, 0.0, 0.0, 1, 0, hvp, 0, tc, 0, 0, 1);
	}
	else if (vyv == 1) {
		prot = zadrktSha(1, kost, d, 0, htk, hvp, 0, tc, 0, 0, prot);
		mu = RaschRTASha(kost, htk, tc, 0.0, hvp, tc, 0, 0, 0, 1);
	}
	else if (vyv == 2) {
		prot = zadrktitom(1, kost, d, vyte, htk, hvp, 0, 0, tc, 0, 0, prot);
		mu = RaschRTAitom(kost, htk, 0.0, 0.0, 1, 0, hvp, 0, tc, 0, 0, 1);
	}
	else if (vyv == 3) { 
		prot = zadrktkvi(1, kost, d, vyte, htk, hvp, 0, 0, tc, 0, 0, prot); 
		mu = RaschRTAkvi(kost, htk, 0.0, 0.0, 1, 0, hvp, 0, tc, 0, 0, 1);
	}
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
	double x = 0.0, y = 0.0, p = 0.0, xt = 0.0, yp = 0.0, hf = 1e0, tc = 22.0 + tev0, altc=1e-6;
	double *srra = NULL, *porv = new double[nfr], srp;
	double *pn = new double[d], *alx = new double[d], *rcf = new double[cfk];
	double *rapo = NULL, **unau, ce, cet, *pp = NULL, lamtem; 
	double *legr = NULL, *prgr=NULL, dvpn = 1448.0*altc / tc / 4e0; //ìåíüøå 1 %
	unau=rasPorpoRazVer(porver, vyfv, 1, vysv, isrp, vpkf);
	k=0; rapo=unau[k]; k++; srra=unau[k]; k++; prgr=unau[k]; k++; 
	legr=unau[k]; k++; pp=unau[k]; srp=pp[0]; if (pp) delete[]pp; k++; pp=unau[k];
	ce = pp[0]; cet = ce; if (pp) delete[]pp; if (unau) delete[]unau;
	j = 0; while (cet > 0.0) { cet = cet - hf; j++; } cei = j; //cout << "cem_srp = " << cei << "\tsrp = " << srp << "\t"; cout << endl; //for (j=0; j<cei; j++) if (j<10) cout << "j = " << j << "\trpr = " << rapo[j] << "\tlegr = " << legr[j] << "\t"; cout << endl;
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
			if (fr == 3) x = RasFracXeffVer150_200(k); //ðàçìåð ÷àñòèöû
			rcf[k + fr*l] = x*(1e-6); y = x*porv[fr]; //ðàçìåð ïîðû //cout << "x = " << x << "\ty = " << y << "\tp = " << porv[fr] << endl;
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
	x = x / altc; printf("sr_okp = %0.10lf\n", x); //îñëàáëåíèå ÊÏ çà ñ÷åò ïîðèñòîé ñòðóêòóðû âåðìèêóëèòà
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
	for (j = 0; j < cemv; j++) { //îïðåäåëÿåì ÝÊÒÏ ìíîãîñëîéíîé ñòåíêè
		tepo = opredTempStenShaFragm(tepo, 2 * kost, ktpvover, vtever, etev, kttkv, cemv, dmkvover, htch, hvozd, qobv[j], etev[j], -hf); //â ñåðåäèíå ñëîÿ
		r = 0.0; for (k = 0; k < ksu; k++) for (f = 0; f<ksu; f++) { p = tepo[k] - tepo[f]; if (p>r) r = p; }
		p = hvozd*(d - hf) + htch*d; r = r / p; t = qobv[j] / r; ecktpver[j] = t;
	}
	delete[]tepo;
	f = 1000; for (j = 0; j<cemv; j++)
	{
		sa = 0.0; sb = 1e0; k = 0;
		do {
			sc = (sa + sb) / 2e0;
			fa = kttkv[j] * sa + ecktpver[j] * (hf - sa) - ektpv[j]; //ýôôåêòèâíûå ÊÒÏ ìíîãîñëîéíîé ñòåíêè è  ïåðåìû÷êè äîëæíû ñðàâíÿòüñÿ
			fb = kttkv[j] * sb + ecktpver[j] * (hf - sb) - ektpv[j]; //÷òîáû íàéòè îòíîñèòåëüíûå äîëè ïëîùàäåé ñå÷åíèÿ ïåðåíîñà îáùåé ÏÒÏ
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
	while ((ra>e) && (j < k)) { //ïîäòÿãèâàåì ïîðèñòîñòü ê çíà÷åíèþ, êîòîðîå çàäàëè èçíà÷àëüíî, âî âðåìÿ ïîäñòðîéêè ÝÊÒÏ
		pc = (pa + pb) / 2e0;
		kc = hps*pc / por;
		po = oprEffDoliTepPerenVer(kc, ksf, pc); tc = po[q]; delete[]po; //ïðè 373 Ê
		tcc = kc*(hf - pc) / pc;
		fc = (hf - tc)*kc / (kc + tcc) - por;
		ka = hps*pa / por;
		po = oprEffDoliTepPerenVer(ka, ksf, pa); ta = po[q]; delete[]po; //îïðåäåëÿåì äîëþ ïëîùàäè ñå÷åíèÿ ïåðåìû÷êè
		tca = ka*(hf - pa) / pa;
		fa = (hf - ta)*ka / (ka + tca) - por;
		kb = hps*pb / por;
		po = oprEffDoliTepPerenVer(kb, ksf, pb); tb = po[q]; delete[]po; //÷åðåç ïåðåìû÷êó òåïëî ðàñïðîñòðàíÿåòñÿ ÷èñòîé òåïëîïðîâîäíîñòüþ
		tcb = kb*(hf - pb) / pb;
		fb = (hf - tb)*kb / (kb + tcb) - por;
		if ((fc*fb > 0.0) && (fa*fc<0.0)) pb = pc; if ((fc*fa>0.0) && (fb*fc < 0.0)) pa = pc;
		j++; ra = fabs(pa - pb);
	}
	dpctvm[q] = tc; cout << "Dol Plo CTP = " << tc << endl; }
	for (j=0; j<cemv; j++) cout << "EC KTP TK ( " << j << " ) = " << ecktpver[j] << endl;
	if (!vy) { pkc[0]=kc; return pkc; } //ñêîððåêòèðîâàííîå çíà÷åíèå òîëùèíû âîçäóøíîé ïðîñëîéêè (ðàçìåð ïîðû), êîãäà ââåëè ïåðåìû÷êó
	else if (vy == 1) return dpctvm;
} //äîëÿ ïëîùàäè, ÷åðåç êîòîðóþ ïðîèñõîäèò ïåðåíîñ òåïëà ÷èñòîé òåïëîïðîâîäíîñòüþ
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
	for (k = 1; k < kost; k++) xi[k] = xi[k - 1] + ht; //ìàññèâ ñåðåäèí êàæäîé èç ñòåíîê ïî òîëùèíå
	for (k = 0; k<kost; k++) teks[k] = koe[0] * xi[k] + koe[1]; //for (k=0; k<kost; k++) cout << "k = " << k << "\txk = " << xi[k] << "\ttex = " << teks[k] << endl; //ëèíåàðèçàöèÿ ïîëÿ òåìïåðàòóð
	qobver = opredKTPTKTochSha(qobv, etev, (teks[0] + teks[kost - 1]) / 2e0, ce); return teks;
}
double KorrZnachVozdProsVermik(double hps, double ksf, double por)
{
	int j = 0, k = 1000, q=0; double pa = 1e-3, pb = 1e0, pc, ra = fabs(pa - pb); dpctv=opredKTPTKTochSha(dpctvm, etev, etev[vtv], cemv);
	double fa, fb, fc, tca, tcb, tcc, ka = hps*pa / por, kb = pb*hps / por, kc, e = tocrasver; //cout << "hps = " << hps << "\tksuf = " << ksf << "\tpor = " << por << endl;
	while ((ra>e) && (j < k)) { //ïîäòÿãèâàåì ïîðèñòîñòü ê çíà÷åíèþ, êîòîðîå çàäàëè èçíà÷àëüíî, âî âðåìÿ ïîäñòðîéêè ÝÊÒÏ
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
	} //ñêîððåêòèðîâàííîå çíà÷åíèå òîëùèíû âîçäóøíîé ïðîñëîéêè (ðàçìåð ïîðû), êîãäà ââåëè ïåðåìû÷êó
	return kc;
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
	k=0; nvyfv[k]=0; k++; nvyfv[k]=1; //ôðàêöèè âåðìèêóëèòà
	k=0; nvysv[k]=0; k++; nvysv[k]=1; k++; nvysv[k]=2; //ñîñòîÿíèÿ âåðìèêóëèòà
	k=0; nvmivmf[k]=1; k++; nvmivmf[k]=2; //ñòàöèîíàðíûå ìåòîäû èçìåðåíèé - 2019 è 2020
	k=0; nvyuv[k]=1; k++; nvyuv[k]=2; k=0; j=-1; b=0; //óêëàäêà âåðìèêóëèòà
	for (kvf=0; kvf<nnvyfv; kvf++) {
		vyfrve=nvyfv[kvf];
		for (jvsv=0; jvsv<nnvysv; jvsv++) {
			vysove=nvysv[jvsv];
			if ((!vyfrve) || (vyfrve==2)) {
				for (qvmi=0; qvmi<nnvmivmf; qvmi++) {
					vymivmf=nvmivmf[qvmi]; 
					muv=new double*[f]; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
					muv=napMasEKTP(vyfrve, vysove, vymivmf, efte, k, h, vyukve, n, k, muv, f, dmko); muu[b]=muv; b++; if (vysove>0) break; } }
			else if (vyfrve==1) {
				for (qvuv=0; qvuv<nnvyuv; qvuv++) {
					vyukve=nvyuv[qvuv];
					muv=new double*[f]; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
					muv=napMasEKTP(vyfrve, vysove, vymivmf, efte, k, h, vyukve, n, k, muv, f, dmko); muu[b]=muv; b++; } } } } 
	for (j=0; j<b; j++) { muv=muu[j];
	d=0; temvct=muv[d]; d++; temvht=muv[d]; d++; tepvt=muv[d]; d++; ktp=muv[d]; d++; ts=muv[d]; d++; po=muv[d]; k=0; nf=po[k]; 
					cf=e; while (cf<nf) { cf=cf+hf; k++; } q=k; 
					d=1; k=0; cfp=tepvt[k]; for (k=d; k<q; k++) { cft=tepvt[k]; if ((cft<=cfp) && (d>0)) { d=-1; break; } cfp=cft; } //for (k=0; k<q; k++) cout << "c = " << c << "\tqo = " << tepvt[k] << "\tts = " << ts[k] << "\t"; //vyvodfile(ts, q, k, nf, "C:\\Users\\Àíäðåé\\Documents\\_Àñïèðàíòóðà\\tmp\\tmp.txt");
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
double **opredtemphc(double *effktp, double *efftem, double *tgv, double *thv, double *qon, int dlma, int n, double h, double *qob, double *tgor, double *thol, double **mu) //n=3 - äëèíà ìàññèâà êîýôôèöèåíòîâ ïðèáëèæàþùåãî ìíîãî÷ëåíà, tgv - òåìïåðàòóðà ãîðÿ÷åé ñòåíêè, thv - òåìïåðàòóðà õîëîäíîé ñòåíêè, dlma - äëèíà ìàññèâà ÝÊÒÏ
{
	int k=0, j=0; double *kq=new double[n], *kho=new double[n], *kgo=new double[n];
	double *tesr=new double[n], ts=0.0, g=0.0, p=0.0, hf=1e0, r=0.0, s=0.0, tego=0.0, teho=0.0;
	if ((!kq) || (!kho) || (!kgo) || (!tesr)) { cout << snmv << endl; j=getchar(); exit(1); }
	for (k=0; k<n; k++) tesr[k]=(tgv[k]+thv[k])/2e0;
	kq=koefPribSha(qon, tesr, n, kq, snmv);
	kgo=koefPribSha(tgv, tesr, n, kgo, snmv);
	kho=koefPribSha(thv, tesr, n, kho, snmv);
	for (k=0; k<dlma; k++) 
	{
		ts=efftem[k];
		g=0.0; p=0.0; for (j=0; j<n; j++) { g=g+pow(ts, p)*kq[j]; p=p+hf; } 
		if (g<0.0) g=0.0; if (qob) qob[k]=g;
		p=effktp[k];
		g=fabs(g*h/p); 
		tego=ts+g/2e0; teho=ts-g/2e0;
		if ((tego<0.0) || (teho<0.0)) 
		{
		if (tego<0.0) 
		{
			r=0.0; s=0.0; for (j=0; j<n; j++) { r=r+pow(ts, s)*kgo[j]; s=s+hf; } 
			if (r<0.0) r=0.0; tego=r; 
		}
		if (teho<0.0) 
		{
			r=0.0; s=0.0; for (j=0; j<n; j++) { r=r+pow(ts, s)*kho[j]; s=s+hf; } 
			if (r<0.0) r=0.0; teho=r; 
		} 
		}
	if (tgor) tgor[k]=tego; if (thol) thol[k]=teho; }
	if (tesr) delete[]tesr; if (kq) delete[]kq; 
	if (kho) delete[]kho; if (kgo) delete[]kgo;
	k=0; mu[k]=thol; k++; mu[k]=tgor; k++; mu[k]=qob; k++; mu[k]=effktp; k++; mu[k]=efftem;
	return mu;
}
double **opredTempHolGor(double *ektp, double *ete, int n, int l, double h0, int v, double **mu, int rmu, int ni, double *efte, int dmkoef) //ìîäåëèðîâàíèå ïðîöåññà òåïëîîáìåíà â îáðàçöå
{ //n - äëèíà ìàññèâà ektp, l - äëèíà ìàññèâîâ temvs, qob, ktp, temvh, temvc, ni - äëèíà efte, qob - ïëîòíîñòü òåïëîâîãî ïîòîêà, êîòîðóþ ìîæåò ñîçäàòü ëàáîðàòîðíàÿ óñòàíîâêà 
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
			qob[k]=g; } //cout << "n = " << n << endl; for (k=0; k<n; k++) cout << "qob = " << qob[k] << "\tektp = " << ektp[k] << "\tte = " << ete[k] << endl; 
		if (koeq) delete[]koeq; for (k=0; k<cemf; k++) { po=muv[k]; delete[]po; } delete[]muv;
		g=0.0; for (k=0; k<n; k++) {
			if (qob[k]>ep) {
				laef=ektp[k]; p=0.0; Thnac=Thna+g*dt; del=hf; etem=ete[k]; ktp[k]=laef;
				while ((del>ep) && (p<nit)) {
					thol=Thnac+p*d; //Tg - ìàññèâ òåìïåðàòóð ãîðÿ÷èõ ñòåíîê, Th - ìàññèâ òåìïåðàòóð õîëîäíûõ ñòåíîê
					tgor=thol+qob[k]*h0/laef;
					del=(2e0*etem-(thol+tgor));
					p=p+hf; } 
				g=g+hf; }
			else { thol=0.0; tgor=0.0; qob[k]=0.0; ktp[k]=0.0; } tsred=(tgor+thol)/2e0; 
			temvc[k]=thol; temvh[k]=tgor; temvs[k]=tsred; }  
	g=0.0; for (k=0; k<n; k++) g=g+hf; //cout << "n = " << g << endl;
	k=1; po=new double[k]; muv=new double*[cemf]; if ((!po) || (!muv)) { cout << snm << endl; k=getchar(); exit(1); } k=0; po[k]=g;
	muv[k]=temvc; k++; muv[k]=temvh; k++; muv[k]=qob; k++; muv[k]=ktp; k++; muv[k]=temvs; k++; muv[k]=po; for (k=0; k<cemf; k++) mu[k]=muv[k]; //for (k=0; k<n; k++) cout << "tc = " << temvc[k] << "\ttg = " << temvh[k] << "\tts = " << temvs[k] << endl; 
	k=0; j=1; mu=vydelPol(k, n, muv, mu, rmu, j); 
	if (temvc) delete[]temvc; if (temvh) delete[]temvh; //if (qob) delete[]qob; 
	if (ktp) delete[]ktp; if (temvs) delete[]temvs; if (po) delete[]po; if (muv) delete[]muv; 
	return mu;
}
void napMasEKTPVerNac(double wmg, double wsi, double wal)
{
	int k=1, nk=0, j=0, jk=0, f=cemdu, qn=0, q=0; 
	double hf=1e0, *po=NULL, s=0.0, *etevv, *isv, *mkovv;
	double **mu=new double*[f], **muv=new double*[f]; if ((!mu) || (!muv)) { cout << snmv << endl; k=getchar(); exit(1); } 
	double nf=0.0, t=0.0, g=0.0, e=1e-1;
	por=novNapMas(tnac, k, j, vmivmf, vyfv, vyuv, vysv, cem);
	muv=napMasEKTP(vyfv, vysv, vmivmf, ete, k, ys0, vyuv, cem, k, muv, f, dmko);
	k=0; j=1; mu=vydelPol(k, cem, muv, mu, f, j);
	k=0; thol=mu[k]; k++; tgor=mu[k]; k++; qob=mu[k]; k++; 
	ektpv=mu[k]; k++; tsred=mu[k]; ete=tsred; k++; po=mu[k]; 
	k=0; nf=po[k]; t=0.0; k=0; while (t<nf) { t=t+hf; k++; } cem=k; //cout << "cem = " << cem << "\tnf = " << nf << "\t";	
	for (k=0; k<cem; k++) cout << "tem = " << ete[k] << "\tktp = " << ektpv[k] << "\t"; cout << endl;
	for (k=0; k<f; k++) { po=muv[k]; delete[]po; } 
	if (!((!vysv) && (!vmivmf) && (!vyfv))) for (k=0; k<f; k++) { po=mu[k]; delete[]po; }
	if (muv) delete[]muv; if (mu) delete[]mu;
	double nf = 0.0, tn=0.0, g=0.0;
	etevv=new double[cemv]; isv=new double[cemv]; mkovv=new double[cemv]; if ((!etevv) || (!isv) || (!mkovv)) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k=0; k<cemv; k++) { etevv[k]=etev[k]; isv[k]=stchsrver[k]; mkovv[k]=mkov[k]; }
	jk=cemv; k=0; nf=po[k]; if (po) delete[]p; k = 0; while (nf>0.0) { nf = nf - hf; k++; } nk = k; cemv = nk; tn = tsredver[0]; tnav = tn; 
	if (etev) delete[]etev; if (mkov) delete[]mkov; if (stchsrver) delete[]stchsrver; 
	etev = new double[cemv]; kttkv = new double[cemv]; mkov = new double[cemv]; stchsrver = new double[cemv]; if ((!stchsrver) || (!etev) || (!kttkv) || (!mkov)) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < cemv; k++) { stchsrver[k] = 0.0; kttkv[k] = 0.0; mkov[k] = 0.0; etev[k] = 0.0; }
	etev[0] = tn; for (k = 1; k<cemv; k++) etev[k] = etev[k - 1] + detev; 
	for (j=0; j<cemv; j++) { f=1; nk=0; for (k=0; k<jk; k++) 
		if (fabs(etev[j]-etevv[k])<hf) { f=-1; nk=k; break; } 
		if (f<0) { stchsrver[j]=isv[nk]; mkov[j]=mkovv[nk]; } else { stchsrver[j]=0.0; mkov[j]=0.0; } } for (k = 0; k<cemv; k++) { g = epsisredver(etev[k], tkuscv, kuscv, dmkooscv, dkoscvt, dkoscvm, dkoscvl); stchsrver[k] = g; } 
	if (etevv) delete[]etevv; if (mkovv) delete[]mkovv; if (isv) delete[]isv; for (k=0; k<cemv; k++) cout << "tem = " << etev[k] << "\tst_ch = " << stchsrver[k] << "\t"; cout << endl;
	kttkv = opredKTPTverKarkVerm(tsredver, ektpv, porver, wsi, wal, wmg, vyfv, dmkov, tnoscv, dtoscv, dmkooscv, kuscv, tkuscv, stchsrver, cemv, vysv, vpkf, isrp); //for (k=0; k<cemv; k++) cout << "tem = " << etev[k] << "\tkttkv = " << kttkv[k] << endl;
}
double **napMasEKTP(int vyfrve, int vysove, int vymeizvemafr, double *te, int v, double h, int vyukve, int nt, int wa, double **mu, int rmu, int dmk)
{
	int k=0, nvm=dmk, j=0, vy=0, n84=0, n207=0, cedu=rmu, q=0, cedumi=cemdum, nn=0;
	double F=13.85*1e-4, nf84=0.0, nf207=0.0, *po=NULL, *tc84=NULL, s=0.0, t=0.0, e=1e-1; 
	double *th84=NULL, hf=1e0, *tp84=NULL, *th207=NULL, *tc207=NULL, *tp207=NULL;
	double *koeft=NULL, *koefh=NULL, *koefc=NULL, *koefq=NULL, *koefs=NULL, **muv=NULL, r=0.0; //cout << "vf = " << vyfrve << "\tvs = " << vysove << "\tvmi = " << vymeizvemafr << "\tvu = " << vyukve << endl; //for (k=0; k<nt; k++) cout << "te ( " << k << " ) = " << te[k] << "\t";
	double *ktp=NULL, *temvs=NULL, *temvh=NULL, *temvc=NULL, *tepv=NULL, *vm=NULL, *ma=NULL; 
	//-------
	if ((!vyfrve) || (vyfrve==2)) { //ôðàêöèÿ 2-0,7 ìì èëè ôðàêöèÿ 1,6-0,35 ìì
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
		if (!vymeizvemafr) { //óñòàíîâêà Netzsch
		if (muv) delete[]muv; muv=new double*[cedu]; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
		muv=oprkoefKTPiskhchao(vymeizvemafr, k, te, h, n207, muv, cedu, nt, dmk); 
		k=0; koefc=muv[k]; k++; koefh=muv[k]; k++; koefq=muv[k]; k++;
		koeft=muv[k]; k++; koefs=muv[k]; k++; po=muv[k]; }
		else if (vymeizvemafr==1) { //äàííûå 2020 ãîäà - ÃÎÑÒ 12170 - ñòàöèîíàðíûé ìåòîä
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
		else if (vymeizvemafr==2) { //äàííûå 2019 ãîäà - ÃÎÑÒ 12170
		koefq=danPoTemTepl2071(temvs, tepv, nvm); koefh=danPoTemTepl2071(temvs, temvh, nvm);
		koefc=danPoTemTepl2071(temvs, temvc, nvm); koeft=danPoTemTepl2071(temvs, ktp, nvm);
		koefs=danPoTemTepl2071(temvs, temvs, nvm); } }
				if (vysove==1) { //ïîñëå ïîâòîðíûõ èçìåðåíèé
			koefc=danPoTemTepl2072(temvs, temvc, nvm); koefh=danPoTemTepl2072(temvs, temvh, nvm); 
			koefq=danPoTemTepl2072(temvs, tepv, nvm); koeft=danPoTemTepl2072(temvs, ktp, nvm); 
			koefs=danPoTemTepl2072(temvs, temvs, nvm); }
				else if (vysove == 2) { //ïîñëå îáæèãà ïðè 1000 °Ñ
			koefc=danPoTemTepl2073(temvs, temvc, nvm); koefh=danPoTemTepl2073(temvs, temvh, nvm); 
			koefq=danPoTemTepl2073(temvs, tepv, nvm); koeft=danPoTemTepl2073(temvs, ktp, nvm); 
			koefs=danPoTemTepl2073(temvs, temvs, nvm); }
				else if (vysove == 3) 
				{ //ïîñëå ïîâòîðíîãî îáæèãà ïðè 1000 °Ñ
			koefc=danPoTemTepl2074(temvs, temvc, nvm); koefh=danPoTemTepl2074(temvs, temvh, nvm); 
			koefq=danPoTemTepl2074(temvs, tepv, nvm); koeft=danPoTemTepl2074(temvs, ktp, nvm); 
			koefs=danPoTemTepl2074(temvs, temvs, nvm); } }
	else if (vyfrve==1) { //ôðàêöèÿ 8-4 ìì
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
		if (vyukve==1) { //ïëîñêî-ïàðàëëåëüíàÿ çàñûïêà
			if (!vysove) { //èñõîäíûé
				koefc=danPoTemTepl840(temvs, temvc, nvm); koefh=danPoTemTepl840(temvs, temvh, nvm); 
				koefq=danPoTemTepl840(temvs, tepv, nvm); koeft=danPoTemTepl840(temvs, ktp, nvm); 
				koefs=danPoTemTepl840(temvs, temvs, nvm); }
			else if (vysove==1) { //ïîñëå ïîâòîðíûõ èçìåðåíèé
				koefc=danPoTemTepl842(temvs, temvc, nvm); koefh=danPoTemTepl842(temvs, temvh, nvm); 
				koefq=danPoTemTepl842(temvs, tepv, nvm); koeft=danPoTemTepl842(temvs, ktp, nvm); 
				koefs=danPoTemTepl842(temvs, temvs, nvm); }
			else if (vysove==2) { //ïîñëå îáæèãà
				koefc=danPoTemTepl845(temvs, temvc, nvm); koefh=danPoTemTepl845(temvs, temvh, nvm); 
				koefq=danPoTemTepl845(temvs, tepv, nvm); koeft=danPoTemTepl845(temvs, ktp, nvm); 
				koefs=danPoTemTepl845(temvs, temvs, nvm); } }
		else if (vyukve==2) { //âåðòèêàëüíàÿ çàñûïêà
			if (!vysove) { //èñõîäíûé
				koefc=danPoTemTepl841(temvs, temvc, nvm); koefh=danPoTemTepl841(temvs, temvh, nvm); 
				koefq=danPoTemTepl841(temvs, tepv, nvm); koeft=danPoTemTepl841(temvs, ktp, nvm); 
				koefs=danPoTemTepl841(temvs, temvs, nvm); }
			else if (vysove == 1) { //ïîñëå ïîâòîðíûõ èçìåðåíèé
				koefc = danPoTemTepl844(temvs, temvc, nvm); koefh = danPoTemTepl844(temvs, temvh, nvm); 
				koefq = danPoTemTepl844(temvs, tepv, nvm); koeft = danPoTemTepl844(temvs, ktp, nvm); 
				koefs = danPoTemTepl844(temvs, temvs, nvm); }
			else if (vysove == 2) { //ïîñëå îáæèãà
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
double **oprkoefKTPiskhchao(int vmiv, int v, double *efte, double h, int n207, double **mu, int rmu, int cem, int dmk) //vmiv - âûáîð ìåòîäà èçìåðåíèé
{
	int nk=0, k=0, j=0, q=0, qn=0, nn=n207, cedumi=cemdum, w=rmu-1, f=rmu; 
	double hf=1e0, nf207=0.0, *po=NULL, *koeft=NULL, *koefh=NULL, s=0.0, t=0.0;
	double *koefc=NULL, *koefs=NULL, *koefq=NULL, **muv=NULL, **muvv=NULL, r=0.0, e=1e-1;
	double *tsv=NULL, *tgv=NULL, *thv=NULL, *qov=NULL, *ktpv=NULL, *mt=NULL, *ts=NULL; 
	if (!vmiv) {  //0 - óñòàíîâêà Netzsch - íåñòàöèîíàðíûé ìåòîä
		k=cedumi; muv=new double*[k]; if (!muv) { cout << snmv << endl; k = getchar(); exit(1); }
		muv=arrTem_Netzsch(muv);
		k=0; po=muv[k]; k++; mt=muv[k];
		k=0; s=e; nf207=po[k]; while (s<nf207) { s=s+hf; k++; } nk=k; //cout << "nk = " << nk << "\t"; //nk=8 - äëèíà ìàññèâà ktpv
		ktpv=arrKTP_Netzsch(); 
		if (muv) delete[]muv; if (po) delete[]po;
		k=rmu; muv=new double*[k]; k=0; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
		muv=opredTempHolGor(ktpv, mt, nk, nn, h, k, muv, rmu, cem, efte, dmk); //cem - äëèíà ìàññèâà efte //for (k=0; k<f; k++) { po=muv[k]; delete[]po; }
		if (ktpv) delete[]ktpv; if (mt) delete[]mt;
		k=0; thv=muv[k]; k++; tgv=muv[k]; k++; qov=muv[k]; k++; ktpv=muv[k]; k++; 
		tsv=muv[k]; k++; po=muv[k]; if (muv) delete[]muv;
		q=0; nf207=po[q]; s=0.0; while (s<nf207) { s=s+hf; q++; } qn=q; 
		koefh=new double[dmk]; koefc=new double[dmk]; koefq=new double[dmk]; koeft=new double[dmk]; koefs=new double[dmk];
		if ((!koefh) || (!koefc) || (!koefq) || (!koeft) || (!koefs)) { cout << snmv << endl; k=getchar(); exit(1); }
		for (k=0; k<dmk; k++) { koefh[k]=0.0; koefc[k]=0.0; koefq[k]=0.0; koeft[k]=0.0; koefs[k]=0.0; }
		koefc=koefPribSha(thv, tsv, qn, koefc, snmv); koefh=koefPribSha(tgv, tsv, qn, koefh, snm); 
		koefq=koefPribSha(qov, tsv, qn, koefq, snmv); koeft=koefPribSha(ktpv, tsv, qn, koeft, snm); 
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
		p=4; q=5; m[k]=(ma[p]*x[p]+ma[q]*x[q])/(x[p]+x[q]); //600, 800 è 1000 °C
		if (ma) delete[]ma; return m;
	}
	else return ma;
}
//---------
double **arrTem_Netzsch(double **mu)
{
	int k = 0, le = 8; double *tem = new double[le], hf = 1e0, nf = 0.0, *po = new double[1];
	if ((!tem) || (!po)) { cout << snm << endl; k = getchar(); exit(1); }
	tem[k]=27.967;  k++; tem[k]=93.833;  k++; tem[k]=192.5; k++; 
	tem[k]=341.667; k++; tem[k]=491.467; k++; tem[k]=641.2; k++; 
	tem[k]=790.933; k++; tem[k]=991.133; k++;
	nf=0.0; for (k=0; k<le; k++) { tem[k]=tem[k]+te0; nf=nf+hf; } k=0; po[k]=nf;
	k=0; mu[k]=po; k++; mu[k]=tem;
	return mu;
}
double *arrKTP_Netzsch()
{
	int k=0, le=8; double *ktp=new double[le]; 
	if (!ktp) { cout << snm << endl; k = getchar(); exit(1); }
	ktp[k] = 7.0*1e-2;   k++; ktp[k] = 8.0*1e-2;   k++; ktp[k] = 103.0*1e-3; k++; 
	ktp[k] = 15.0*1e-2;  k++; ktp[k] = 207.0*1e-3; k++; ktp[k] = 283.0*1e-3; k++; 
	ktp[k] = 373.0*1e-3; k++; ktp[k] = 477.0*1e-3; 
	return ktp;
}
double *arrKTP_2020()
{
	int k=0, n=6; double *a = new double[n]; if (!a) { cout << snm << endl; k = getchar(); exit(1); }
	a[k] = 0.175566644058715; k++; a[k] = 0.176801537812368; k++; a[k] = 0.179324717653617; k++;
	a[k] = 0.211768953068592; k++; a[k] = 0.237194543621728; k++; a[k] = 0.237231989760775;
	return a;
}
double *arrTem1_2020()
{
	int k=0, n=6; double *a=new double[n]; if (!a) { cout << snm << endl; k = getchar(); exit(1); }
	a[k] = 585.0; k++; a[k] = 600.0;  k++; a[k] = 585.0; k++; 
	a[k] = 800.0; k++; a[k] = 1000.0; k++; a[k] = 1000.0;
	for (k = 0; k < n; k++) a[k] = a[k] + te0;
	return a;
}
double **arrTem2_2020(double **mu)
{
	int k=1, n=6; double *a=new double[n], nf=0.0, hf=1e0, *po=new double[k];
	if ((!a) || (!po)) { cout << snm << endl; k = getchar(); exit(1); }
	for (k = 0; k < n; k++) nf = nf + hf; k=0; po[k]=nf;
	a[k] = 119.75; k++; a[k] = 138.0; k++; a[k] = 129.5; k++; 
	a[k] = 200.0;  k++; a[k] = 273.0; k++; a[k] = 261.0;
	for (k=0; k<n; k++) a[k]=a[k]+te0;
	k=0; mu[k]=po; k++; mu[k]=a;
	return mu;
}
double *arrTem3_2020()
{
	int k=0, n=6; double *a=new double[n]; if (!a) { cout << snm << endl; k = getchar(); exit(1); }
	a[k] = 377.0; k++; a[k] = 396.0; k++; a[k] = 383.5; k++; 
	a[k] = 548.0; k++; a[k] = 703.0; k++; a[k] = 697.25;
	for (k = 0; k < n; k++) a[k] = a[k] + te0;
	return a;
}
//---------
double *danPoTemTepl840(double *temvs, double *temvh, int n) //Çàñûïêà ïëîñêî-ïàðàëëåëüíàÿ, èñõîäíûé
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
double *danPoTemTepl841(double *temvs, double *temvh, int n) //Çàñûïêà âåðòèêàëüíàÿ, èñõîäíûé
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
double *danPoTemTepl842(double *temvs, double *temvh, int n) //Çàñûïêà ïëîñêî-ïàðàëëåëüíàÿ, ïîâòîðû
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
double *danPoTemTepl843(double *temvs, double *temvh, int n) //Çàñûïêà âåðòèêàëüíàÿ, ïîñëå îáæèãà ïðè 1000 °Ñ
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
double *danPoTemTepl844(double *temvs, double *temvh, int n) //Çàñûïêà âåðòèêàëüíàÿ, ïîâòîðû
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
double *danPoTemTepl845(double *temvs, double *temvh, int n) //Çàñûïêà ïëîñêî-ïàðàëëåëüíàÿ, ïîñëå îáæèãà ïðè 1000 °Ñ
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
double *danPoTemTepl2071(double *temvs, double *tepv, int n) //Çàñûïêà èñõîäíàÿ, ôðàêöèÿ 2-0,7 ìì
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
double *danPoTemTepl2072(double *temvs, double *tepv, int n) //Ôðàêöèÿ 2-0,7 ìì (ïîâòîðíûå èçìåðåíèÿ)
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
double *danPoTemTepl2073(double *temvs, double *tepv, int n) //Ôðàêöèÿ 2-0,7 ìì, ïîñëå îáæèãà ïðè 1000 °Ñ
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
double *danPoTemTepl2074(double *temvs, double *tepv, int n) //Ôðàêöèÿ 2-0,7 ìì, ïîñëå ïîâòîðíîãî îáæèãà ïðè 1000 °Ñ
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
	for (k=0; k<n; k++) temcol[k]=temcol[k] + te0;
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
	for (k=0; k<n; k++) temh[k]=temh[k]+te0;
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
	for (k=0; k<n; k++) temcol[k]=temcol[k]+te0;
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
	for (k=0; k<n; k++) temh[k]=temh[k]+te0;
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
	double *po=NULL, nf=0.0, t=0.0, hf=1e0, *vm=NULL, e=1e-1, r=0.0, s=0.0;
	k=w; po=mu[k]; k=0; nf=po[k]; 
	t=e; while (t<nf) { k++; t=t+hf; } nk=k; mui=new int[nk]; if (!mui) { cout << snmv << endl; k=getchar(); exit(1); }
	q=0; for (k=0; k<nk; k++) {
		x=1; for (m=0; m<w; m++) { 
			po=mu[m]; if ((po[k]<e) && (x>0)) { x=-1; break; } } 
		if (x>0) { mui[k]=k; q++; } else mui[k]=-1; } qn=q;
		if (fl>0) { m=1; po=mu[m]; for (k=0; k<nk; k++) if (po[k]>templa) { mui[k]=-1; qn--; break; } } //if (fl>0) { cout << "qn = " << qn << "\tnf = " << nf << "\tw = " << w << "\tn = " << n << "\t"; } 
	for (m=0; m<w; m++) {
	vm=new double[qn]; if (!vm) { cout << snmv << endl; k=getchar(); exit(1); } //for (k=0; k<qn; k++) vm[k]=0.0;
	q=0; po=mu[m]; for (k=0; k<nk; k++) {
		x=mui[k]; if (x>=0) { vm[q]=po[x]; q++; } } muv[m]=vm; } qn=q; qk=q; if (mui) delete[]mui; 
	nf=0.0; for (k=0; k<qn; k++) nf=nf+hf; 
	k=1; po=new double[k]; if (!po) { cout << snmv << endl; k=getchar(); exit(1); } k=0; po[k]=nf; muv[w]=po; //if (fl>0) { q=w; po=muv[q]; q=0; t=po[q]; cout << "r = " << t << endl; s=e; q=0; while (s<t) { s=s+hf; q++; } qk=q; q=0; po=muv[q]; for (q=0; q<qk; q++) cout << "temc = " << po[q] << "\t"; cout << endl; q=1; po=muv[q]; for (q=0; q<qk; q++) cout << "temh = " << po[q] << "\t"; cout << endl; q=2; po=muv[q]; for (q=0; q<qk; q++) cout << "qo = " << po[q] << "\t"; cout << endl; q=3; po=muv[q]; for (q=0; q<qk; q++) cout << "ktp = " << po[q] << "\t"; cout << endl; q=4; po=muv[q]; for (q=0; q<qk; q++) cout << "tems = " << po[q] << "\t"; cout << endl; q=w; po=mu[q]; q=0; t=po[q]; cout << "r = " << t << endl; s=e; q=0; while (s<t) { s=s+hf; q++; } qk=q; q=0; po=mu[q]; for (q=0; q<qk; q++) cout << "temc = " << po[q] << "\t"; cout << endl; q=1; po=mu[q]; for (q=0; q<qk; q++) cout << "temh = " << po[q] << "\t"; cout << endl; q=2; po=mu[q]; for (q=0; q<qk; q++) cout << "qo = " << po[q] << "\t"; cout << endl; q=3; po=mu[q]; for (q=0; q<qk; q++) cout << "ktp = " << po[q] << "\t"; cout << endl; q=4; po=mu[q]; for (q=0; q<qk; q++) cout << "tems = " << po[q] << "\t"; cout << endl;  cout << mu << "\t" << muv << "\t"; }
	for (k=0; k<f; k++) { vm=mu[k]; //if ((fl>0) && (vm) && (k<w)) { for (q=0; q<nk; q++) cout << "k = "<< k << "\tq = " << q << "\tp = " << vm[q] << "\t"; cout << endl; } 
	if ((fl<0) && (vm)) delete[]vm; } //if (fl>0) { cout << "qn = " << qn << "\tn = " << nf << "\t"; } 
	return muv;
}
