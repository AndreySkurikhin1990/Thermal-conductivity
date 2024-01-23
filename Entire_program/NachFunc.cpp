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
//-----
void zadrktNac(int);
void initarrves(int, int, int, double, int, int, int, int, int, int, int, int, int, double, double, int, double, double, int);
double *koefoslab(double, double, double, double *, int, double *, char *);
double NovNapMas(int, int, int, int, int, int, int, int); 
double RaschAlphaTvKar(int, int, int, int, int, int, int, double, char *, int, double, int, double *, double *, double *, int, 
	int, int, double *, int, double *, double *, double *, double *, int, double, double, double, double, double, double, double, 
	int, double, double, double, double *, double, double, int, double *, double *, int, int, int); 
double **napMasEKTPVerNac(double, double, double, double, int, int, int, int, int, int, int, int, char *, int, double, double, 
	double, double *, int, double *, double *, int, double *, double *, int, int, double *, double *, double *, int, double); //double **napMasEKTPkviNac(double, double, double, double, int, int, char *, double, int, double, double, double *, int, double *, double *, int, double *, double *, int, double *, double *, double *, int);
double **NapMasVozd(char *);
double **zadNacZnaRazMas(double *, int, double *, double *, int, double *, double *, int, double, double, int, double, double, 
	double, double, int, int, char *);
double BolTochRasAlpha(int, double, double, double, double *, double *, double *, int, int, char *, double *, int, int, int, double *, 
	double *, double *, double *, int, double, double, double, double, double, double, double, int, double, double, double, double, 
	double *, double, double, int, double *, double *, int, int, int);
double opredKTPTKToch(double *, double *, double, int);
void progpoTem(int, double, char *, double, double, double, double *, double *, double *, double *, double *, double *, 
	double *, int, int, int, int, int, int, int, int, int, double *, double *, int, double *, double *, int, double *, 
	double *, double, double, int, int, int, int, int, int, int, double, double, double, double *, double, double, double, double);
double oprSrRazPor(int, int, double, int, int, int, int, int, int, char *);
double oprMaxRazPor(int, int, double, int, int, int, int, int, int, char *);
double oprKoefPogMod(int);
double *provnaOtrits(double *, double *, int, double *, double *, int, double *, int, double *, double *, double *, double *, 
	double, char *, double);
double oprKoefPogMod(int);
double *opredKoefOtrEdin(int, int, char *, double, double *, int, int);
double *reshnewtrafs(double ***, double ***, double *, double *, double *, double *, int, double *, double *, int, double *, double *, int, double, 
	double, double *, double, double *, int, double, double, double, char *, double, double *, double, double, int, double, double *, 
	double, double, double, double, double, int, int, double *, int, double, int, int, int, int, double, double, double, double);
double **RaschRTA(int, double, double, double, int, double, int, double, int, int, int, int, double *, double, int, double, double, 
	double, double, double, char *, double, double *, double *, int, double *, double *, int, double, double *, double *, int, double, double);
double *RasIzlSerStenNacSok(double *, double *, double *, double ***, double ***, double, double, double, double, double, double, double, int, double *, 
	double *, int, double *, double *, char *, double *, int, double, double, double, double, double, double, double, int, int, int,
	double *, int, double, int, int, int, int, int);
double ***zadrkt(int, char *, int, double *, double *, double, double, double, double, int);
double opredTempStenFragm(double, int, double *, double *, double *, double *, int, int, double, double, double, double);
double *KoefPoglRosselNac(double *, int, double, double, int, int);
double **opredTempLPSten(double, int, double *, double *, double *, double *, int, int, double, double, double, int, char *, double);
double *OprTempNewtRafsNach(int, double, double, double, double *, double *, double *, int, int, char *, double *, int, int, int, double *, 
	double *, double *, double *, int, double, double, double, double, double, double, double, int, double, double, double, double, 
	double *, double, double, double *, double *, int, double *, double *, int, double *, int, double, int, int); //double *gotVelKoefOtr(int, int, int, char *, int, int, int, int, int); double *gotVelKoefTeplRos(int, int, int, char *, int, int, int, int, int);
double **opredTempLPStenPerem(double, int, double *, double *, double *, double *, int, int, double, double, double, int, char *, double, double *);
double gotVeldkosp(int, int, int, int, int, int, int, double);
double opredTolTverChasObFigObPor(double, double, double, double, double, int, double, double);
double **opredTolVozdPros(int, int, int, int, int, int, int, int, int, int, int, int, char *, int, double, double *, double *, 
	double *, double *, double *, double *, double *, double *, double *, int, int, double, double);
double **gotZnachPeregor(int, int, int, int, char *, int, int);
double **napMasEKTPShaNac(double, double, double, double, int, int, int, int, char *, double, int, double, double, double, 
	double, double *, int, double *, double *, int, double *, double *, int, double *, double *, double *, int, double);
double **napMasEKTPitomNac(double, double, double, double, int, int, char *, int, double, double, double, double *, int, double *, 
	double *, int, double *, double *, int, double *, double *, double *, int, double);
double **napMasEKTPkviNac(double, double, double, double, int, int, char *, double, int, double, double, double *, int, 
	double *, double *, int, double *, double *, int, double *, double *, double *, int, double);
double oprDkoalMod();
//-----
void zadrktNac(int vybves)
{
int cemn=11, cem=cemn, cP=14; //число поколений = 12
double hf=1e0, dh=1e-3, ys0=3e1*1e-3, y0sha=65.0*1e-3;
double dkosp=hf; //dkosp - дополнительный коэффициент для КО с учетом пористой структуры
double dkoalmod=hf; //dkoal - определение (получение) КП из КО
int k=0;
int ks=10; //число стенок = 10
int vm=0;
int vrsh=3, vybsha=3, vystsha=1; //выбор марки шамота
int vmivmf=0, nnnvmimfv=vmivmf; //выбор метода измерений для фракции 2-0,7 мм: 0 - нестационарный, 1 - стационарный
int vyfv=0, nnnvyfv=vyfv; //выбор фракции: 0 - фракция 2-0,7 мм, 1 - фракция 8-4 мм, 2 - фракция 1,6-0,35 мм
int vyuv=1, nnnvyuv=vyuv; //выбор укладки: 1 - плоскопараллельная, 2 - вертикальная
int vysv=0, nnnvysv=vysv; //выбор состояния: 0 - исходное, 1 - после повторных измерений, 2 - после прокаливания при 1000 град С
int isrp=0, vpkf=0, vpmf=0; //выбор пористости малой фракции, крупной фракции и формы пор
int vybitom=2, nivybmar=vybitom, kvm=0; //выбор марки ИТОМ: 0 - 440, 1 - 620, 2 - 860, 3 - 1000
int vybkvi=6, nkvybmar=vybkvi; //выбор марки КВИ 4 - 400 ... 10 - 1000
double tocras=1e-7; 
if (!vybves) { ys0=y0sha; vm=vybsha; }
if (vybves==1) vm=0; 
if (vybves==2) vm=vybitom; 
if (vybves==3) vm=vybkvi; 
if ((vybves<0) || (vybves>3)) { cout << "Nepravilno ukazan nomer veschestva!" << endl; k=getchar(); exit(1); }
//---------
if (!vybves) { dkoalmod=0.65; dkosp=0.67; }
else if (vybves==1) { dkoalmod=0.44; dkoalmod=oprDkoalMod(); dkosp=1.3; }
else if (vybves==2) { dkoalmod=0.64; dkosp=1.3; }
else if (vybves==3) { dkoalmod=0.64; dkosp=1.3; } 
if (!vybves) initarrves(vm, vybves, cem, ys0, vyfv, vyuv, vmivmf, vpkf, vpmf, vysv, isrp, vystsha, vrsh, dh, tocras, ks, dkoalmod, dkosp, cP);
if (vybves==1) {
	int kvf=0, jvsv=0, qvmi=0, qvuv=0;
	int nnvyfv=2, nnvysv=3, nnvmivmf=3, nnvyuv=2;
	double *nvyfv=new double[nnvyfv], *nvysv=new double[nnvysv], *nvmivmf=new double[nnvmivmf];
	double *nvyuv=new double[nnvyuv]; 
	if ((!nvyfv) || (!nvysv) || (!nvmivmf) || (!nvyuv)) { cout << "No memory" << endl; k=getchar(); exit(1); } 
	k=0; nvyfv[k]=k; k++; nvyfv[k]=k; //фракции вермикулита
	k=0; nvysv[k]=k; k++; nvysv[k]=k; k++; nvysv[k]=k; //состояния вермикулита
	k=0; nvmivmf[k]=k; k++; nvmivmf[k]=k; k++; nvmivmf[k]=k; //стационарные методы измерений - 2019 и 2020
	k=0; nvyuv[k]=k+1; k++; nvyuv[k]=k+1; //укладка вермикулита
	for (kvf=nnnvyfv; kvf<nnvyfv; kvf++) {
		vyfv=nvyfv[kvf];
		for (jvsv=nnnvysv; jvsv<nnvysv; jvsv++) {
			vysv=nvysv[jvsv];
			if ((!vyfv) || (vyfv==2)) {
				for (qvmi=nnnvmimfv; qvmi<nnvmivmf; qvmi++) {
					vmivmf=nvmivmf[qvmi]; 
	initarrves(vm, vybves, cem, ys0, vyfv, vyuv, vmivmf, vpkf, vpmf, vysv, isrp, vystsha, vrsh, dh, tocras, ks, dkoalmod, dkosp, cP); 
				if (vysv) break; } }
			if (vyfv==1) {
				for (qvuv=(nnnvyuv-1); qvuv<nnvyuv; qvuv++) {
					vyuv=nvyuv[qvuv]; //cout << "vyuv = " << vyuv << "\t";
	initarrves(vm, vybves, cem, ys0, vyfv, vyuv, vmivmf, vpkf, vpmf, vysv, isrp, vystsha, vrsh, dh, tocras, ks, dkoalmod, dkosp, cP); } } } } 
	if (nvyfv) delete[]nvyfv; if (nvysv) delete[]nvysv; 
	if (nvmivmf) delete[]nvmivmf; if (nvyuv) delete[]nvyuv; }
if (vybves==2) {
	int knvybmar=3;
	for (kvm=nivybmar; kvm<=knvybmar; kvm++)
		initarrves(kvm, vybves, cem, ys0, vyfv, vyuv, vmivmf, vpkf, vpmf, vysv, isrp, vystsha, vrsh, dh, tocras, ks, dkoalmod, dkosp, cP); }
if (vybves==3) {
	int kvybmar=10;
	for (kvm=nkvybmar; kvm<=kvybmar; kvm++)
		initarrves(kvm, vybves, cem, ys0, vyfv, vyuv, vmivmf, vpkf, vpmf, vysv, isrp, vystsha, vrsh, dh, tocras, ks, dkoalmod, dkosp, cP); }
}
double **initNovMas(int koel, char *snm)
{ 
	int j=7; double *tx=NULL, *qx=NULL, *lx=NULL, *kx=NULL, *dtx=NULL, *oox=NULL, **mu=NULL, *temras=NULL, r=0.0;
	tx=new double[koel]; qx=new double[koel]; lx=new double[koel]; 
	kx=new double[koel]; dtx=new double[koel]; oox=new double[koel]; mu=new double*[j];  temras=new double[koel];
	if ((!tx) || (!qx) || (!lx) || (!kx) || (!dtx) || (!oox) || (!mu)) { cout << snm << endl; j=getchar(); exit(1); }
	for (j=0; j<koel; j++) { tx[j]=r; qx[j]=r; lx[j]=r; kx[j]=r; dtx[j]=r; oox[j]=r; temras[j]=r; }
	j=0; mu[j]=tx; j++; mu[j]=qx; j++; mu[j]=lx; j++; mu[j]=kx; j++; mu[j]=dtx; j++; mu[j]=oox; j++; mu[j]=temras; 
	return mu;
}
void initarrves(int vybmar, int vybves, int cem, double ys0, int vyfv, int vyuv, int vmimfv, int vpkf, int vpmf, int vysv, int isrp, 
	int vystsha, int vrsh, double dh, double tocras, int kost, double dkoalmod, double dkosp, int chPo)
{
	double te0=273.15, tnosc=3e2, dtosc=1e2, tnac=2e2+te0, dete=1e2, tnd=6e2, dtd=2e2, e=1e-6, t=e, *dpctp=NULL; //dpctp - доля площади перемычки, через которую переносится чистой теплопроводностью ТК
	double wmgo=0.0, wsio=0.0, walo=0.0, *tkusc=NULL, *dkoscm=NULL, *dkosct=NULL, *vsm=NULL, *poi=NULL; 
	int k=0, j=0, i=0, vybsha=vybmar, fl=1, dkoscl=6, dmkoosc=14, dmko=4, dmkvoz=28, dso=11; //вывод значений КТП ТК по различным моделям
	double smgo=0.0, salo=0.0, ssio=0.0, hf=1e0, pksvv=0.0, **mu=NULL, *kttk=NULL, *kusc=NULL;
	double por=0.0, *efftem=NULL, *tevoz=NULL, *ktpvoz=NULL, *Prvoz=NULL, *ektpm=NULL, *tgorm=NULL;
	double *tholm=NULL, *qobm=NULL, *stch=NULL, srp=0.0, s=0.0, vs=0.0, vss=0.0;
	double qobsh=0.0, nf=0.0, hk=ys0/2e0, dko=0.0, *htcm=NULL, tc=22.0+te0; 
	int vogf=0, vfqi=1, vgz=1;
	char *snm=new char[dso]; if (!snm) { cout << "No memory"; k=getchar(); exit(1); } 
	for (k=0; k<dso; k++) snm[k]='\0';
	k=0; snm[k]='N';  k++; snm[k]='o'; k++; snm[k]='_'; k++; snm[k]='m'; k++; snm[k]='e'; k++; 
		 snm[k]='m';  k++; snm[k]='o'; k++; snm[k]='r'; k++; snm[k]='y'; k++; snm[k]='!'; k++; snm[k]='\0';
	por=NovNapMas(vybves, vybmar, vmimfv, vyfv, vyuv, vysv, vpmf, vpkf); 
	dkoalmod=oprKoefPogMod(vybves); dkosp=hf; dko=dkoalmod*dkosp; 
	mu=NapMasVozd(snm); k=0; ktpvoz=mu[k]; k++; tevoz=mu[k]; k++; Prvoz=mu[k]; if (mu) delete[]mu; //КТП, число Прандтля и температура для них для воздуха
	efftem=new double[cem]; dkoscm=new double[dkoscl]; dkosct=new double[dkoscl]; 
	kusc=new double[dmkoosc]; tkusc=new double[dmkoosc]; //dkoscm - отклонение СЧ (или ПС) от литературных данных 
	if ((!dkoscm) || (!dkosct) || (!kusc) || (!tkusc) || (!efftem)) { cout << snm << endl; k=getchar(); exit(1); }
	mu=zadNacZnaRazMas(efftem, cem, dkoscm, dkosct, dkoscl, kusc, tkusc, dmkoosc, tnd, dtd, vybves, tnac, dete, tnosc, dtosc, vybmar, vyfv, snm);
	k=0; efftem=mu[k]; k++; dkoscm=mu[k]; k++; dkosct=mu[k]; k++; kusc=mu[k]; k++; tkusc=mu[k]; k++; vsm=mu[k];
	k=0; wmgo=vsm[k]; k++; wsio=vsm[k]; k++; walo=vsm[k]; if (mu) delete[]mu; 
	if (!vybves) { 
		mu=napMasEKTPShaNac(wmgo, wsio, walo, por, vybves, vybsha, vystsha, vrsh, snm, ys0, dmko, tnac, dete, 
		tnosc, dtosc, efftem, cem, tkusc, kusc, dmkoosc, dkosct, dkoscm, dkoscl, ktpvoz, tevoz, Prvoz, dmkvoz, dko); }
	if (vybves==1) { 
		mu=napMasEKTPVerNac(wmgo, wsio, walo, por, vybves, vyfv, vyuv, vmimfv, vysv, vpmf, vpkf, isrp, snm, dmko, ys0, 
		tnosc, dtosc, efftem, cem, tkusc, kusc, dmkoosc, dkosct, dkoscm, dkoscl, fl, ktpvoz, tevoz, Prvoz, dmkvoz, dko); }
	if (vybves==2) 
		mu=napMasEKTPitomNac(wmgo, wsio, walo, por, vybves, vybmar, snm, dmko, ys0, tnosc, dtosc, efftem,
		cem, tkusc, kusc, dmkoosc, dkosct, dkoscm, dkoscl, ktpvoz, tevoz, Prvoz, dmkvoz, dko); 
	if (vybves==3) 
		mu=napMasEKTPkviNac(wmgo, wsio, walo, por, vybves, vybmar, snm, ys0, dmko, tnosc, dtosc, efftem, cem, 
		tkusc, kusc, dmkoosc, dkosct, dkoscm, dkoscl, ktpvoz, tevoz, Prvoz, dmkvoz, dko);
	k=0; efftem=mu[k]; k++; ektpm=mu[k]; k++; tgorm=mu[k]; k++; tholm=mu[k]; 
	k++; qobm=mu[k];   k++; kttk=mu[k];  k++; stch=mu[k];  k++; poi=mu[k]; 
	if (mu) delete[]mu;
	j=0; nf=poi[j]; t=e; while (t<nf) { t=t+hf; j++; } cem=j; //for (k=0; k<cem; k++) printf("sch ( %d ) = %0.20lf\t", k, stch[k]);
	if (!vybves) {
	k=1; por=poi[k]; k++; vs=poi[k]; k++; vss=poi[k]; 
	j=0; t=e; while (t<vs) { j++; t=t+hf; } vybsha=j; j=0; t=e; while (t<vs) { j++; t=t+hf; } vystsha=j; }
	if (poi) delete[]poi; cout << "vf = " << vyfv << "\tvs = " << vysv << "\tvm = " << vmimfv << "\tvu = " << vyuv << "\tcemj = " << cem << "\t"; 
	progpoTem(cem, ys0, snm, por, dh, tocras, qobm, efftem, kttk, stch, ektpm, ktpvoz, tevoz, dmkvoz, vybves, vybmar, vyfv, 
	vysv, isrp, vpkf, vrsh, vystsha, kusc, tkusc, dmkoosc, dkoscm, dkosct, dkoscl, tholm, tgorm, dkoalmod, dkosp, kost, chPo,
	vyuv, vmimfv, vogf, vfqi, vpmf, wmgo, wsio, walo, Prvoz, tnosc, dtosc, tnac, dete);
	if (efftem) delete[]efftem; if (ektpm) delete[]ektpm; if (tgorm) delete[]tgorm; 
	if (tholm) delete[]tholm; if (qobm) delete[]qobm; if (kttk) delete[]kttk; 
	if (stch) delete[]stch; if (ktpvoz) delete[]ktpvoz; if (tevoz) delete[]tevoz; 
	if (Prvoz) delete[]Prvoz; if (tkusc) delete[]tkusc; if (dkoscm) delete[]dkoscm; 
	if (dkosct) delete[]dkosct; if (kusc) delete[]kusc; if (vsm) delete[]vsm; 
	if (snm) delete[]snm; 
}
double **zadNacZnaRazMas(double *efftem, int cem, double *dkoscm, double *dkosct, int dkoscl, double *kusc, double *tkusc, 
	int dmkoosc, double tnd, double dtd, int vybves, double tnac, double dete, double tnosc, double dtosc, int vybmar, int vyfv, char *snm)
{
	int k=6, nu=6, nw=6;
	double wmgo=0.0, wsio=0.0, walo=0.0, smgo=0.0, salo=0.0, ssio=0.0;
	double tm=0.0, ko=1e-2, hf=1e0, **mu=new double*[nu], *sov=new double[nw];
	if ((!mu) || (!sov)) { cout << snm << endl; k=getchar(); exit(1); }
	k=0; dkosct[k]=tnd; for (k=1; k<dkoscl; k++) dkosct[k]=dkosct[k-1]+dtd;
	if (!vybves) { //если выбран шамот
	k=0; dkoscm[k]=8.15;  k++; dkoscm[k]=6.87;  k++; dkoscm[k]=7.55; 
	k++; dkoscm[k]=14.84; k++; dkoscm[k]=20.09; k++; dkoscm[k]=26.3; }
	if (vybves==1) { //если выбран вермикулит
	k=0; dkoscm[k]=4.68;  k++; dkoscm[k]=4.69; k++; dkoscm[k]=5.65; 
	k++; dkoscm[k]=13.17; k++; dkoscm[k]=20.2; k++; dkoscm[k]=27.81;
		 ssio=36.8*ko; salo=13.2*ko; smgo=21.7*ko;
		 if (!vyfv) { ssio=45.8*ko; salo=14.0*ko; smgo=29.0*ko; } //2-0.7
		 else if (vyfv==1) { ssio=46.7*ko; salo=12.9*ko; smgo=28.2*ko; } //8-4
	wmgo=smgo; walo=salo; wsio=ssio; }
	if (vybves==2) { //если выбран ИТОМ
		if (!vybmar) { //ИТОМ-440
		salo=26e0*ko; smgo=22e0*ko; ssio=52e0*ko; //для трехкомпонентной смеси
		walo=23e0*ko; wmgo=19e0*ko; wsio=49e0*ko; //для многокомпонентной смеси
		k=0; dkoscm[k]=6.07;  k++; dkoscm[k]=5.36;  k++; dkoscm[k]=6.19; 
		k++; dkoscm[k]=13.48; k++; dkoscm[k]=19.93; k++; dkoscm[k]=27.69; } 
	if (vybmar==1) { //ИТОМ-620
		salo=29e0*ko; smgo=16e0*ko; ssio=55e0*ko;
		walo=26e0*ko; wmgo=15e0*ko; wsio=5e1*ko; 
		k=0; dkoscm[k]=6.41;  k++; dkoscm[k]=5.53; k++;  dkoscm[k]=6.32; 
		k++; dkoscm[k]=13.56; k++; dkoscm[k]=19.86; k++; dkoscm[k]=27.66; } 
	if (vybmar==2) { //ИТОМ-860
		salo=33e0*ko; smgo=11e0*ko; ssio=56e0*ko;
		walo=3e1*ko; wmgo=1e1*ko; wsio=52e0*ko; 
		k=0; dkoscm[k]=7.28;  k++; dkoscm[k]=6.1;   k++; dkoscm[k]=6.83; 
		k++; dkoscm[k]=14.02; k++; dkoscm[k]=19.86; k++; dkoscm[k]=27.22; } 
	if (vybmar==3) { //ИТОМ-1000
		salo=35e0*ko; smgo=9e0*ko; ssio=56e0*ko; walo=33e0*ko; wmgo=7e0*ko; wsio=53e0*ko; 
		k=0; dkoscm[k]=7.45; k++;  dkoscm[k]=6.24;  k++; dkoscm[k]=6.95; 
		k++; dkoscm[k]=14.15; k++; dkoscm[k]=19.89; k++; dkoscm[k]=27.09; } 
	if ((vybmar<0) || (vybmar>3)) { cout << "Takoy marki ITOM net!" << endl; k = getchar(); exit(1); } }
	if (vybves==3) { //если выбран КВИ
		if (vybmar==4) { walo=33e0*ko; wmgo=15e0*ko; wsio=52e0*ko; salo=25e0*ko; smgo=11e0*ko; ssio=4e1*ko; }
		if (vybmar==5) { walo=34e0*ko; wmgo=11e0*ko; wsio=54e0*ko; salo=28e0*ko; smgo=8e0*ko; ssio=44e0*ko; }
		if (vybmar==6) { walo=36e0*ko; wmgo=9e0*ko; wsio=55e0*ko; salo=3e1*ko; smgo=7e0*ko; ssio=45e0*ko; }
		if (vybmar==7) { walo=37e0*ko; wmgo=8e0*ko; wsio=55e0*ko; salo=31e0*ko; smgo=6e0*ko; ssio=45e0*ko; }
		if (vybmar==8) { walo=38e0*ko; wmgo=7e0*ko; wsio=55e0*ko; salo=31e0*ko; smgo=5e0*ko; ssio=45e0*ko; }
		if (vybmar==9) { walo=39e0*ko; wmgo=7e0*ko; wsio=55e0*ko; salo=32e0*ko; smgo=5e0*ko; ssio=45e0*ko; }
		if (vybmar==10) { walo=39e0*ko; wmgo=6e0*ko; wsio=55e0*ko; salo=32e0*ko; smgo=4e0*ko; ssio=45e0*ko; }
	k=0; dkoscm[k]=3e0;  k++; dkoscm[k]=6e0; k++; dkoscm[k]=7e0; 
	k++; dkoscm[k]=15e0; k++; dkoscm[k]=2e1; k++; dkoscm[k]=26e0;
		if (vybmar>6) { k=1; dkoscm[k]=7e0; } if (vybmar>9) { k=2; dkoscm[k]=8e0; } }
	walo=walo/(walo+wmgo+wsio); wmgo=wmgo/(walo+wmgo+wsio); wsio=wsio/(walo+wmgo+wsio);
	salo=salo/(salo+smgo+ssio); smgo=smgo/(salo+smgo+ssio); ssio=ssio/(salo+smgo+ssio);
	for (k=0; k<dkoscl; k++) dkoscm[k]=hf-dkoscm[k]*ko; //for (k=0; k<dkoscl; k++) cout << "T = " << dkosct[k] << "\tktp = " << dkoscm[k] << "\t"; cout << endl;
	//-----
	k=0; efftem[k]=tnac; for (k=1; k<cem; k++) efftem[k]=efftem[k-1]+dete;
	k=0; tkusc[k]=tnosc; for (k=1; k<dmkoosc; k++) tkusc[k]=tkusc[k-1]+dtosc;
	kusc=koefoslab(wmgo, wsio, walo, tkusc, dmkoosc, kusc, snm); //for (k=0; k<dmkoosc; k++) cout << "T_kut ( " << k << " ) = " << tkusc[k] << "\t" << kusc[k] << "\t"; cout << endl; k=getchar();
	k=0; sov[k]=wmgo; k++; sov[k]=wsio; k++; sov[k]=walo; k++; sov[k]=smgo; k++; sov[k]=salo; k++; sov[k]=ssio;
	k=0; mu[k]=efftem; k++; mu[k]=dkoscm; k++; mu[k]=dkosct; k++; mu[k]=kusc; k++; mu[k]=tkusc; k++; mu[k]=sov;
	return mu;
}
void progpoTem(int cem, double ys0, char *snm, double por, double dh, double tocras, double *qobm, double *efftem, double *ktptk, 
	double *stch, double *ektpm, double *ktpvo, double *vte, int dmkvoz, int vybves, int vybmar, int vyfv, 
	int vysv, int isrp, int vpkf, int vrsh, int vystsha, double *kusc, double *tkusc, int dmkoosc, double *dkoscm, 
	double *dkosct, int dkoscl, double *tholm, double *tgorm, double dkoal, double dkosp, int kost, int chP, int vyuv, 
	int vmimfv, int vogf, int vfqi, int vpmf, double wmgo, double wsio, double walo, double *Prvoz, double tnosc, 
	double dtosc, double tnac, double dete)
{
	FILE *fo=NULL;
	double hk=0.0, h=0.0, hvna=0.0, *tx=NULL, *qx=NULL, *lx=NULL, tem=0.0, thol=0.0, tgor=0.0, ektp=0.0;
	double *kx=NULL, *dtx=NULL, *oox=NULL, *temras=NULL, hf=1e0, e=1e-10, tsred=0.0, ot=0.0, ozl=hf, *ozlmin=new double[cem];
	double qobsh=0.0, *poi=NULL, mrp=0.0, htch=0.0, dkusce=0.0, dkoscee=0.0, *laok=new double[cem], laizv=0.0;
	double nzpa=1e0*hf, nzpb=nzpa, kzpa=1e1*hf, kzpb=kzpa, p1=hf*nzpa, p2=hf*nzpb, pka=kzpa*hf, pkb=kzpb*hf;
	double *dpctp=NULL, *htcm=NULL, dkov=dkoal*dkosp, koef=1e-3, yn=14.0*koef, yk=16.0*koef, ko=1e3;
	int k=0, vtn=0, vtk=cem-1, vt=0, knxt=0, nnxt=0, j=0, vp=1, dmko=0, q=0;
	int frko=1, frktr=1, vtko=0, kolPer=0, vgztp=1, vrdko=1, zp=1, zv=1, zt=1, qk=10; 
	double **mu=NULL, hvot=0.0, hvo=0.0, p=0.0, r=0.0, *Rmp=NULL, *rts=NULL, ka=0.0, kb=0.0, *ktr=NULL;
	double laiz=0.0, pab=0.0, hmi=1e-3, koe=1e-6, dpctpe=0.0;
	p=0.0; for (k=0; k<kost; k++) p=p+hf; 
	if (vybves==3) { yn=13.0*koef; yk=17.0*koef; } //yk=ys0; yn=0.0; 
	if (yk<hf) yk=yk*ko; if (yn<hf) yn=yn*ko; j=0; r=e; while (r<yn) { r=r+hf; j++; } nnxt=j; 
	j=nnxt; r=yn+e; while (r<yk) { r=r+hf; j++; } knxt=j; dmko=knxt-nnxt+1; 
	cout << "kx = " << knxt << "\tnx = " << nnxt << "\t" << "\tp = " << p << "\tvtn = " << vtn << "\tvybves = " << vybves << "\t";
	hvna=0.0; for (k=0; k<nnxt; k++) hvna=hvna+hf;
	rts=new double[dmko]; 
	if ((!rts) || (!ozlmin)) { cout << snm << endl; k=getchar(); exit(1); }
	for (k=vtn; k<cem; k++) { ozlmin[k]=1e3; laok[k]=0.0; } for (k=nnxt; k<knxt; k++) rts[k-nnxt]=0.0;
	mu=initNovMas(cem, snm);
k=0; tx=mu[k]; k++; qx=mu[k]; k++; lx=mu[k]; k++; kx=mu[k]; k++; dtx=mu[k]; k++; oox=mu[k]; k++; temras=mu[k]; 
if (mu) delete[]mu;
hvo=oprSrRazPor(vybves, vybmar, por, vyfv, vysv, isrp, vpkf, vrsh, vystsha, snm);
mrp=oprMaxRazPor(vybves, vybmar, por, vyfv, vysv, isrp, vpkf, vrsh, vystsha, snm); 
if (hvo>hmi) hvo=hvo*koe; if (mrp>hmi) mrp=mrp*koe; cout << "por = " << por << "\thvo = " << hvo << "\tmrp = " << mrp << "\t"; //k=getchar();
p1=nzpa*hf; p2=nzpb*hf; q=0; cout << "p1 = " << p1 << "\tp2 = " << p2 << "\t";
	while (q<qk) {
			for (vt=vtn; vt<=vtk; vt++) { //пробегаем по температуре 
				tgor=tgorm[vt]; thol=tholm[vt];
				ka=(thol-tgor)/ys0; kb=tgor; hk=yn; if (hk>hmi) hk=hk/ko;
					for (j=nnxt; j<=knxt; j++) { tem=ka*hk+kb; k=j-nnxt; rts[k]=tem; hk=hk+dh; cout << "rts = " << rts[k] << "\t"; 
					} tem=efftem[vt]; hk=ys0/2e0;
	poi=provnaOtrits(kusc, tkusc, dmkoosc, dkoscm, dkosct, dkoscl, efftem, cem, qobm, tholm, tgorm, ektpm, ys0, snm, tem); 
k=0; dkusce=poi[k]; k++; dkoscee=poi[k]; k++; qobsh=poi[k]; k++; thol=poi[k];	
k++; tgor=poi[k];	k++; ektp=poi[k];	 k++; tsred=poi[k]; if (poi) delete[]poi; 
cout << "qo = " << qobsh << "\tth = " << thol << "\ttg = " << tgor << "\tektp = " << ektp << "\tht = " << htch << "\tts = " << tsred << "\tvybves = " << vybves << endl;;
if (frko==1) Rmp=opredKoefOtrEdin(vybves, vybmar, snm, dkov, rts, dmko, vtko); //if (!frko) Rmp=gotVelKoefOtr(vt, vybves, dmko, snm, vtko, vyfv, vysv, vmimfv, vybmar); 
for (k=nnxt; k<=knxt; k++) cout << "Koef_Otr ( " << k-nnxt << " ) = " << Rmp[k-nnxt] << "\ttem = " << rts[k-nnxt] << "\t"; 
	if (zp) { 
		fo=fopen("Koefficient_Otrazheniya.txt", "at"); 
	if (!fo) { cout << "File is not open!" << endl; k=getchar(); exit(1); }
	fprintf(fo, "vybves = %d\tvybmar = %d\tvyuv = %d\tvysv = %d\tvyfv = %d\tvmi = %d\tvmrko = %d\vt = %d\n", vybves, vybmar, vyuv, vysv, vyfv, vmimfv, vtko, vt);
	fprintf(fo, "Temperatura\n");
	for (k=nnxt; k<knxt; k++) fprintf(fo, "%0.15lf\t", rts[k-nnxt]); 
	fprintf(fo, "Koefficient otrazheniya\n");
	for (k=nnxt; k<knxt; k++) fprintf(fo, "%0.20lf\t", Rmp[k-nnxt]); 
	fprintf(fo, "\n"); 
	fclose(fo); }
if (frktr==1) ktr=KoefPoglRosselNac(rts, dmko, dkosp, dkoal, vybves, vybmar); //if (!frktr) ktr=gotVelKoefTeplRos(vt, vybves, dmko, snm, vyfv, vyuv, vysv, vmimfv, vybmar); 
for (k=nnxt; k<=knxt; k++) cout << "ktr ( " << (k-nnxt) << " ) = " << ktr[k-nnxt] << "\ttem = " << rts[k-nnxt] << "\t"; cout << endl;
	if (zv) {
		fo=fopen("Koefficient_Rosselanda.txt", "at"); 
	if (!fo) { cout << "File is not open!" << endl; k=getchar(); exit(1); }
	fprintf(fo, "vybves = %d\tvybmar = %d\tvyuv = %d\tvysv = %d\tvyfv = %d\tvmi = %d\tvt = %d\n", vybves, vybmar, vyuv, vysv, vyfv, vmimfv, vt);
	fprintf(fo, "Temperatura\n"); 
	for (k=nnxt; k<knxt; k++) fprintf(fo, "%0.15lf\t", rts[k-nnxt]); 
	fprintf(fo, "Koefficient Rosselanda\n"); 
	for (k=nnxt; k<knxt; k++) fprintf(fo, "%0.20lf\t", ktr[k-nnxt]); 
	fprintf(fo, "\n"); 
	fclose(fo); }
vp=0;
if (!vgztp) mu=gotZnachPeregor(vyfv, vysv, vyuv, vmimfv, snm, vybves, vybmar);
else mu=opredTolVozdPros(vrsh, vystsha, kost, vybves, vybmar, vyfv, vysv, vyuv, isrp, vpmf, vmimfv, vpkf, snm, cem, por, ktpvo, 
	vte, efftem, ektpm, tgorm, tholm, qobm, ktptk, stch, dmkvoz, vogf, p1, p2); 
k=0; dpctp=mu[k]; k++; htcm=mu[k]; if (mu) delete[]mu; 
	if (!vrdko) dkosp=gotVeldkosp(vt, vybves, vybmar, vyfv, vysv, vyuv, vmimfv, (p1+p2)/2e0);
	else dkosp=RaschAlphaTvKar(vybves, vybmar, vystsha, vyfv, vysv, isrp, vpkf, por, snm, kost, dpctpe, vrsh, ktpvo, vte, efftem, cem, 
dmkvoz, kost, ktptk, dmkoosc, kusc, tkusc, dkoscm, dkosct, dkoscl, thol, tgor, ys0, dkoal, dkosp, hk, qobsh, chP, dkusce, dkoscee, tem, 
stch, p1, p2, vp, Rmp, rts, dmko, vfqi, vogf);
	printf("dkosp = %0.20lf\tp1 = %0.1lf\tp2 = %0.1lf\tvt = %d\n", dkosp, p1, p2, vt); 
	if (zt) { 
		fo=fopen("Dopolnitelniy_koefficient.txt", "at"); 
	fprintf(fo, "vybves = %d\tvybmar = %d\tvyuv = %d\tvysv = %d\tvyfv = %d\tvmi = %d\tp1 = %0.1lf\tp2 = %0.1lf\tdkosp = %0.20lf\tvrdko = %d\tvt = %d\n\n", 
		vybves, vybmar, vyuv, vysv, vyfv, vmimfv, p1, p2, dkosp, vrdko, vt);
	fclose(fo); }
	htch=opredKTPTKToch(htcm, efftem, tem, cem); if (htch<0.0) htch=0.0;
	dpctpe=opredKTPTKToch(dpctp, efftem, tem, cem); if (dpctpe<0.0) dpctpe=0.0; if (dpctpe>hf) dpctpe=hf;
poi=OprTempNewtRafsNach(kost, hvo, por, tem, ktpvo, vte, efftem, cem, dmkvoz, snm, ktptk, vybves, vybmar, dmkoosc, kusc,
tkusc, dkoscm, dkosct, dkoscl, thol, tgor, ys0, dkoal, dkosp, hk, qobsh, chP, htch, dkusce, dkoscee, dpctpe, stch, p1, p2, Rmp, 
rts, dmko, ektpm, ktr, dmko, dpctp, vfqi, dh, kolPer, vogf);
		k=0; ozl=fabs(poi[k]); k++; lx[vt]=poi[k]; k++; laizv=poi[k]; if (poi) delete[]poi; 
		if (ozl<ozlmin[vt]) { laok[vt]=lx[vt]; laiz=laizv; ozlmin[vt]=ozl; } 
		ot=fabs(ektp-lx[vt])*1e2/ektp; 
		qx[vt]=qobsh; kx[vt]=hk; 
		cout << "hk = " << hk << "\tte = " << tem  << "\tvt = " << vt << "\tqx = " << qx[vt] << "\tp1 = " << p1 << "\tp2 = " << p2 << 
					"\tlx = " << lx[vt] << "\tosh = " << ot << endl;
				//-----
	fo=fopen("Koefficient_Teploprovodnosti.txt", "at"); 
	fprintf(fo, "При температуре T = %0.2lf\tp1 = %0.1lf\tp2 = %0.1lf\tdkosp = %0.3lf\tvyfv = %d\tvysv = %d\tvyuv = %d\tvmimf = %d\tvybves = %d\tvybmar = %d\t", 
	efftem[vt], p1, p2, dkosp, vyfv, vysv, vyuv, vmimfv, vybves, vybmar); 
	fprintf(fo, "Ошибка: %lf\tКТП: %lf\tВыбор КО = %d\tВыбор фигуры = %d\tВыбор метода = %d\tЭКТП = %lf\n", 
	ot, lx[vt], vtko, vogf, vfqi, ektpm[vt]);
	fclose(fo); 
	if (ktr) delete[]ktr; if (Rmp) delete[]Rmp; 
	if (dpctp) delete[]dpctp; if (htcm) delete[]htcm; }
			p1=p1+hf; p2=p2+hf; q++; cout << "lam_ras = " << laok[vt] << "\totn_otkl = " << ozlmin[vt] << endl; 
	if ((p1>(pka-e)) || (p2>(pkb-e))) break; }
	if (tx) delete[]tx; if (qx) delete[]qx; if (lx) delete[]lx; if (kx) delete[]kx; 
	if (dtx) delete[]dtx; if (oox) delete[]oox; if (temras) delete[]temras; if (rts) delete[]rts; 
	if (laok) delete[]laok; if (ozlmin) delete[]ozlmin;
}
double *provnaOtrits(double *kusc, double *tkusc, int dmkoosc, double *dkoscm, double *dkosct, int dkoscl, double *efftem, 
	int cem, double *qobm, double *tholm, double *tgorm, double *ektpm, double tolob, char *snm, double tc)
{
	int k=0, ce=7; 
	double hf=1e0, *vm=new double[ce], tsred=0.0;
	if (!vm) { cout << snm << endl; k=getchar(); exit(1); }
	double dkusce=opredKTPTKToch(kusc, tkusc, tc, dmkoosc); if (dkusce>hf) dkusce=hf; if (dkusce<0.0) dkusce=0.0;
	double dkoscee=opredKTPTKToch(dkoscm, dkosct, tc, dkoscl); if (dkoscee>hf) dkoscee=hf; if (dkoscee<0.0) dkoscee=0.0;
	double qobsh=hf, thol=hf, tgor=hf, ektp=hf, t=hf;
	qobsh=opredKTPTKToch(qobm, efftem, tc, cem); 
	thol=opredKTPTKToch(tholm, efftem, tc, cem); 
	tgor=opredKTPTKToch(tgorm, efftem, tc, cem); 
	ektp=opredKTPTKToch(ektpm, efftem, tc, cem); 
	tsred=opredKTPTKToch(efftem, efftem, tc, cem);
	if (ektp>0.0) {
		if (qobsh>0.0) {
			if (thol>0.0) {
				if (tgor<0.0) { t=qobsh*tolob/ektp; tgor=thol+t; } } 
			else { t=qobsh*tolob/ektp; if (tgor>t) thol=tgor-t; } } 
		else { if ((tgor>0.0) && (thol>0.0) && (tgor>thol)) qobsh=ektp*(tgor-thol)/tolob; else qobsh=0.0; } } 
	else { if ((tgor>0.0) && (thol>0.0) && (tgor>thol) && (qobsh>0.0)) ektp=qobsh*tolob/(tgor-thol); else ektp=0.0; }
	if (tsred<0.0) tsred=(tgor+thol)/2e0; 
	k=0; vm[k]=dkusce; k++; vm[k]=dkoscee; k++; vm[k]=qobsh; k++; vm[k]=thol;   
	k++; vm[k]=tgor;   k++; vm[k]=ektp;    k++; vm[k]=tsred;
	return vm;
}
double *OprTempNewtRafsNach(int kost, double hvp, double por, double tc, double *ktpvo, double *vte, double *ete, int cem, int dmkvoz, 
	char *snm, double *ktptk, int vybves, int vybmar, int dmkoosc, double *kusc, double *tkusc, double *dkoscm, double *dkosct, 
	int dkoscl, double thol, double tgor, double y0, double dkoal, double dkosp, double hk, double qobsh, int kolPok, double htk, 
	double dkusce, double dkoscee, double dpctpe, double *stch, double p1, double p2, double *Refm, double *rts, int dmrts, 
	double *laefm, double *ktr, int dmktr, double *dpctp, int vfqi, double dx, int kolPer, int vogf)
{
	int j=0, nk=j, l=0, k=0, zf=0, iz=1, v=0, c=0, w=0, kevm=5, vp=0, n=2*kost, moi=-1; //w - расчет температуры, iz - изменение
	double ***protv=NULL, ***prots=NULL, hf=1e0, hlr0=hf, hrl0=0.0, vyzn=0.0, ra=fabs(hlr0-hrl0), tol=0.0, *ptr=NULL;
	double *Ts=NULL, *Tss=NULL, **mu=NULL, **mauk=NULL, *po=NULL, *sislr=NULL, *sisrl=NULL, *qxi=NULL; 
	double d=0.0, ka=d, kb=d, e=1e-9, mo=-hf, tx=tgor-fabs(tgor-thol)*hk/y0, tora=1e-4, tma=d, tmi=d, ts=0.0;
	double *hrl=NULL, *hlr=NULL, *te=new double[n], r=0.0, *uvp=NULL, lambda=0.0, tsred=0.0, *ktpr=new double[n];
	if ((!te) || (!ktpr)) { cout << snm << endl; k=getchar(); exit(1); }
	for (j=0; j<kost; j++) d=d+hf; tol = d*htk + (d - hf)*hvp;
	if (vybves>=4) { cout << "Net takogo veschestva!" << endl; k = getchar(); exit(1); } 
	mauk=RaschRTA(kost, htk, ka, kb, iz, hvp, v, tc, c, w, vybves, cem, ete, hk, vybmar, dkusce, dkoscee, thol, tgor, y0, snm, dkosp, 
		kusc, tkusc, dmkoosc, dkoscm, dkosct, dkoscl, dkoal, Refm, rts, dmrts, p1*hvp, p2*hvp); 
	double *Refa=NULL, *Traa=NULL, *Aba=NULL, *Ref=NULL, *Tra=NULL, *Ab=NULL, *Reft=NULL, *Trat=NULL, *Abt=NULL;
	k=0; Refa=mauk[k]; k++; Traa=mauk[k]; k++; Aba=mauk[k];  k++; Ref=mauk[k];
	k++; Tra=mauk[k];  k++; Ab=mauk[k];   k++; Reft=mauk[k]; k++; Trat=mauk[k]; k++; Abt=mauk[k]; k++; nk=k; 
	k=0; protv=zadrkt(kost, snm, kolPok, Ref, Tra, htk, hvp, p1*hvp, p2*hvp, k); //vyv - выбор вещества: 0 - шамот, 1 - вермикулит
	k++; prots=zadrkt(kost, snm, kolPok, Reft, Trat, htk, hvp, p1*hvp, p2*hvp, k); 
	mu=opredTempLPStenPerem(tc, n, ktpvo, vte, ete, ktptk, cem, dmkvoz, htk, hvp, qobsh, kost, snm, mo, dpctp); 
	k=0; Ts=mu[k]; k++; Tss=mu[k]; if (mu) delete[]mu; 
	j=0; for (k=0; k<n; k++) { if (!(k%2)) te[k]=Ts[j]; else { te[k]=Tss[j]; j++; } } 
	if (Ts) delete[]Ts; if (Tss) delete[]Tss; 
	k=0; tmi=te[k]; tma=te[k];
	for (k=0; k<n; k++) { if (tmi>te[k]) tmi=te[k]; if (tma<te[k]) tmi=te[k]; } 
	hlr0=opredKTPTKToch(ktr, rts, tc, dmktr); if (hlr0<0.0) hlr0=0.0; 
	hlr0=fabs(hlr0*opredTempStenFragm(tc, n, ktpvo, vte, ete, ktptk, cem, dmkvoz, htk, hvp, qobsh, mo)/tol); 
	r=0.0; k=0; vp=2; 
	hlr=RasIzlSerStenNacSok(Ref, Tra, Ab, protv, prots, tc, hlr0, hrl0, por, dpctpe, htk, hvp, kost, ete, ktptk, cem, 
		ktpvo, vte, snm, stch, dmkvoz, qobsh, hk, p1, p2, tx, dkosp, dkoal, vybves, vybmar, vp, ktr, dmktr, r, k, moi, 
		vfqi, kolPok, vogf); 
	r=0.0; k=0; vp=3;
	hrl=RasIzlSerStenNacSok(Ref, Tra, Ab, protv, prots, tc, hlr0, hrl0, por, dpctpe, htk, hvp, kost, ete, ktptk, cem, 
		ktpvo, vte, snm, stch, dmkvoz, qobsh, hk, p1, p2, tx, dkosp, dkoal, vybves, vybmar, vp, ktr, dmktr, r, k, moi, 
		vfqi, kolPok, vogf);
	uvp=reshnewtrafs(protv, prots, hlr, hrl, Tra, Ab, kost, ktptk, ete, cem, ktpvo, vte, dmkvoz, htk, hvp, te, tora, 
		laefm, kevm, tma, tmi, y0, snm, qobsh, Ref, hlr0, hrl0, vybves, por, stch, hk, p1, p2, dkosp, dkoal, vybmar, 
		vp, ktr, dmktr, dpctpe, vfqi, kolPok, kolPer, vogf, tol, tc, thol, tgor); 
	v=4; for (w=0; w<v; w++) { 
		for (k=0; k<kost; k++) { 
		ptr=prots[w][k]; if (ptr) delete[]ptr; 
		ptr=protv[w][k]; if (ptr) delete[]ptr; } }
	if (protv) delete[]protv; if (prots) delete[]prots; 
	for (k=0; k<nk; k++) { po=mauk[k]; if (po) delete[]po; } if (mauk) delete[]mauk; 
	if (hlr) delete[]hlr; if (hrl) delete[]hrl;
	if (te) delete[]te; if (ktpr) delete[]ktpr; 
	return uvp;
}
double oprDkoalMod()
{
	int cems=3*4, k=0, cemc=3*3, cemp=cems;
	double *kpps=new double[cems], *krps=new double[cems], *kprs=new double[cems];
	double *kppp=new double[cemp], *krpp=new double[cemp], *kprp=new double[cemp];
	double *kppc=new double[cemc], *krpc=new double[cemc], *kprc=new double[cemc];
	double rs=0.0, ss=0.0, rp=0.0, sp=0.0, rc=0.0, sc=0.0;
	if ((!kpps) || (!kprs) || (!krps) || (!kppp) || (!kprp) || (!krpp) || (!kppc) || (!kprc) || (!krpc)) 
	{ cout << "No memory!" << endl; k=getchar(); exit(1); }
k=0; kprs[k]=3923.0; k++; kprs[k]=4038.0; k++; kprs[k]=5304.0;
k++; kprs[k]=2539.0; k++; kprs[k]=1944.0; k++; kprs[k]=1463.0;
k++; kprs[k]=2492.0; k++; kprs[k]=1755.0; k++; kprs[k]=1235.0;
k++; kprs[k]=2380.0; k++; kprs[k]=1537.0; k++; kprs[k]=1197.0;
//-----
k=0; kpps[k]=31288.0; k++; kpps[k]=34236.0; k++; kpps[k]=27945.0;
k++; kpps[k]=17847.0; k++; kpps[k]=17307.0; k++; kpps[k]=16748.0;
k++; kpps[k]=10270.0; k++; kpps[k]=11025.0; k++; kpps[k]=9793.0;
k++; kpps[k]=10229.0; k++; kpps[k]=7652.0; k++; kpps[k]=8170.0;
//-----
k=0; krps[k]=39517.0; k++; krps[k]=29924.0; k++; krps[k]=20562.0;
k++; krps[k]=30504.0; k++; krps[k]=26066.0; k++; krps[k]=22953.0;
k++; krps[k]=13790.0; k++; krps[k]=15047.0; k++; krps[k]=14918.0;
k++; krps[k]=15844.0; k++; krps[k]=10855.0; k++; krps[k]=11533.0;
//-----
rs=0.0; ss=0.0; for (k=0; k<cems; k++) { rs=rs+kpps[k]; ss=ss+kpps[k]+krps[k]; }
rs=rs/ss; //cout << "rs = " << rs << endl;
//-----
k=0; kppp[k]=33691.0; k++; kppp[k]=33212.0; k++; kppp[k]=37391.0;
k++; kppp[k]=22772.0; k++; kppp[k]=24099.0; k++; kppp[k]=25633.0;
k++; kppp[k]=11632.0; k++; kppp[k]=13803.0; k++; kppp[k]=13405.0;
k++; kppp[k]=13456.0; k++; kppp[k]=9617.0;  k++; kppp[k]=10491.0;
//-----
k=0; krpp[k]=28209.0; k++; krpp[k]=24030.0; k++; krpp[k]=22468.0;
k++; krpp[k]=19985.0; k++; krpp[k]=19284.0; k++; krpp[k]=18137.0;
k++; krpp[k]=12406.0; k++; krpp[k]=12299.0; k++; krpp[k]=11334.0;
k++; krpp[k]=12657.0; k++; krpp[k]=8895.0;  k++; krpp[k]=9219.0;
//-----
k=0; kprp[k]=12769.0; k++; kprp[k]=10477.0; k++; kprp[k]=10811.0;
k++; kprp[k]=10015.0; k++; kprp[k]=10531.0; k++; kprp[k]=9893.0;
k++; kprp[k]=6368.0;  k++; kprp[k]=7427.0;  k++; kprp[k]=6909.0;
k++; kprp[k]=7576.0;  k++; kprp[k]=5266.0;  k++; kprp[k]=5785.0;
//-----
rp=0.0; sp=0.0; for (k=0; k<cemp; k++) { rp=rp+kppp[k]; sp=sp+kppp[k]+krpp[k]; }
rp=rp/sp; //cout << "rp = " << rp << endl;
//-----
k=0; kprc[k]=2869.0; k++; kprc[k]=1867.0; k++; kprc[k]=1383.0;
k++; kprc[k]=4691.0; k++; kprc[k]=2815.0; k++; kprc[k]=2182.0;
k++; kprc[k]=2509.0; k++; kprc[k]=1765.0; k++; kprc[k]=1215.0;
//-----
k=0; kppc[k]=21116.0; k++; kppc[k]=20128.0; k++; kppc[k]=22313.0;
k++; kppc[k]=13900.0; k++; kppc[k]=10572.0; k++; kppc[k]=11753.0;
k++; kppc[k]=9775.0;  k++; kppc[k]=7559.0;  k++; kppc[k]=6916.0;
//-----
k=0; krpc[k]=21616.0; k++; krpc[k]=21361.0; k++; krpc[k]=18674.0;
k++; krpc[k]=9093.0;  k++; krpc[k]=14297.0; k++; krpc[k]=12083.0;
k++; krpc[k]=18883.0; k++; krpc[k]=10308.0; k++; krpc[k]=12087.0;
//-----
rc=0.0; sc=0.0; for (k=0; k<cemc; k++) { rc=rc+kppc[k]; sc=sc+kppc[k]+krpc[k]; }
rc=rc/sc; //cout << "rc = " << rc << endl;
if (krps) delete[]krps; if (kprs) delete[]kprs; if (kpps) delete[]kpps;
if (krpp) delete[]krpp; if (kprp) delete[]kprp; if (kppp) delete[]kppp;
if (krpc) delete[]krpc; if (kprc) delete[]kprc; if (kppc) delete[]kppc;
rc=(rc+rp+rs)/3e0;
return rc;
}