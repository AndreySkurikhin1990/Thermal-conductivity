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
double **initarrves(int, int, int, double, int, int, int, int, int, int, int, int, int, double, double, int, double, double, int);
double *koefoslab(double, double, double, double *, int, double *, char *);
double NovNapMas(int, int, int, int, int, int, int, int); 
double RaschAlphaTvKar(int, int, int, int, int, int, int, double, char *, int, double, int, double *, double *, double *, int, 
	int, int, double *, int, double *, double *, double *, double *, int, double, double, double, double, double, double, double, 
	int, double, double, double, double *, double, double, int, double *, double *, int, int); //double **napMasEKTPShaNac(double, double, double, int, int, int, int, char *, double, int, double, double, double, double, double *, int, double *, double *, int, double *, double *, int, double *, double *, double *, int);//double **napMasEKTPitomNac(double, double, double, int, int, char *, int, double, double, double, double *, int, double *, double *, int, double *, double *, int, double *, double *, double *, int);
double **napMasEKTPVerNac(double, double, double, double, int, int, int, int, int, int, int, int, char *, int, double, double, 
	double, double *, int, double *, double *, int, double *, double *, int, int, double *, double *, double *, int, double); //double **napMasEKTPkviNac(double, double, double, double, int, int, char *, double, int, double, double, double *, int, double *, double *, int, double *, double *, int, double *, double *, double *, int);
double **NapMasVozd(char *);
double **zadNacZnaRazMas(double *, int, double *, double *, int, double *, double *, int, double, double, int, double, double, 
	double, double, int, int, char *);
double BolTochRasAlpha(int, double, double, double, double *, double *, double *, int, int, char *, double *, int, int, int, double *, 
	double *, double *, double *, int, double, double, double, double, double, double, double, int, double, double, double, double, 
	double *, double, double, int, double *, double *, int, int);
double opredKTPTKToch(double *, double *, double, int);
void progpoTem(int, double, char *, double, double, double, double *, double *, double *, double *, double *, double *, double *, double *, 
	int, int, int, int, int, int, int, int, int, double *, double *, int, double *, double *, int, double *, double *, double, double, 
	int, int, double *);
double oprSrRazPor(int, int, double, int, int, int, int, int, int, char *);
double oprMaxRazPor(int, int, double, int, int, int, int, int, int, char *);
void vyvodfile(double *, int, int, double, double, double, int, int);
char **napolStrok(int, int, char);
char **napNazvFile(int, int, char *, char);
double **opredTolVozdPros(int, int, int, int, int, int, int, int, int, int, int, int, char *, double, double, double, int, double, 
	double, double, int, double *, double *, int, double *, double *, int, double, double *, double *, double *, double, double, 
	double *, double *, double *, double *, double *, double *, double *, int);
double oprKoefPogMod(int);
double *provnaOtrits(double *, double *, int, double *, double *, int, double *, double *, int, double *, double *, double *, 
	double *, double, char *, double, double *);
double oprKoefPogMod(int);
double opredTolTverChasObFig(double, double, double, double, double);
double *opredKoefOtrEdin(int, int, char *, double, double *, int, int);
double *reshnewtrafs(double ***, double ***, double *, double *, double *, double *, int, double *, double *, int, double *, double *, int, double, 
	double, double *, double, double *, int, double, double, double, char *, double, double *, double, double, int, double, double *, 
	double, double, double, double, double, int, int, double *, int, double *, int, int, int);
double **RaschRTA(int, double, double, double, int, double, int, double, int, int, int, int, double *, double, int, double, double, 
	double, double, double, char *, double, double *, double *, int, double *, double *, int, double, double *, double *, int, double, double);
double *RasIzlSerStenNacSok(double *, double *, double *, double ***, double ***, double, double, double, double, double, double, double, int, double *, 
	double *, int, double *, double *, char *, double *, int, double, double, double, double, double, double, double, int, int, int,
	double *, int, double, int, int, int, int);
double ***zadrkt(int, char *, int, double *, double *, double, double, double, double, int);
double opredTempStenFragm(double, int, double *, double *, double *, double *, int, int, double, double, double, double);
double *KoefPoglRosselNac(double *, int, double, double, int, int);
double **opredTempLPSten(double, int, double *, double *, double *, double *, int, int, double, double, double, int, char *, double);
double OprTempNewtRafsNach(int, double, double, double, double *, double *, double *, int, int, char *, double *, int, int, int, double *, 
	double *, double *, double *, int, double, double, double, double, double, double, double, int, double, double, double, double, 
	double *, double, double, double *, double *, int, double *, double *, int, double *, int, double, int);
double *gotVelKoefOtr(int, int, int, char *, int);
double *gotVelKoefTeplRos(int, int, int, char *);
//-----
void zadrktNac(int vybves)
{
int cemn=11, cem=cemn, cP=12;
double por=0.0, hf=1e0, dh=1e-3, y0=3e1*1e-3, y0sha=65.0*1e-3;
double dkosp=hf, *kttk=NULL, **mu=NULL; //dkosp - дополнительный коэффициент для КО с учетом пористой структуры
double dkoalmod=hf; //dkoal - определение (получение) КП из КО
int j=0, k=0, q=0, f=6, qk=0, cei=0, N=13, ks=10, vm=0;
int vrsh=3, vybsha=3, vystsha=1; //выбор марки шамота
int vmivmf=0; //выбор метода измерений для фракции 2-0,7 мм: 0 - нестационарный, 1 - стационарный
int vyfv=0; //выбор фракции: 0 - фракция 2-0,7 мм, 1 - фракция 8-4 мм, 2 - фракция 1,6-0,35 мм
int vyuv=2; //выбор укладки: 1 - плоскопараллельная, 2 - вертикальная
int vysv=0; //выбор состояния: 0 - исходное, 1 - после повторных измерений, 2 - после прокаливания при 1000 град С
int isrp=0, vpkf=0, vpmf=0; //выбор пористости малой фракции, крупной фракции и формы пор
int vybitom=0; //выбор марки ИТОМ: 0 - 440, 1 - 620, 2 - 860, 3 - 1000
int vybkvi=6; //выбор марки КВИ 4 - 400 ... 10 - 1000
double tocras=1e-7, d=tocras, ys0=y0, r=0.0; if (!vybves) ys0=y0sha;
double dpct=0.314; //dpct - доля площади перемычки, через которую переносится чистой теплопроводностью ТК
if (!vybves) vm=vybsha; 
else if (vybves==1) vm=0; 
else if (vybves==2) vm=vybitom; 
else if (vybves==3) vm=vybkvi; 
else { cout << "Nepravilno ukazan nomer veschestva!" << endl; k=getchar(); exit(1); }
//---------
double t=0.0, nf=t, hn=t, ko=t, e=1e-8, x=t, g=t, hx=t, p=t; 
d=0.0; for (j=0; j<ks; j++) d=d+hf; 
//--------------
if (!vybves) { dkoalmod=0.65; dkosp=0.67; }
else if (vybves==1) { dkoalmod=0.44; dkosp=1.3; }
else if (vybves==2) { dkoalmod=0.64; dkosp=1.3; }
else if (vybves==3) { dkoalmod=0.64; dkosp=1.3; } //cout << "dkoal = " << dkoalmod << "\t"; 
mu=initarrves(vm, vybves, cem, ys0, vyfv, vyuv, vmivmf, vpkf, vpmf, vysv, isrp, vystsha, vrsh, dh, tocras, ks, dkoalmod, dkosp, cP);
}
double **initNovMas(int koel, char *snm)
{ int j=7; double *tx=NULL, *qx=NULL, *lx=NULL, *kx=NULL, *dtx=NULL, *oox=NULL, **mu=NULL, *temras=NULL, r=0.0;
	tx=new double[koel]; qx=new double[koel]; lx=new double[koel]; 
	kx=new double[koel]; dtx=new double[koel]; oox=new double[koel]; mu=new double*[j];  temras=new double[koel];
	if ((!tx) || (!qx) || (!lx) || (!kx) || (!dtx) || (!oox) || (!mu)) { cout << snm << endl; j=getchar(); exit(1); }
	for (j=0; j<koel; j++) { tx[j]=r; qx[j]=r; lx[j]=r; kx[j]=r; dtx[j]=r; oox[j]=r; temras[j]=r; }
	j=0; mu[j]=tx; j++; mu[j]=qx; j++; mu[j]=lx; j++; mu[j]=kx; j++; mu[j]=dtx; j++; mu[j]=oox; j++; mu[j]=temras; 
	return mu;
}
double **initarrves(int vybmar, int vybves, int cem, double ys0, int vyfv, int vyuv, int vmimfv, int vpkf, int vpmf, int vysv, int isrp, 
	int vystsha, int vrsh, double dh, double tocras, int kost, double dkoalmod, double dkosp, int chPo)
{
	double te0=273.15, tnosc=3e2, dtosc=1e2, tnac=2e2+te0, dete=1e2, tnd=6e2, dtd=2e2, e=1e-6, t=e, *dpctp=NULL;
	double wmgo=0.0, wsio=0.0, walo=0.0, *tkusc=NULL, *dkoscm=NULL, *dkosct=NULL, *vsm=NULL, *poi=NULL; 
	int k=0, j=0, i=0, vybsha=vybmar, fl=1, dkoscl=6, dmkoosc=14, dmko=4, dmkvoz=28, dso=11; //вывод значений КТП ТК по различным моделям
	double smgo=0.0, salo=0.0, ssio=0.0, hf=1e0, pksvv=0.0, **mu=NULL, *kttk=NULL, *kusc=NULL;
	double por=0.0, *efftem=NULL, *tevoz=NULL, *ktpvoz=NULL, *Prvoz=NULL, *ektpm=NULL, *tgorm=NULL;
	double *tholm=NULL, *qobm=NULL, *stch=NULL, srp=0.0, s=0.0;
	double qobsh=0.0, nf=0.0, hk=ys0/2e0, dko=0.0, *htcm=NULL, tc=22.0+te0;  
	char *snm=new char[dso]; if (!snm) { cout << "No memory"; k=getchar(); exit(1); } 
	for (k=0; k<dso; k++) snm[k]='\0';
	k=0; snm[k]='N';  k++; snm[k]='o'; k++; snm[k]='_'; k++; snm[k]='m'; k++; snm[k]='e'; k++; 
		 snm[k]='m';  k++; snm[k]='o'; k++; snm[k]='r'; k++; snm[k]='y'; k++; snm[k]='!'; k++; snm[k]='\0';
	por=NovNapMas(vybves, vybmar, vmimfv, vyfv, vyuv, vysv, vpmf, vpkf); //cout << "poris = " << por << "\t";
	dkoalmod=oprKoefPogMod(vybves); dkosp=hf; dko=dkoalmod*dkosp; //cout << "dkoal = " << dkoalmod << "\t"; if (dko<0.0) dko=0.0; if (dko>hf) dko=hf;
	mu=NapMasVozd(snm); k=0; ktpvoz=mu[k]; k++; tevoz=mu[k]; k++; Prvoz=mu[k]; if (mu) delete[]mu; //КТП, число Прандтля и температура для них для воздуха
	efftem=new double[cem]; dkoscm=new double[dkoscl]; dkosct=new double[dkoscl]; 
	kusc=new double[dmkoosc]; tkusc=new double[dmkoosc]; //dkoscm - отклонение СЧ (или ПС) от литературных данных 
	if ((!dkoscm) || (!dkosct) || (!kusc) || (!tkusc) || (!efftem)) { cout << snm << endl; k=getchar(); exit(1); }
	mu=zadNacZnaRazMas(efftem, cem, dkoscm, dkosct, dkoscl, kusc, tkusc, dmkoosc, tnd, dtd, vybves, tnac, dete, tnosc, dtosc, vybmar, vyfv, snm);
	k=0; efftem=mu[k]; k++; dkoscm=mu[k]; k++; dkosct=mu[k]; k++; kusc=mu[k]; k++; tkusc=mu[k]; k++; vsm=mu[k];
	k=0; wmgo=vsm[k]; k++; wsio=vsm[k]; k++; walo=vsm[k]; if (mu) delete[]mu; //if (!vybves) mu=napMasEKTPShaNac(wmgo, wsio, walo, por, vybves, vybsha, vystsha, vrsh, snm, ys0, dmko, tnac, dete, tnosc, dtosc, efftem, cem, tkusc, kusc, dmkoosc, dkosct, dkoscm, dkoscl, ktpvoz, tevoz, Prvoz, dmkvoz); //k=3; vsm=mu[k]; j=0; vs=vsm[j]; j++; vss=vsm[j]; j++; por=vsm[j]; //j=0; t=e; while (t<vs) { j++; t=t+hf; } vybsha=j; j=0; t=e; while (t<vs) { j++; t=t+hf; } vystsha=j; 
	if (vybves==1) mu=napMasEKTPVerNac(wmgo, wsio, walo, por, vybves, vyfv, vyuv, vmimfv, vysv, vpmf, vpkf, isrp, snm, dmko, ys0, 
		tnosc, dtosc, efftem, cem, tkusc, kusc, dmkoosc, dkosct, dkoscm, dkoscl, fl, ktpvoz, tevoz, Prvoz, dmkvoz, dko); //if (vybves==2) mu=napMasEKTPitomNac(wmgo, wsio, walo, por, vybves, vybmar, snm, dmko, ys0, tnosc, dtosc, efftem, cem, tkusc, kusc, dmkoosc, dkosct, dkoscm, dkoscl, ktpvoz, tevoz, Prvoz, dmkvoz); //if (vybves==3) mu=napMasEKTPkviNac(wmgo, wsio, walo, por, vybves, vybmar, snm, ys0, dmko, tnosc, dtosc, efftem, cem, tkusc, kusc, dmkoosc, dkosct, dkoscm, dkoscl, ktpvoz, tevoz, Prvoz, dmkvoz);
	if ((vybves>3) || (vybves<0)) { cout << "Oshibka v vybore veschestva!" << endl; k=getchar(); exit(1); }
	k=0; efftem=mu[k]; k++; ektpm=mu[k]; k++; tgorm=mu[k]; k++; tholm=mu[k]; k++; qobm=mu[k]; 
	k++; poi=mu[k];    k++; kttk=mu[k];  k++; stch=mu[k]; if (mu) delete[]mu;
	k=0; nf=poi[k]; t=e; while (t<nf) { t=t+hf; k++; } cem=k; if (poi) delete[]poi; //cout << "nf = " << nf << "\t";
	mu=opredTolVozdPros(vrsh, vystsha, kost, vybves, vybmar, vyfv, vysv, vyuv, isrp, vpmf, vmimfv, vpkf, snm, wmgo, wsio, walo, dmko, 
		ys0, tnosc, dtosc, cem, tkusc, kusc, dmkoosc, dkosct, dkoscm, dkoscl, por, ktpvoz, tevoz, Prvoz, tnac, dete, efftem, ektpm, 
		tgorm, tholm, qobm, kttk, stch, dmkvoz); 
	k=0; dpctp=mu[k]; k++; htcm=mu[k]; if (mu) delete[]mu;
	dkosp=1.0001190948;
	progpoTem(cem, ys0, snm, por, dh, tocras, qobm, efftem, kttk, stch, dpctp, ektpm, ktpvoz, tevoz, dmkvoz, vybves, vybmar, vyfv, 
	vysv, isrp, vpkf, vrsh, vystsha, kusc, tkusc, dmkoosc, dkoscm, dkosct, dkoscl, tholm, tgorm, dkoalmod, dkosp, kost, chPo, htcm);
	if (tkusc) delete[]tkusc; if (dkoscm) delete[]dkoscm; if (dkosct) delete[]dkosct; 
	if (kusc) delete[]kusc; if (vsm) delete[]vsm; if (snm) delete[]snm; 
	return mu;
}
double **zadNacZnaRazMas(double *efftem, int cem, double *dkoscm, double *dkosct, int dkoscl, double *kusc, double *tkusc, 
	int dmkoosc, double tnd, double dtd, int vybves, double tnac, double dete, double tnosc, double dtosc, int vybmar, int vyfv, char *snm)
{
	int k=6, nu=6, nw=6;
	double wmgo=0.0, wsio=0.0, walo=0.0, smgo=0.0, salo=0.0, ssio=0.0;
	double tm=0.0, ko=1e-2, hf=1e0, **mu=new double*[nu], *sov=new double[nw];
	if ((!mu) || (!sov)) { cout << snm << endl; k=getchar(); exit(1); }
	k=0; dkosct[k]=tnd; for (k=1; k<dkoscl; k++) dkosct[k]=dkosct[k-1]+dtd;
	//-----если выбран шамот
	if (!vybves) {
	k=0; dkoscm[k]=8.15;  k++; dkoscm[k]=6.87;  k++; dkoscm[k]=7.55; k++; 
		 dkoscm[k]=14.84; k++; dkoscm[k]=20.09; k++; dkoscm[k]=26.3; }
	//-----если выбран вермикулит
	if (vybves==1) {
	k=0; dkoscm[k]=4.68;  k++; dkoscm[k]=4.69; k++; dkoscm[k]=5.65; k++; 
		 dkoscm[k]=13.17; k++; dkoscm[k]=20.2; k++; dkoscm[k]=27.81;
		 ssio=36.8*ko; salo=13.2*ko; smgo=21.7*ko;
		 if (!vyfv) { ssio=45.8*ko; salo=14.0*ko; smgo=29.0*ko; } //2-0.7
		 else if (vyfv==1) { ssio=46.7*ko; salo=12.9*ko; smgo=28.2*ko; } //8-4
	wmgo=smgo; walo=salo; wsio=ssio; }
	//-----если выбран ИТОМ
	if (vybves==2) { 
		if (!vybmar) {
		salo=26e0*ko; smgo=22e0*ko; ssio=52e0*ko; //для трехкомпонентной смеси
		walo=23e0*ko; wmgo=19e0*ko; wsio=49e0*ko; //для многокомпонентной смеси
		k=0; dkoscm[k]=6.07;  k++; dkoscm[k]=5.36;  k++; dkoscm[k]=6.19; k++; 
			 dkoscm[k]=13.48; k++; dkoscm[k]=19.93; k++; dkoscm[k]=27.69;
	} //ИТОМ-440
	else if (vybmar==1) {
		salo=29e0*ko; smgo=16e0*ko; ssio=55e0*ko;
		walo=26e0*ko; wmgo=15e0*ko; wsio=5e1*ko; 
		k=0; dkoscm[k]=6.41;  k++; dkoscm[k]=5.53; k++;  dkoscm[k]=6.32; k++; 
			 dkoscm[k]=13.56; k++; dkoscm[k]=19.86; k++; dkoscm[k]=27.66;
	} //ИТОМ-620
	else if (vybmar==2) {
		salo=33e0*ko; smgo=11e0*ko; ssio=56e0*ko;
		walo=3e1*ko; wmgo=1e1*ko; wsio=52e0*ko; 
		k=0; dkoscm[k]=7.28;  k++; dkoscm[k]=6.1;   k++; dkoscm[k]=6.83; k++; 
			 dkoscm[k]=14.02; k++; dkoscm[k]=19.86; k++; dkoscm[k]=27.22;
	} //ИТОМ-860
	else if (vybmar==3) {
		salo=35e0*ko; smgo=9e0*ko; ssio=56e0*ko; walo=33e0*ko; wmgo=7e0*ko; wsio=53e0*ko; 
		k=0; dkoscm[k]=7.45; k++;  dkoscm[k]=6.24;  k++; dkoscm[k]=6.95; k++; 
			 dkoscm[k]=14.15; k++; dkoscm[k]=19.89; k++; dkoscm[k]=27.09;
	} //ИТОМ-1000
	else { cout << "Takoy marki ITOM net!" << endl; k = getchar(); exit(1); } }
	//-----если выбран КВИ
	if (vybves==3) { 
	k=0; dkoscm[k]=3e0;  k++; dkoscm[k]=6e0; k++; dkoscm[k]=7e0; k++; 
		 dkoscm[k]=15e0; k++; dkoscm[k]=2e1; k++; dkoscm[k]=26e0;
		if (vybmar>6) { k=1; dkoscm[k]=7e0; } if (vybmar>9) { k=2; dkoscm[k]=8e0; } }
	if ((vybves>=4) || (vybves<0)) { cout << "Takogo veschestva net!" << endl; k=getchar(); exit(1); }
	for (k=0; k<dkoscl; k++) dkoscm[k]=hf-dkoscm[k]*ko;
	//-----
	k=0; efftem[k]=tnac; for (k=1; k<cem; k++) efftem[k]=efftem[k-1]+dete;
	k=0; tkusc[k]=tnosc; for (k=1; k<dmkoosc; k++) tkusc[k]=tkusc[k-1]+dtosc;
	kusc=koefoslab(wmgo, wsio, walo, tkusc, dmkoosc, kusc, snm);
	k=0; sov[k]=wmgo; k++; sov[k]=wsio; k++; sov[k]=walo; k++; sov[k]=smgo; k++; sov[k]=salo; k++; sov[k]=ssio;
	k=0; mu[k]=efftem; k++; mu[k]=dkoscm; k++; mu[k]=dkosct; k++; mu[k]=kusc; k++; mu[k]=tkusc; k++; mu[k]=sov;
	return mu;
}
void vyvodfile(double *ao, int sao, int zf, double hvv, double p1, double p2, int vybves, int vybmar)
{
	int k=0, j=k+1, q=j+1, l=q+1, i=0, dlar=0; 
	char pbf='A', **mauk=napolStrok(vybves, vybmar, pbf);
	char *snm=mauk[k], *sfno=mauk[j]; 
	char **mu=napNazvFile(vybves, vybmar, snm, pbf);
	char *szfov=mu[k], *sppov=mu[j], *skpov=mu[q], *sdvov=mu[l]; 
	FILE *fo = fopen(szfov, "at");
	if (!fo) { cout << sfno << endl; k = getchar(); exit(1); }
	if (zf < 1) {
		fprintf(fo, "\nhv = %0.10lf\np1 = %lf\np2 = %lf\n", hvv, p1, p2); 
	} for (k = 0; k < sao; k++) fprintf(fo, "%0.15lf\n", ao[k]); fprintf(fo, "\n"); 
	for (k=0; k<=j; k++) { snm=mauk[k]; if (snm) delete[]snm; } if (mauk) delete[]mauk; 
	for (k=0; k<=l; k++) { snm=mu[k]; if (snm) delete[]snm; } if (mu) delete[]mu; 
	fclose(fo);
}
void progpoTem(int cem, double ys0, char *snm, double por, double dh, double tocras, double *qobm, double *efftem, double *ktptk, 
	double *stch, double *dpctp, double *ektpm, double *ktpvo, double *vte, int dmkvoz, int vybves, int vybmar, int vyfv, 
	int vysv, int isrp, int vpkf, int vrsh, int vystsha, double *kusc, double *tkusc, int dmkoosc, double *dkoscm, 
	double *dkosct, int dkoscl, double *tholm, double *tgorm, double dkoal, double dkosp, int kost, int chP, double *htcm)
{
	FILE *fo=NULL;
	double hk=0.0, h=0.0, hvna=0.0, *tx=NULL, *qx=NULL, *lx=NULL, tem=0.0, thol=0.0, tgor=0.0, ektp=0.0;
	double *kx=NULL, *dtx=NULL, *oox=NULL, *temras=NULL, hf=1e0, e=1e-10, p1=hf, p2=hf, pk=1e1, tsred=0.0;
	double qobsh=0.0, *poi=NULL, mrp=0.0, htch=0.0, dkusce=0.0, dkoscee=0.0, dpctpe=0.0, nzpa=2e0, nzpb=2e0;
	int k=0, vtn=0, vtk=cem, vt=0, knxt=0, nnxt=0, j=0, vp=1, vpp=1, dmko=0;
	int frko=0, frktr=0, vtko=1, vykoot=1, vfqi=0, kolPer=2; //задаем диапазон температур
	double **mu=NULL, hvot=0.0, hvo=0.0, p=0.0, r=0.0, ko=hf, *Rmp=NULL, *rts=NULL, ka=0.0, kb=0.0, *ktr=NULL;
	p=0.0; for (k=0; k<kost; k++) p=p+hf;
	if (ys0<hf) ko=1e3;
	j=nnxt; r=e; while (r<(ys0*ko)) { r=r+hf; j++; } knxt=j; dmko=knxt-nnxt;
	hvna=0.0; for (k=0; k<nnxt; k++) hvna=hvna+hf;
	rts=new double[dmko]; if (!rts) { cout << snm << endl; k=getchar(); exit(1); }
	/*vp=0; dkosp=RaschAlphaTvKar(vybves, vybmar, vystsha, vyfv, vysv, isrp, vpkf, por, snm, kost, dpctpe, vrsh, ktpvoz, tevoz, efftem, cem, 
		dmkvoz, kost, kttk, dmkoosc, kusc, tkusc, dkoscm, dkosct, dkoscl, thol, tgor, ys0, dkoal, dkosp, hk, qobsh, chPo, dkusce, 
		dkoscee, tc, stch, p1, p2, vp, Rmp, rts, dmko, vfqi); cout << "dkosp = " << dkosp << endl; */
	mu=initNovMas(dmko, snm);
k=0; tx=mu[k]; k++; qx=mu[k]; k++; lx=mu[k]; k++; kx=mu[k]; k++; dtx=mu[k]; k++; oox=mu[k]; k++; temras=mu[k]; 
if (mu) delete[]mu;
hvo=oprSrRazPor(vybves, vybmar, por, vyfv, vysv, isrp, vpkf, vrsh, vystsha, snm);
mrp=oprMaxRazPor(vybves, vybmar, por, vyfv, vysv, isrp, vpkf, vrsh, vystsha, snm);
htch=opredTolTverChasObFig(dpctpe, p, hvo, por, mrp);
			for (vt=vtn; vt<vtk; vt++) { //пробегаем по температуре 
				tgor=tgorm[vt]; thol=tholm[vt];
				ka=(thol-tgor)/ys0; kb=tgor; hk=hvna;
					for (j=nnxt; j<knxt; j++) { tem=ka*hk+kb; k=j-nnxt; rts[k]=tem; hk=hk+dh; } 
//-----//for (k=0; k<dmko; k++) cout << "T ( " << k << " ) = " << rts[k] << "\t"; cout << endl;
if (frko==1) Rmp=opredKoefOtrEdin(vybves, vybmar, snm, dkoal*dkosp, rts, dmko, vtko);
else if (!frko) Rmp=gotVelKoefOtr(vt, vybves, dmko, snm, vykoot);
if (frktr==1) ktr=KoefPoglRosselNac(rts, dmko, dkosp, dkoal, vybves, vybmar); 
else if (!frktr) ktr=gotVelKoefTeplRos(vt, vybves, dmko, snm); 
//for (k=0; k<dmko; k++) cout << "Ref ( " << k << " ) = " << Rmp[k] << "\ttem = " << rts[k] << "\tktr = " << ktr[k] << "\t"; cout << endl;
	/*fo=fopen("C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\2022\\Koeff_Tepl_Ros.txt", "at"); 
	fprintf(fo, "Коэффициент теплопроводности по Росселанду при T = %lf\n", efftem[vt]); 
	for (k=0; k<dmko; k++) fprintf(fo, "%0.15lf\n", ktr[k]); fprintf(fo, "\n"); 
	fclose(fo);
	fo=fopen("C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\2022\\Temp_Koeff_Tepl_Ros.txt", "at"); 
	fprintf(fo, "Температура T = %lf\n", efftem[vt]); 
	for (k=0; k<dmko; k++) fprintf(fo, "%0.15lf\n", rts[k]); fprintf(fo, "\n"); 
	fclose(fo);
	fo=fopen("C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\2022\\Koeff_Otrazh.txt", "at"); 
	fprintf(fo, "Температура T = %lf\n", efftem[vt]); 
	for (k=0; k<dmko; k++) fprintf(fo, "%0.15lf\n", Rmp[k]); fprintf(fo, "\n"); 
	fclose(fo);*/
	p1=nzpa*hf; p2=nzpb*hf;
					while (p1<pk) {
						hk=hvna; tem=efftem[vt]; 
						poi=provnaOtrits(kusc, tkusc, dmkoosc, dkoscm, dkosct, dkoscl, dpctp, efftem, cem, qobm, tholm, 
							tgorm, ektpm, ys0, snm, tem, htcm); 
						k=0; dkusce=poi[k]; k++; dkoscee=poi[k]; k++; dpctpe=poi[k]; k++; qobsh=poi[k]; 
						k++; thol=poi[k];	k++; tgor=poi[k];	 k++; ektp=poi[k];	 k++; htch=poi[k]; 
						if (htch<e) htch=poi[k]; 
						k++; tsred=poi[k];	if (poi) delete[]poi; 
						cout << "qobsh = " << qobsh << "\tth = " << thol << "\ttg = " << tgor << "\tktp = " << ektp << "\tht = " << htch << 
	"\tts = " << efftem[vt] << "\ttsr = " << tsred << endl;
					for (j=nnxt; j<knxt; j++) { //пробегаем по координате 
				k=j-nnxt; kx[k]=hk; 
	if (!vpp) { vp=1; qx[k]=BolTochRasAlpha(kost, hvo, por, tem, ktpvo, vte, efftem, cem, dmkvoz, snm, ktptk, vybves, vybmar, dmkoosc, kusc, tkusc, 
		dkoscm, dkosct, dkoscl, thol, tgor, ys0, dkoal, dkosp, hk, qobsh, chP, htch, dkusce, dkoscee, dpctpe, stch, p1, p2, vp, Rmp, rts, dmko, 
		vfqi); lx[k]=(tgor-thol)/ys0; if (fabs(lx[k])>e) lx[k]=qx[k]/lx[k]; else lx[k]=0.0; } //vp=1 - поиск ПТПИ
	else if (vpp==1) { 
		lx[k]=OprTempNewtRafsNach(kost, hvo, por, tem, ktpvo, vte, efftem, cem, dmkvoz, snm, ktptk, vybves, vybmar, dmkoosc, kusc,
		tkusc, dkoscm, dkosct, dkoscl, thol, tgor, ys0, dkoal, dkosp, hk, qobsh, chP, htch, dkusce, dkoscee, dpctpe, stch, p1, p2, Rmp, 
		rts, dmko, ektpm, ktr, dmko, dpctp, vfqi, dh, kolPer);
		qx[k]=lx[k]*(tgor-thol)/ys0; } //vpp==1 - инициация функции расчета ПТПИ методом Ньютона-Рафсона
				cout << "hk = " << hk << "\tte = " << tem  << "\tk = " << k << "\tqx = " << qx[k] << "\tp1 = " << p1 << "\tp2 = " << p2 << 
					"\tlx = " << lx[k] << endl;
				hk=hk+dh; }
					fo=fopen("C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\2022\\KTP_Tepl_Pot.txt", "at"); 
					fprintf(fo, "При T = %lf\tp1 = %0.10lf\tp2 = %0.10lf\thvo = %0.10lf\tOPTP = %0.10lf\n", efftem[vt], p1, p2, hvo, qobsh); 
					for (k=nnxt; k<knxt; k++) fprintf(fo, "%0.10lf\n", qx[k-nnxt]); fprintf(fo, "\n"); 
					fprintf(fo, "Координата\n"); 
					for (k=nnxt; k<knxt; k++) fprintf(fo, "%0.10lf\n", kx[k-nnxt]); fprintf(fo, "\n"); 
					fclose(fo); //k=knxt-nnxt; for (j=0; j<k; j++) cout << "j = " << j << "\tqx = " << qx[j] << "\thvo = " << hvo << "\t"; cout << endl; 
			p1=p1+hf; p2=p2+hf; }
	if (ktr) delete[]ktr; if (Rmp) delete[]Rmp; }
	if (tx) delete[]tx; if (qx) delete[]qx; if (lx) delete[]lx; if (kx) delete[]kx; 
	if (dtx) delete[]dtx; if (oox) delete[]oox; if (temras) delete[]temras; if (rts) delete[]rts; 
}
double *provnaOtrits(double *kusc, double *tkusc, int dmkoosc, double *dkoscm, double *dkosct, int dkoscl, double *dpctp, double *efftem, 
	int cem, double *qobm, double *tholm, double *tgorm, double *ektpm, double tolob, char *snm, double tc, double *htcm)
{
	int k=9; 
	double hf=1e0, *vm=new double[k], tsred=0.0;
	if (!vm) { cout << snm << endl; k=getchar(); exit(1); }
	double dkusce=opredKTPTKToch(kusc, tkusc, tc, dmkoosc); if (dkusce>hf) dkusce=hf; if (dkusce<0.0) dkusce=0.0;
	double dkoscee=opredKTPTKToch(dkoscm, dkosct, tc, dkoscl); if (dkoscee>hf) dkoscee=hf; if (dkoscee<0.0) dkoscee=0.0;
	double dpctpe=opredKTPTKToch(dpctp, efftem, tc, cem); if (dpctpe>hf) dpctpe=hf; if (dpctpe<0.0) dpctpe=0.0;
	double htc=opredKTPTKToch(htcm, efftem, tc, cem); if (dpctpe<0.0) dpctpe=0.0;
	double qobsh=hf, thol=hf, tgor=hf, ektp=hf, t=hf;
	qobsh=opredKTPTKToch(qobm, efftem, tc, cem); 
	thol=opredKTPTKToch(tholm, efftem, tc, cem); 
	tgor=opredKTPTKToch(tgorm, efftem, tc, cem); 
	ektp=opredKTPTKToch(ektpm, efftem, tc, cem); 
	tsred=opredKTPTKToch(efftem, efftem, tc, cem);
	if (ektp>0.0) {
		if (qobsh>0.0) {
			if (thol>0.0) {
				if (tgor<0.0) { t=qobsh*tolob/ektp; tgor=thol+t; }
			} else { t=qobsh*tolob/ektp; if (tgor>t) thol=tgor-t; }
		} else { if ((tgor>0.0) && (thol>0.0) && (tgor>thol)) qobsh=ektp*(tgor-thol)/tolob; else qobsh=0.0; }
	} else { if ((tgor>0.0) && (thol>0.0) && (tgor>thol) && (qobsh>0.0)) ektp=qobsh*tolob/(tgor-thol); else ektp=0.0; }
	if (tsred<0.0) tsred=(tgor+thol)/2e0;
	k=0; vm[k]=dkusce; k++; vm[k]=dkoscee; k++; vm[k]=dpctpe; k++; vm[k]=qobsh;  
	k++; vm[k]=thol;   k++; vm[k]=tgor;    k++; vm[k]=ektp;   k++; vm[k]=htc;
	k++; vm[k]=tsred;
	return vm;
}
double oprKoefPogMod(int vybves)
{
	double al=0.0, sal=0.0, hf=1e0, ko=0.0;
	int kon=0, vymo=0, mvyko=0, mvymo=0, vyfr=0, mvyfr=0, k=0;
	int nvymo=0, nvyfr=0, nvyko=0;
	mvymo=1; mvyko=2; mvyfr=1; if (vybves==1) { mvymo=3; mvyko=3; mvyfr=4; }
	if ((vybves>3) || (vybves<0)) { cout << "Oshibka v vybore veschestva!" << endl; k=getchar(); exit(1); }
	for (vymo=nvymo; vymo<mvymo; vymo++) {
		if ((vymo>0) && (vybves==1)) nvyfr=1;
		for (vyfr=nvyfr; vyfr<mvyfr; vyfr++) {
			for (kon=nvyko; kon<mvyko; kon++) { 
if (!vybves)
	if (!vymo) {
		if (!vyfr) {
		if (!kon) al=0.6259373;
		if (kon==1) al=0.678323; } }
if (vybves==1) { 
	if (!vymo) { //модель ППП
		if (!vyfr) {
			if (!kon) al=0.544286;
			if (kon==1) al=0.580206;
			if (kon==2) al=0.624646; }
		if (vyfr==1) {
			if (!kon) al=0.532591;
			if (kon==1) al=0.555503;
			if (kon==2) al=0.585634; }
		if (vyfr==2) {
			if (!kon) al=0.483897;
			if (kon==1) al=0.528810;
			if (kon==2) al=0.541856; }
		if (vyfr==3) {
			if (!kon) al=0.515290;
			if (kon==1) al=0.519519;
			if (kon==2) al=0.532284; } }
	if (vymo==1) { //модель шаров
		if (vyfr==1) {
			if (!kon) al=0.369115;
			if (kon==1) al=0.399028;
			if (kon==2) al=0.421853; }
		if (vyfr==2) {
			if (!kon) al=0.426855;
			if (kon==1) al=0.422853;
			if (kon==2) al=0.396290; }
		if (vyfr==3) {
			if (!kon) al=0.392310;
			if (kon==1) al=0.413473;
			if (kon==2) al=0.414668; } } 
	if (vymo==2) { //модель цилиндров
		if (vyfr==1) { 
			if (!kon) al=0.494147;
			if (ko==1) al=0.485141;
			if (kon==2) al=0.544403; }
		if (vyfr==2) { 
			if (!kon) al=0.604525;
			if (kon==1) al=0.425123;
			if (kon==2) al=0.493076; }
		if (vyfr==3) {
			if (!kon) al=0.341095;
			if (kon==1) al=0.423071;
			if (kon==2) al=0.363935; } } }
if (vybves==2)
	if (!vymo) {
		if (!vyfr) {
		if (!kon) al=0.628873;
		if (kon==1) al=0.649737; } }
	sal=sal+al; ko=ko+hf; } } }
	sal=sal/ko;
return sal;
}
double OprTempNewtRafsNach(int kost, double hvp, double por, double tc, double *ktpvo, double *vte, double *ete, int cem, int dmkvoz, 
	char *snm, double *ktptk, int vybves, int vybmar, int dmkoosc, double *kusc, double *tkusc, double *dkoscm, double *dkosct, 
	int dkoscl, double thol, double tgor, double y0, double dkoal, double dkosp, double hk, double qobsh, int kolPok, double htk, 
	double dkusce, double dkoscee, double dpctpe, double *stch, double p1, double p2, double *Refm, double *rts, int dmrts, 
	double *laefm, double *ktr, int dmktr, double *dpctp, int vfqi, double dx, int kolPer)
{
	int j=0, nk=j, l=0, k=0, zf=0, iz=1, v=0, c=0, w=0, kevm=5, vp=0, n=2*kost; //w - расчет температуры, iz - изменение
	double ***protv=NULL, ***prots=NULL, hf=1e0, hlr0=hf, hrl0=0.0, vyzn=0.0, ra=fabs(hlr0-hrl0), tol=0.0, *ptr=NULL;
	double *Ts=NULL, *Tss=NULL, **mu=NULL, **mauk=NULL, *po=NULL, *sislr=NULL, *sisrl=NULL, *qxi=NULL;
	double d=0.0, ka=d, kb=d, e=1e-9, mo=-hf, tx=tgor-fabs(tgor-thol)*hk/y0, tora=1e-4, tma=d, tmi=d, ts=0.0;
	double *hrl=NULL, *hlr=NULL, *te=new double[n], r=0.0, *tem=NULL, lambda=0.0, tsred=0.0, *ktpr=new double[n];
	if ((!te) || (!ktpr)) { cout << snm << endl; k=getchar(); exit(1); }
	for (j=0; j<kost; j++) d=d+hf; tol = d*htk + (d - hf)*hvp;
	if (vybves>=4) { cout << "Net takogo veschestva!" << endl; k = getchar(); exit(1); }
	mauk=RaschRTA(kost, htk, ka, kb, iz, hvp, v, tc, c, w, vybves, cem, ete, hk, vybmar, dkusce, dkoscee, thol, tgor, y0, snm, dkosp, 
		kusc, tkusc, dmkoosc, dkoscm, dkosct, dkoscl, dkoal, Refm, rts, dmrts, p1*hvp, p2*hvp);
	double *Refa=NULL, *Traa=NULL, *Aba=NULL, *Ref=NULL, *Tra=NULL, *Ab=NULL, *Reft=NULL, *Trat=NULL, *Abt=NULL;
	k=0; Refa=mauk[k]; k++; Traa=mauk[k]; k++; Aba=mauk[k];  k++; Ref=mauk[k];
	k++; Tra=mauk[k];  k++; Ab=mauk[k];   k++; Reft=mauk[k]; k++; Trat=mauk[k]; k++; Abt=mauk[k]; k++; nk=k; //for (k=0; k<kost; k++) cout << "R = " << Ref[k] << "\tT = " << Tra[k] << "\tA = " << Ab[k] << endl;
	k=0; protv=zadrkt(kost, snm, kolPok, Ref, Tra, htk, hvp, p1*hvp, p2*hvp, k); //vyv - выбор вещества: 0 - шамот, 1 - вермикулит
	k++; prots=zadrkt(kost, snm, kolPok, Reft, Trat, htk, hvp, p1*hvp, p2*hvp, k);
	mu=opredTempLPSten(tc, n, ktpvo, vte, ete, ktptk, cem, dmkvoz, htk, hvp, qobsh, kost, snm, mo);
	k=0; Ts=mu[k]; k++; Tss=mu[k]; if (mu) delete[]mu; 
	j=0; for (k=0; k<n; k++) { if (!(k%2)) te[k]=Ts[j]; else { te[k]=Tss[j]; j++; } }
	k=0; tmi=Ts[k]; tma=Ts[k];
	for (k=0; k<kost; k++) { if (tmi>Ts[k]) tmi=Ts[k]; if (tmi>Tss[k]) tmi=Tss[k]; 
	if (tma<Ts[k]) tmi=Ts[k]; if (tma<Tss[k]) tmi=Tss[k]; } //for (k=0; k<cem; k++) cout << "ktr ( " << k << " ) = " << ktr[k] << "\ttem = " << ete[k] << "\t"; cout << endl;
	hlr0=opredKTPTKToch(ktr, rts, tc, dmktr); if (hlr0<0.0) hlr0=0.0; //cout << "KTP_Ros = " << hlr0 << endl;
	hlr0=fabs(hlr0*opredTempStenFragm(tc, n, ktpvo, vte, ete, ktptk, cem, dmkvoz, htk, hvp, qobsh, mo)/tol); cout << "hlr0 = " << hlr0 << "\t";
	r=0.0; k=0; vp=2; 
	hlr=RasIzlSerStenNacSok(Ref, Tra, Ab, protv, prots, tc, hlr0, hrl0, por, dpctpe, htk, hvp, kost, ete, ktptk, cem, ktpvo, vte, snm, stch, 
		dmkvoz, qobsh, hk, p1, p2, tx, dkosp, dkoal, vybves, vybmar, vp, ktr, dmktr, r, k, -1, vfqi, kolPok); //for (k=0; k<kost; k++) cout << "Ts ( " << k << " ) = " << Ts[k] << "\tTss = " << Tss[k] << "\thlr = " << hlr[k] << "\t"; cout << endl; 
	r=0.0; k=0; vp=3;
	hrl=RasIzlSerStenNacSok(Ref, Tra, Ab, protv, prots, tc, hlr0, hrl0, por, dpctpe, htk, hvp, kost, ete, ktptk, cem, ktpvo, vte, snm, stch, 
		dmkvoz, qobsh, hk, p1, p2, tx, dkosp, dkoal, vybves, vybmar, vp, ktr, dmktr, r, k, -1, vfqi, kolPok); //for (k=0; k<kost; k++) cout << "hlr ( " << k << " ) = " << hlr[k] << "\thrl = " << hrl[k] << "\t"; cout << endl; 
	tem=reshnewtrafs(protv, prots, hlr, hrl, Tra, Ab, kost, ktptk, ete, cem, ktpvo, vte, dmkvoz, htk, hvp, te, tora, laefm, kevm, tma, tmi, 
		y0, snm, qobsh, Ref, hlr0, hrl0, vybves, por, stch, hk, p1, p2, dkosp, dkoal, vybmar, vp, ktr, dmktr, dpctp, vfqi, kolPok, kolPer); //for (k=0; k<(2*kost); k++) cout << "tem ( " << k << " ) = " << te[k] << "\t";
	if (tem) { k=0; tmi=tem[k]; tma=tem[k];
	for (k=0; k<kost; k++) { if (tmi>tem[k]) tmi=tem[k]; 
	if (tma<tem[k]) tma=tem[k]; } tsred=(tmi+tma)/2e0; 
	lambda=fabs(tma-tmi)/tol; 
	if (fabs(lambda)>e) lambda=qobsh/lambda; else lambda=0.0;  
	cout << "ts = " << tsred << "\tlam = " << lambda << "\ttmi = " << tmi << "\ttma = " << tma << endl; 
	for (k=0; k<n; k++) ktpr[k]=0.0; //for (k=0; k<n; k++) cout << "tem ( " << k << " ) = " << tem[k] << "\t"; cout << endl;
	ts=0.0; for (k=0; k<n; k++) ts=ts+tem[k];
	if (fabs(lambda)>e) {
	for (k=1; k<n; k++) { lambda=(tem[k]-tem[k-1])/dx; 
	if (fabs(lambda)>e) ktpr[k]=qobsh/lambda; else ktpr[k]=0.0; }
	lambda=0.0; d=0.0; for (k=0; k<n; k++) if (ktpr[k]>e) { lambda=lambda+ktpr[k]; d=d+hf; }
	lambda=lambda/d; } else lambda=0.0; }
	v=4; for (w=0; w<v; w++) { for (k=0; k<kost; k++) { 
		ptr=prots[w][k]; if (ptr) delete[]ptr; ptr=protv[w][k]; if (ptr) delete[]ptr; } }
	if (protv) delete[]protv; if (prots) delete[]prots; if (Ts) delete[]Ts; if (Tss) delete[]Tss; 
	for (k=0; k<nk; k++) { po=mauk[k]; if (po) delete[]po; } if (mauk) delete[]mauk; 
	if (te) delete[]te; if (tem) delete[]tem; if (ktpr) delete[]ktpr; //k=getchar();
	return lambda;
}