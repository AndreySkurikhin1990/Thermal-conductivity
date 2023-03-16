#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <cmath>
#include <windows.h>
#include <iostream>
using namespace std;
const double pi = acos(-1e0), porex=95.0*1e-2, tocras=1e-8;
const double porist=28.0*1e-2, tem0=273.15, pgl=1e5, lgl=1e-9;
const int dsnftk=50, hko=1000, dmvtk=28, dmsrps=23, cemtk=11;
struct raspporporazmm { double rprm; struct raspporporazmm *slel; };
double *ktpvotk=NULL, *vtetk=NULL, *Prvo=NULL, *stch=NULL;
char *sntk=NULL, *sstk=NULL, *szfvtk=NULL, *szfstk=NULL, *szfstkk=NULL;
int dmsrpv=18, dmsrpi=17, dmrprs3=31, dmrprs4=32, dmrprs5=25, dmrprs6=33, dmkvi=1;
double ReshNeyaUravMaxKoTepl(double, double, double); //Максимальный
double ReshNeyaUravMinKoTepl(double, double, double); //Минимальный
double ReshNeyaUravKoTeplLiechtenecker(double, double, double); //Лихтенэкер
double opredRdeltaSrVzaPro(double, double, double, int); //Взаимопроникающие компоненты
double KoTeplFrickEllipsoidChasti(double, double, double); //Фрик - эллипсоиды
double ReshNeyaUravKoTeplMaxwell(double, double, double); //Максвелл
double ReshNeyaUravKoTeplBurgFrik(double, double , double); //Бургер-Фрик
double ReshNeyaUravKoTeplShumanVoss(double, double, double); //Шуман-Восс
double ReshNeyaUravKoTeplBruggman(double, double, double); //Бруггман
double ReshNeyaUravKoTeplGoringChirchillSpher(double, double, double); //Горринг-Черчилль (сферы)
double ReshNeyaUravKoTeplGoringChirchillCyl(double, double, double); //Горринг-Черчилль (цилиндры)
double ReshNeyaUravKoTeplGoringChirchillCylPlotSl(double, double, double, double); //Горринг-Черчилль (плотный слой)
double ReshNeyaUravKoTeplMaxwAcken(double, double, double); //Максвелл-Эйкен
double ReshNeyaUravKoTeplSwift(double, double); //Свифт
double MetodHashinaShtrikmanaMin(double, double, double); //Хашин-Штрикман (минимальный)
double MetodHashinaShtrikmanaMax(double, double, double); //Хашин-Штрикман (максимальный)
double RdeltaMin(double, double, double, double,double, double);
double RdeltaMax(double, double, double, double, double, double, double);
double urovOtsechen(double, double, double);
double PoiskPShumanVoss(double);
double PoiskDelta(double, double);
double opredFiBolN(double, double);
double utochLamTverKark(double, double, double, double, double, double, double, double, double);
double raschTeplVzaimoPronKomp(double, double, double, double);
double *DulnevZern(double, double *, double *, int, double *, double *, double, double, double, double, int, double, double, int, double *, double *, double *);
double *opredKTPTverKarkSha(double *, double *, double, double, double, double, double, double, int, double *, double *, double *, int, int, int);
double *opredKTPTverKarkShamota(double *, double *, double, double, double, double, int, double, double, int, double *, double *, double *, double *, double *, int, double);
double opredTeploprovVozd(double);
double DulnevSloForSha(double, double, double, double, double, double);
double MoleSostTeplVozd(double, double, double);
double opredFiBolNN(double, double);
double *opredTvChaKoeTepSrInSpSha(int, double, double *, double *, int, double *, double *, double);
double DulKoefTep1Adi(double, double, double);
double DulKoefTep1Izoterm(double, double, double);
double DulKoefTep1OdolevMatr(double, double, double);
double DulKoefTep1OdolevStatSm(double, double, double);
double DulKoefTep1AdiCyl(double, double, double, double);
double DulKoefTep1IzotermCyl(double, double, double, double);
double *BezUchetaRaspedPorSha(double, double *, double *, double *, int, int, double *, double);
double *UchetRapsredPorPoRazmSha(double, double *, double *, double *, int, int, double *, int, double *, double *);
double *UchetRapsredPorPoRazmVer(double, double *, double *, double *, int, int, double *, int, double *, double *, int, int, int);
double *DulnevKoefTepShamot(int, double, double *, int, double *, double *, double *, double, int, double *, double *);
double *DulnevKoefTepVerm(int, double, double *, int, double *, double *, double *, double, int, double *, double *, int, int, int, int);
void zapisvfile(double *, int, char *);
double *oprKTPTverKarkVermi(double *, double *, double, double, double, double, int, int, double, double, int, double *, double *, double *, int, int, double *, double *, int, double, int);
double urovPodder(double, double, double);
double sign(double);
double DulKoefTep1CylPer(double, double, double, double);
void initpervznach(int, double *, int);
double opredKoefAkkomodLandau(double);
double opredPlotnVozd(double);
double *arrTempAir();
double *arrPlotnAir();
double *rasPorpoRazSha(double, int, int);
double opredDulnLam1(double, double, double, double, double, double);
double *MetodDulnevaSigalovoyPolyDispVer(double *, int, double, double *, int, double);
double *MetodDulnevaSigalovoyDopVer(double *, int, double, double *, int, double);
double *MetodDulnevaSigalovoyBezIspVer(double *, int, double, double *, int, double);
double DulnevSloForVer(double, double, double, double, double, double);
void NapMasVozdSha(double *, double *, int);
void NapMasPr();
double *opredKTPTverKarkVerm(double *, double *, double, double, double, double, int, int, double, double, int, double *, double *, double *, int, int, int, int, int);
double opredUrovPodderM03(double);
double *VasFraySha(double *, int);
double SysPlastShakhPor(double, double, double);
double MetodStarostina(double, double, double);
double MetodRusselNeprVozd(double, double, double);
double MetodRusselNeprTverTelo(double, double, double);
double *VasFrayObsh(double *, double *, double, int);
double Poiskkx(double);
void osvpamtk();
double *opredKTPTverKarkitom(double *, double *, double, double, double, double, double, double, int, double *, double *, double *, int, int, int);
double *oprKTPTverKarkitom(double *, double *, double, double, double, double, int, double, double, int, double *, double *, double *, int, double *, double *, int, int, int, double);
double *DopFunctOpredKTPTverChas(double *, double *, int, double, double *, int, int, double *, double);
double **vybFunRasPorpoRazSha(double, int, int);
double *oprKTPTverKarkkvi(double *, double *, double, double, double, double, int, double, double, int, double *, double *, double *, int, double *, double *, int, int, int, double);
void napStrTveKar(int);
double opredKTPTKTochSha(double *, double *, double, int);
double urovPod(double);
double *proverkakvi(double *, double *, double *, double, double, int);
double **rasPorpoRazkvi(int);
double **rasPorpoRazitom(int);
double **rasPorpoRazVer(double, int, int, int, int, int);
double **opredSredKTPTK();
double **initarrver(int, int, int, int, int);
//------------------------
void napStrTveKar(int vy)
{ int k=0; char l1, l2, l3, l4='_';
sntk[k]='D'; k++; sntk[k]=':'; k++; sntk[k]='\\'; k++; sntk[k]='\\'; k++; sntk[k]='_';  k++; sntk[k]='А';  k++;
sntk[k]='с'; k++; sntk[k]='п'; k++; sntk[k]='и';  k++; sntk[k]='р';  k++; sntk[k]='а';  k++; sntk[k]='н';  k++;
sntk[k]='т'; k++; sntk[k]='у'; k++; sntk[k]='р';  k++; sntk[k]='а';  k++; sntk[k]='\\'; k++; sntk[k]='\\'; k++;
sntk[k]='t'; k++; sntk[k]='m'; k++; sntk[k]='p';  k++; sntk[k]='\\'; k++; sntk[k]='\\'; k++; sntk[k]='\0';
k=0;
sstk[k]='C';  k++; sstk[k]=':';  k++; sstk[k]='\\'; k++; sstk[k]='\\'; k++; sstk[k]='U';  k++; sstk[k]='s';  k++;
sstk[k]='e';  k++; sstk[k]='r';  k++; sstk[k]='s';  k++; sstk[k]='\\'; k++; sstk[k]='\\'; k++; sstk[k]='А';  k++;
sstk[k]='н';  k++; sstk[k]='д';  k++; sstk[k]='р';  k++; sstk[k]='е';  k++; sstk[k]='й';  k++; sstk[k]='\\'; k++;
sstk[k]='\\'; k++; sstk[k]='D';  k++; sstk[k]='o';  k++; sstk[k]='c';  k++; sstk[k]='u';  k++; sstk[k]='m';  k++;
sstk[k]='e';  k++; sstk[k]='n';  k++; sstk[k]='t';  k++; sstk[k]='s';  k++; sstk[k]='\\'; k++; sstk[k]='\\'; k++;
sstk[k]='_';  k++; sstk[k]='А';  k++; sstk[k]='с';  k++; sstk[k]='п';  k++; sstk[k]='и';  k++; sstk[k]='р';  k++;
sstk[k]='а';  k++; sstk[k]='н';  k++; sstk[k]='т';  k++; sstk[k]='у';  k++; sstk[k]='р';  k++; sstk[k]='а';  k++;
sstk[k]='\\'; k++; sstk[k]='\\'; k++; sstk[k]='t';  k++; sstk[k]='m';  k++; sstk[k]='p';  k++; sstk[k]='\\'; k++;
sstk[k]='\\'; k++; sstk[k]='\0';
if (vy==0) { l1='V'; l2='e'; l3='r'; } 
else if (vy==1) { l1='S'; l2='h'; l3='a'; } 
else if (vy==2) { l1='I'; l2='t'; l3='o'; l4='m'; }
else if (vy==3) { l1='K'; l2='V'; l3='I'; }
k=0;
szfstkk[k]='O'; k++; szfstkk[k]='u'; k++; szfstkk[k]='t'; k++; szfstkk[k]='p'; k++; szfstkk[k]='u'; k++; szfstkk[k]='t'; k++;
szfstkk[k]='_'; k++; szfstkk[k]='T'; k++; szfstkk[k]='K'; k++; szfstkk[k]='_'; k++; szfstkk[k]=l1; k++; szfstkk[k]=l2; k++;
szfstkk[k]=l3; k++; szfstkk[k]=l4; k++; szfstkk[k]='.'; k++; szfstkk[k]='t'; k++; szfstkk[k]='x'; k++; szfstkk[k]='t'; k++; szfstkk[k]='\0';
strcpy(szfstk,sstk); strcat(szfstk,szfstkk); szfstk[strlen(szfstk)+1]='\0'; }
void osvpamtk()
{ delete []stch; delete []ktpvotk; delete []vtetk; delete []Prvo; delete []sntk; 
delete []sstk; delete []szfvtk; delete []szfstk; delete []szfstkk; }
void initpervznach(int lentem, double *stchn, int vyb)
{ int k, u=lentem; ktpvotk=new double[dmvtk]; vtetk=new double[dmvtk]; 
if ((ktpvotk) && (vtetk)) NapMasVozdSha(ktpvotk, vtetk, dmvtk); else { cout << "No memory!" << endl; k=getchar(); exit(1); }
Prvo=new double[dmvtk]; if (!Prvo) { cout << "No memory!" << endl; k=getchar(); exit(1); } else NapMasPr();
sntk=new char[dsnftk]; sstk=new char[dsnftk];
szfvtk=new char[2*dsnftk]; szfstk=new char[2*dsnftk]; szfstkk=new char[dsnftk];
if ((!sntk) || (!sstk) || (!szfvtk) || (!szfstk) || (!szfstkk)) { cout << "No memory!" << endl; k=getchar(); exit(1); } 
for (k=0; k<dsnftk; k++) { sntk[k]='\0'; sstk[k]='\0'; szfstkk[k]='\0'; } 
for (k=0; k<(2*dsnftk); k++) { szfvtk[k]='\0'; szfstk[k]='\0'; } 
napStrTveKar(vyb); 
stch=new double[u]; if (!stch) { cout << "No memory!" << endl; k=getchar(); exit(1); } for (k=0; k<u; k++) stch[k]=stchn[k]; }
//Управляющие функции
double *opredKTPTverKarkSha(double *tem, double *laefm, double por, double wsio, double walo, double wmgo, double tena, double dtos, int dmkusc, double *kuscs, double *tkuscs, double *stchs, int lentem, int vrsh, int vystsha)
{ int u=lentem, k=1, dmi=0;
initpervznach(u, stchs, k); 
double **mu=vybFunRasPorpoRazSha(por, vrsh, vystsha), *raspr=mu[0];
double *srra=mu[1], *legr=mu[2], *prgr=mu[3], *ms=mu[4], srp=ms[0];
double *mm=mu[5], marp=mm[0], t=tocras, ht=1e0; 
dmi=0; while (t<marp) { dmi++; t=t+ht; } cout << "dmi = " << dmi << endl;
double *ktptks=opredKTPTverKarkShamota(tem, laefm, por, wsio, walo, wmgo, u, tena, dtos, dmkusc, kuscs, tkuscs, stchs, srra, raspr, dmi, srp);
osvpamtk(); if (srra) delete[]srra; if (raspr) delete[]raspr; if (legr) delete[]legr; if (prgr) delete[]prgr;
if (mm) delete[]mm; if (ms) delete[]ms; return ktptks; }
double *opredKTPTverKarkShamota(double *tem, double *laefm, double por, double wsio, double walo, double wmgo, int u, double tena, double dtos, int dmkusc, double *kuscs, double *tkuscs, double *stchs, double *srra, double *rpr, int dms, double srp)
{ int k, cvp=22, cve=cvp+7*2, cved=cve+1, f, na=0, j, q=0, *no=new int[cved];
double ht=1e0, *srk=new double[u], **sr=new double*[cved], *lavo=new double[u], ep=tocras, *srzn=new double[u], y, pormax=5e-1, lvmi=ht;
if ((!no) || (!sr) || (!lavo) || (!srzn) || (!srk)) { cout << "No memory!"; k=getchar(); exit(1);} for (k=0; k<cved; k++) no[k]=0;
double s=0.0, p=0.0, x0=fabs(srp)/pow(por,1e0/3e0)/2e0; for (k=0; k<u; k++) lavo[k]=opredTeploprovVozd(tem[k]);
for (k=na; k<cved; k++) { srk=new double[u]; if (!srk) { cout << "No memory!"; j=getchar(); exit(1);} for (j=0; j<u; j++) srk[j]=0.0;
    if (k==na) srk=DulnevZern(por, tem, srk, na, laefm, lavo, wsio, walo, wmgo, x0, u, tena, dtos, dmkusc, kuscs, tkuscs, stchs);
	else if ((k<cvp) && (k>na)) srk=opredTvChaKoeTepSrInSpSha(k, por, tem, laefm, u, srk, lavo, x0);
    else if ((k<cve) && (k>=cvp)) srk=DulnevKoefTepShamot(k-cvp, por, tem, u, srk, laefm, lavo, srp, (cve-cvp)/2, srra, rpr);
	else if ((k<cved) && (k>=cve)) srk=VasFraySha(tem, u);
    sr[k]=srk; }
q=0; p=0.0; k=0; for (j=0; j<cved; j++) if (fabs(sr[j][k])>lvmi) { p=p+ht; no[q]=j; q++; } 
for (k=0; k<u; k++) { s=0.0; for (j=0; j<q; j++) s=s+sr[no[j]][k]; 
srzn[k]=s/p; cout << "kttk = " << srzn[k] << "\tlaef = " << laefm[k] << "\t"; 
} cout << endl; //j=0; for (k=0; k<cved; k++) { if (k==no[j]) { srk=sr[k]; zapisvfile(srk, u, szfstk); j++; } } 
cout  << "Nomera nenulevyh elementov" << endl; for (k=0; k<q; k++) cout << no[k]+1 << " "; cout << endl << "q = " << q << endl; 
for (k=0; k<cved; k++) { srk=sr[k]; if (srk) delete []srk; }
if (sr) delete []sr; if (lavo) delete []lavo; if (no) delete []no; 
return srzn; }
double *opredKTPTverKarkVerm(double *tem, double *laefm, double por, double wsio, double walo, double wmgo, int vfv, int n, double tena, double dtos, int dmkusc, double *kuscv, double *tkuscv, double *stchv, int lentem, int vysove, int vpkf, int isrp, int fl)
{ int m=lentem, k=0, dmi=0; initpervznach(m, stchv, k); 
double **mu=rasPorpoRazVer(por, vfv, k, vysove, isrp, vpkf); //0 - старые, 1 - новые значения
double *raspr=NULL, *srra=NULL, e=1e-6, koef=1e6;
double *prgr=NULL, *legr=NULL, *ms=NULL, srp=0.0, *mm=NULL, marp=0.0, t=tocras, ht=1e0;
k=0; raspr=mu[k]; k++; srra=mu[k]; k++; prgr=mu[k]; 
k++; legr=mu[k];  k++; ms=mu[k];   k++; mm=mu[k]; 
k=0; srp=ms[k]; marp=mm[k];
t=e; dmi=0; while (t<(marp*koef)) { dmi++; t=t+ht; } //cout << "dmi = " << dmi << endl;
double *ktptsv=oprKTPTverKarkVermi(tem, laefm, por, wsio, walo, wmgo, vfv, m, tena, dtos, dmkusc, kuscv, tkuscv, stchv, vysove, vpkf, srra, raspr, dmi, srp, fl);
osvpamtk(); if (raspr) delete[]raspr; if (srra) delete[]srra; if (legr) delete[]legr; 
if (prgr) delete[]prgr; if (mm) delete[]mm; if (ms) delete[]ms;
if (mu) delete[]mu; return ktptsv; }
double *oprKTPTverKarkVermi(double *tem, double *laefm, double por, double wsio, double walo, double wmgo, int vfv, int n, double tena, double dtos, int dmkusc, double *kuscv, double *tkuscv, double *stchv, int vysove, int vpkf, double *srra, double *raspr, int dmv, double srp, int fl)
{ 
int k, cvp=22, cve=cvp+7*2, cved=cve+4, na=0, j=0, q=0, *no=new int[cved], m=dmsrpv, f=(cve-cvp)/2, vv=0, u=n, w=0, pf=0;
double *srk=NULL, **sr=new double*[cved], *lavo=new double[n], ep=tocras, rmi=0.0, rma=0.0;
double *srzn=new double[n], s=0.0, p=0.0, srch=0.0, r, ht=1e0, lvmi=ht, ko=1e-3;
if ((!no) || (!sr) || (!lavo) || (!srzn)) { cout << "No memory!"; j=getchar(); exit(1); } 
for (k=0; k<cved; k++) { no[k]=0; sr[k]=NULL; }
for (k=0; k<n; k++) { lavo[k]=opredTeploprovVozd(tem[k]); //cout << "te = " << tem[k] << "\tlavo = " << lavo[k] << "\t";
if (lvmi>lavo[k]) lvmi=lavo[k]; }
if (!vfv) { rma=2e0; rmi=7e-1; }
else if (vfv==1) { rma=8e0; rmi=4e0; }
else if (vfv==2) { rma=1.6; rmi=0.35; } srch=(rmi+rma)*ko/2e0; 
w=0; p=0.0; for (k=na; k<cved; k++) {
srk=new double[n]; if (!srk) { cout << "No memory!"; j=getchar(); exit(1); } for (j=0; j<n; j++) srk[j]=0.0;
	if (k==na) srk=DulnevZern(por, tem, srk, na, laefm, lavo, wsio, walo, wmgo, srch, n, tena, dtos, dmkusc, kuscv, tkuscv, stchv);
	else if ((k<cvp) && (k>na)) srk=opredTvChaKoeTepSrInSpSha(k, por, tem, laefm, n, srk, lavo, srch);
    else if ((k<cve) && (k>=cvp)) srk=DulnevKoefTepVerm(k, por, tem, n, srk, laefm, lavo, srp, f, srra, raspr, 2, vv, dmv, cvp);
	else if ((k<cved) && (k>=cve)) srk=DopFunctOpredKTPTverChas(srk, tem, n, por, laefm, k-cve, vfv, lavo, srch); //if (((k+1)==30) || ((k+1)==35)) for (j=0; j<n; j++) cout << "j = " << j << "\tk = " << k << "\tsrk = " << srk[j] << endl; cout << endl;
	q=1; for (j=0; j<n; j++) if ((srk[j]<lvmi) && (q>0)) { q=-1; break; } 
	pf=1; for (j=1; j<n; j++) if ((srk[j]<srk[j-1]) && (pf>0)) { pf=-1; break; }
	if ((q>0) && (pf>0)) { sr[w]=srk; no[w]=k; w++; p=p+ht; } else delete[]srk; } 
	cved=w; for (j=0; j<n; j++) { 
		s=0.0; for (k=0; k<cved; k++) { w=no[k]; r=sr[k][j]; //cout << "\tr = " << r << "\tw = " << w+1 << endl; 
		s=s+r; } if (s>lvmi) { srzn[j]=s/p; //cout << "s = " << s/p << "\tj = " << j << endl; 
		} } 
	if (fl<0) { cout << "cved = " << cved << "\tl = " << lvmi << "\tsrp = " << srp << endl; 
	for (k=0; k<cved; k++) { w=no[k]; srk=sr[k]; if (srk) { for (j=0; j<n; j++) cout << "k = " << w+1 << "\ttem = " << tem[j] << "\tlam_tk = " << srk[j] << endl; cout << endl; } } 
	cout  << "Nomera nenulevyh elementov" << endl; for (k=0; k<cved; k++) cout << no[k]+1 << " "; cout << endl; 
	for (k=0; k<n; k++) cout << "T ( " << k << " ) = " << tem[k] << "\tlam_tk = " << srzn[k] << "\tlam_e = " << laefm[k] << endl; }
for (k=0; k<cved; k++) { srk=sr[k]; if (srk) delete []srk; }
if (sr) delete []sr; if (lavo) delete []lavo; if (no) delete []no; //k=getchar();
return srzn; }
double *opredKTPTverKarkitom(double *tem, double *laefm, double por, double wsio, double walo, double wmgo, int vvi, int n, double tena, double dtos, int dmkusc, double *kusci, double *tkusci, double *stchi, int lentem, int vybves, int vpv)
{ int u=lentem, k=3, dmi=0; 
initpervznach(u, stchi, 2); 
double **mu=rasPorpoRazitom(vvi), *raspr=mu[0], *srra=mu[1], *legr=mu[2];
double *prgr=mu[3], *ms=mu[4], *mm=mu[5], srp=ms[0], marp=mm[0];
double *ktptki=NULL, e=1e-1, t=e, ht=1e0; 
dmi=0; t=e; while (t<marp) { t=t+ht; dmi++; } //cout << "dmi = " << dmi << endl;
if (vybves) { if ((vvi>3) || (vvi<0)) { cout << "Marka ITOM ne naydena!" << endl; u=getchar(); exit(1); } } //for (k=0; k<m; k++) cout << "rpr ( " << k << " ) = " << raspr[k] << "\t srp ( " << k << ") = " << srra[k] << endl; 
ktptki=oprKTPTverKarkitom(tem, laefm, por, wsio, walo, wmgo, u, tena, dtos, dmkusc, kusci, tkusci, stchi, vvi, srra, raspr, vybves, vpv, dmi, srp);
osvpamtk(); if (raspr) delete[]raspr; if (srra) delete[]srra; if (legr) delete[]legr; 
if (prgr) delete[]prgr; if (ms) delete[]ms; if (mm) delete[]mm; if (mu) delete[]mu; //cout << "dm = " << dm << endl; 
return ktptki; }
double *oprKTPTverKarkitom(double *tem, double *laefm, double por, double wsio, double walo, double wmgo, int u, double tena, double dtos, int dmkusc, double *kusci, double *tkusci, double *stchi, int vvi, double *srra, double *raspr, int vybves, int vpv, int dmi, double srp)
{ int k, cvp=1+17+4, cve=cvp+12+2, cved=cve+1, f=(cve-cvp)/2, na=0, j, q=0, *no=new int[cved], n=u;
double *srk=NULL, **sr=new double*[cved], *lavo=new double[u], ep=tocras, *srzn=new double[u], y, pormax=83e-2, ht=1e0, lvmi=ht;
if ((!no) || (!sr) || (!lavo) || (!srzn)) { cout << "No memory!"; k=getchar(); exit(1);} for (k=0; k<cved; k++) no[k]=0;
double s=0.0, p=0.0, x0=fabs(srp)/pow(por,1e0/3e0)/2e0;
for (k=0; k<n; k++) { lavo[k]=opredTeploprovVozd(tem[k]); if (lvmi>lavo[k]) lvmi=lavo[k]; }
for (k=na; k<cved; k++) {
srk=new double[u]; if (!srk) { cout << "No memory!"; j=getchar(); exit(1);}
for (j=na; j<u; j++) srk[j]=0.0;
	if (k==na) srk=DulnevZern(por, tem, srk, na, laefm, lavo, wsio, walo, wmgo, x0, u, tena, dtos, dmkusc, kusci, tkusci, stchi);
    else if ((k<cvp) && (k>na)) { srk=opredTvChaKoeTepSrInSpSha(k, por, tem, laefm, u, srk, lavo, x0); }
    else if ((k<cve) && (k>=cvp)) { srk=DulnevKoefTepVerm(k, por, tem, u, srk, laefm, lavo, srp, f, srra, raspr, vybves, vpv, dmi, cvp); }
	else if ((k<cved) && (k>=cve)) { srk=VasFrayObsh(laefm, lavo, por, u); }
    sr[k]=srk; }
q=0; p=0.0; k=0; for (j=0; j<cved; j++) if (fabs(sr[j][k])>lvmi) { p=p+ht; no[q]=j; q++; }
for (k=0; k<u; k++) {
	s=0.0;
	for (j=0; j<q; j++) if (fabs(sr[no[j]][k])>lvmi) s=s+sr[no[j]][k]; srzn[k]=s/p; 
	cout << "kttk = " << srzn[k] << "\tlaef = " << laefm[k] << "\t"; } cout << endl; //j=0; for (k=0; k<cved; k++) { if (k==no[j]) { srk=sr[k]; zapisvfile(srk, u, szfstk); j++; } } 
cout  << "Nomera nenulevyh elementov" << endl; for (k=0; k<q; k++) cout << no[k]+1 << " "; cout << endl << "q = " << q << endl; 
for (k=0; k<cved; k++) { srk=sr[k]; if (srk) delete []srk; } 
if (sr) delete []sr; if (lavo) delete []lavo; if (no) delete []no; //k=getchar();
return srzn; }
double *opredKTPTverKarkkvi(double *tem, double *laefm, double por, double wsio, double walo, double wmgo, int vvi, int n, double tena, double dtos, int dmkusc, double *kusci, double *tkusci, double *stchkvi, int lentem, int vybves, int vpv)
{ int u=lentem, k=3, w=k, dmi=0;
initpervznach(u, stchkvi, w);
double **mu=rasPorpoRazkvi(vvi), *raspr=mu[0], *srra=mu[1], *legr=mu[2], *prgr=mu[3], *ms=mu[4], *mm=mu[5];
double srp=ms[0], marp=mm[0], *ktptkkvi=NULL, e=1e-1, t=e, ht=1e0; 
dmi=0; t=e; while (t<marp) { t=t+ht; dmi++; } //cout << "dm = " << dm << endl;
if (!vybves) { if ((vpv>10) || (vpv<4)) { cout << "Net takoy marki keramovermikulita!" << endl; k=getchar(); exit(1); } }
ktptkkvi=oprKTPTverKarkkvi(tem, laefm, por, wsio, walo, wmgo, u, tena, dtos, dmkusc, kusci, tkusci, stchkvi, vvi, srra, raspr, vybves, vpv, dmi, srp);
osvpamtk(); if (raspr) delete[]raspr; if (srra) delete[]srra; if (legr) delete[]legr; 
if (prgr) delete[]prgr; if (ms) delete[]ms; if (mm) delete[]mm;
return ktptkkvi; }
double *oprKTPTverKarkkvi(double *tem, double *laefm, double por, double wsio, double walo, double wmgo, int u, double tena, double dtos, int dmkusc, double *kusci, double *tkusci, double *stchk, int vvi, double *srra, double *raspr, int vybves, int vpv, int dmi, double srp)
{ int k, na=0, cvp=22, cve=cvp+7*2, cved=cve+1, f=(cve-cvp)/2, j, q=0, *no=new int[cved], w; 
double *srk, **sr=new double*[cved], *lavo=new double[u], ep=tocras, *srzn=new double[u], y, ht=1e0, lvmi=ht;
if ((!no) || (!sr) || (!lavo) || (!srzn)) { cout << "No memory!"; k=getchar(); exit(1);} for (k=0; k<cved; k++) no[k]=0;
double s=0.0, p=0.0, x0=fabs(srp)/pow(por,1e0/3e0)/2e0;
for (k=0; k<u; k++) { lavo[k]=opredTeploprovVozd(tem[k]); if (lvmi>lavo[k]) lvmi=lavo[k]; }
for (k=0; k<cved; k++) {
srk=new double[u]; if (!srk) { cout << "No memory!"; j=getchar(); exit(1);}
for (j=0; j<u; j++) srk[j]=0.0;
	if (!k) srk=DulnevZern(por, tem, srk, na, laefm, lavo, wsio, walo, wmgo, x0, u, tena, dtos, dmkusc, kusci, tkusci, stchk);
    else if ((k<cvp) && (k>0)) { srk=opredTvChaKoeTepSrInSpSha(k, por, tem, laefm, u, srk, lavo, x0); }
    else if ((k<cve) && (k>=cvp)) { srk=DulnevKoefTepVerm(k, por, tem, u, srk, laefm, lavo, srp, f, srra, raspr, vybves, vpv, dmi, cvp); }
	else if ((k<cved) && (k>=cve)) { srk=VasFrayObsh(laefm, lavo, por, u); //for (j=0; j<u; j++) cout << "k = " << k << "\tVas Fra Ob = " << srk[j] << endl; 
	} sr[k]=srk; } 
q=0; p=0.0; k=0; for (j=0; j<cved; j++) if (fabs(sr[j][k])>lvmi) { p=p+ht; no[q]=j; q++; } 
for (k=0; k<u; k++) { s=0.0; for (j=0; j<q; j++) s=s+sr[no[j]][k]; srzn[k]=s/p; 
cout << "kttk = " << srzn[k] << "\tlaef = " << laefm[k] << "\t"; 
} cout << endl;
//for (k=0; k<q; k++) { w=no[k]; f=1; for (j=0; j<u; j++) if (fabs(sr[w][j])<lvmi) f=-1; if (f>0) { for (j=0; j<u; j++) { srk[j]=sr[w][j]; cout << "nomer = " << w+1 << "\tkttk = " << srk[j] << endl; } cout << endl; } } //j=0; for (k=0; k<cved; k++) { if (k==no[j]) { srk=sr[k]; zapisvfile(srk, u, szfstk); j++; } } 
cout  << "Nomera nenulevyh elementov" << endl; for (k=0; k<q; k++) cout << no[k]+1 << " "; cout << endl << "q = " << q << endl; 
for (k=0; k<cved; k++) { srk=sr[k]; delete []srk; } delete []sr; delete []lavo; delete []no; //k=getchar();
return srzn; }
//Функции расчета КТП твердого каркаса и вспомогательные для них
double opredFiBolN(double x, double y)
{ int n=11, k; double *y1=new double[n], *y01=new double[n], *fib=new double[n], fb01, fb1, fibo; 
if ((!y1) || (!y01) || (!fib)) {cout << "No memory" << endl; k=getchar(); exit(1); } 
k=0;
y1[k]=0.0;    k++; y1[k]=47e-3;  k++; y1[k]=113e-3; k++; y1[k]=195e-3; k++; y1[k]=3e-1;   k++; 
y1[k]=413e-3; k++; y1[k]=563e-3; k++; y1[k]=725e-3; k++; y1[k]=832e-3; k++; y1[k]=932e-3; k++; y1[k]=1e0;
k=0;
y01[k]=0.0;  k++; y01[k]=261e-3; k++; y01[k]=411e-3; k++; y01[k]=526e-3; k++; y01[k]=618e-3; k++;
y01[k]=7e-1; k++; y01[k]=774e-3; k++; y01[k]=837e-3; k++; y01[k]=9e-1;   k++; y01[k]=958e-3; k++; y01[k]=1e0;
fib[0]=0.0; for (k=1; k<n; k++) fib[k]=fib[k-1]+1e-1; 
int f1=1, f2=1, q1=0, q2=0;
for (k=0; k<n; k++) {
	if ((f1>0) && (y<y01[k])) { q1=k; f1=0; break; }
	if ((f2>0) && (y<y1[k])) { q2=k; f2=0; break; } }
if (!q1) q1=1; if (!q2) q2=1;
fb01=fib[q1-1]+(fib[q1]-fib[q1-1])*(y-y01[q1-1])/(y01[q1]-y01[q1-1]);
fb1=fib[q2-1]+(fib[q2]-fib[q2-1])*(y-y1[q2-1])/(y1[q2]-y1[q2-1]);
fibo=fb1+(fb01-fb1)*(x-1e0)/(1e-1-1e0);
delete []y1; delete []y01; delete []fib;
return fibo; }
void NapMasPr()
{ int k=0; Prvo[k]=0.707; k++; Prvo[k]=0.705; k++; Prvo[k]=0.703; k++; Prvo[k]=0.701; k++; Prvo[k]=.699; k++;
Prvo[k]=0.698; k++; Prvo[k]=0.696; k++; Prvo[k]=0.694; k++; Prvo[k]=0.692; k++; Prvo[k]=0.69; k++;
Prvo[k]=0.688; k++; Prvo[k]=0.686; k++; Prvo[k]=0.684; k++; Prvo[k]=0.682; k++; Prvo[k]=0.681; k++;
Prvo[k]=0.68; k++; Prvo[k]=0.677; k++; Prvo[k]=0.674; k++; Prvo[k]=0.676; k++; Prvo[k]=0.678; k++;
Prvo[k]=0.687; k++; Prvo[k]=0.699; k++; Prvo[k]=0.706; k++; Prvo[k]=0.713; k++; Prvo[k]=0.717; k++;
Prvo[k]=0.719; k++; Prvo[k]=0.722; k++; Prvo[k]=0.724; }
double opredPrVozd(double T)
{ int n=dmvtk; double Pr=opredKTPTKTochSha(Prvo, vtetk, T, n);
return Pr; }
double DulnevSloForSha(double m2, double T, double eps, double lavo, double lam1, double srch)
{ double r=srch, d=2.0*r, epf=1e-2; double Nk=pow(pow(m2,2.0)-10.0*m2+9.0,0.5); Nk=(m2+3.0+Nk)/2.0/m2; /*y2=3.3e-3*pow(1-m2,-2./9.); pud=0; hsl=30e-3; rona=0.2; ro1=rona/(1-m2); y2=y2*pow(pud+9.8*ro1*(1-m2)*hsl,1./3.); y1=1e-20; y2=y1;*/
double eta=1e-4, y1=1e-3, y2=y1/pow(eta,0.5); /*y1=(10.0+50.0)*1e-4/2.0;*/
double y3=2.0*pow(Nk-1.0,0.5)/Nk, y4=y3/pow(1.0-m2,1.0/3.0); /*hsh=2e-3/2.0;*/
double hsh=1e-10, fb, fbn, fbnn;
fbn=opredFiBolN(y2,y2/y3); fbnn=opredFiBolNN(lavo/lam1,m2);
if (fbn>1e0) if (fbnn<1e0) fb=fbnn; else fb=0.0; else if (fbnn>1e0) fb=fbn; else 
if (fabs(fbn-fbnn)<epf) fb=(fbn+fbnn)/2.0; else if (fbn>fbnn) fb=fbnn; else fb=fbn;
double A=pow(y2,2.0)-pow(y1,2.0), F=pow(1.0-pow(y2,2.0),0.5), D=pow(1.0-pow(y3,2.0),0.5), E=pow(y4,2.0)-pow(y3,2.0);
double delsrsz=d*(hsh+1e0/Nk), lamszm=MoleSostTeplVozd(T,delsrsz,lavo), epspr=eps/(2.0-eps), sig=5.668e-8;
double lamszl=4.0*epspr*sig*delsrsz*pow(T,3.0), lamsz=lamszm+lamszl, nusz=lamsz/lam1, w=(lavo/lamszm-nusz*D)/(1.0-nusz);
double y12=pow(y1,2.0)/(hsh/2.0+(1.0-hsh/2.0)*fb), DF=fabs(w-D)/fabs(w-F), nug=nusz, numz=MoleSostTeplVozd(T,hsh*r,lavo)/lam1;
double nu2sp=nug; DF=(D-F+w*log(DF))*2.0*nug/(1.0-nug); 
double AF=A/(1.0-hsh/2.0-F+hsh/2.0/numz), ADF=1.0/(AF+DF);
ADF=1.0/(D/pow(y3,2.0)+ADF); ADF=ADF+nu2sp*E+y12;
return (ADF*lam1/pow(y4,2.0)); }
double opredFiBolNN(double nu, double m2)
{ double cb=1e2, ca=-1e1, ra=1e2, ep=tocras, cc, fa, fb, fc; 
int h=0, kit=100;
while ((ra>ep) && (h<kit)) {
cc=(ca+cb)/2.0;
fa=2.0*pow(ca,3.0)-3.0*pow(ca,2.0)+1.0-m2;
fb=2.0*pow(cb,3.0)-3.0*pow(cb,2.0)+1.0-m2;
fc=2.0*pow(cc,3.0)-3.0*pow(cc,2.0)+1.0-m2;
if ((fc*fb>0) && (fa*fc<0)) cb=cc; if ((fc*fa>0) && (fb*fc<0)) ca=cc; 
ra=fabs(fa-fb); h++; }
double fi1=fabs(cc/(2.0-cc)); fi1=pow(fi1,0.5);
return pow(nu,2.0)*(7.5-11.0*nu+4.5*pow(nu,2.0))*(1.0-fi1)+fi1; }
double MoleSostTeplVozd(double T, double de, double lamg)
{double gam=7e0/5e0, H=101325.0, H0=1e5, kB=1.380658e-23, n=H/kB/T, d=(6e-1*21e-2+65e-2*79e-2)*1e-10, sig=pi*pow(d,2.0)/4.0;
double dli=1e0/pow(2e0,0.5)/n/sig, Pr=opredPrVozd(T), a=opredKoefAkkomodLandau(T); 
double cz=8.42e-3, Ty=113.0, la0=cz/H0/(1e0+Ty/T), Kn=la0*H0/H/de, B=4e0*gam/(gam+1e0)*(2e0-a)/a*Kn/Pr;
return lamg/(1e0+B); }
double opredTeploprovVozd(double temp)
{ int n=dmvtk; double ktptsp=opredKTPTKTochSha(ktpvotk, vtetk, temp, n);
return ktptsp; }
double urovPodder(double pro, double ref, double urpo)
{ double fl; int f=1; 
if (pro<0) f=0; else if (pro<(urpo*ref)) f=0; else f=1;
if (f) fl=pro; else fl=0.0; return fl; }
double urovOtsechen(double pro, double ref, double urot)
{ double fl; int f=1; 
if (pro<0.0) f=-1; else if (pro>(urot*ref)) f=-1; else f=1;
if (f>0) fl=pro; else fl=0.0; return fl; }
double ReshNeyaUravMaxKoTepl(double lae, double lan, double po) //lan - КТП непрерывной фазы (воздух), ищем КТП твердого скелета
{ double ladb=pgl, lada=lgl, ladc, ra=1e2, ep=tocras, fa, fb, fc; int h=0; //1
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0;    
fa=lan*po+(1e0-po)*lada-lae;
fb=lan*po+(1e0-po)*ladb-lae;
fc=lan*po+(1e0-po)*ladc-lae;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(fa-fb); h++; } //cout << "lam_max = " << ladc << " laef = " << lae << " lavo = " << lan << endl;
return ladc; }
double ReshNeyaUravMinKoTepl(double lae, double lan, double po)
{ double ladb=pgl, lada=lgl, ladc=1e0, ra=1e2, ep=tocras, fa, fb, fc; int h=0; //2
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2.0;    
fa=lan*lada/(po*lada+(1.0-po)*lan)-lae;
fb=lan*ladb/(po*ladb+(1.0-po)*lan)-lae;
fc=lan*ladc/(po*ladc+(1.0-po)*lan)-lae;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(fa-fb); h++; } //cout << "lam_min = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplLiechtenecker(double lae, double lan, double po)
{ double ladb=pgl, lada=lgl, ladc=1e5, ra=1e5, ep=tocras, fa, fb, fc; int h=0; //3
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2.;
fa=pow(lada,(1.-po))*pow(lan,po)-lae;
fb=pow(ladb,(1.-po))*pow(lan,po)-lae;
fc=pow(ladc,(1.-po))*pow(lan,po)-lae;
if ((fc*fb>0) && (fa*fc<0)) ladb=ladc; if ((fc*fa>0) && (fb*fc<0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Lie = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplMaxwell(double lae, double lan, double po)
{ double ladb=pgl, lada=lgl, ladc=1e0, ra=1e2, ep=tocras, fa, fb, fc; int h=0; //4
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0;
fa=lan*(lada+2e0*lan-2e0*(1e0-po)*(lan-lada))/(lada+2e0*lan+(1e0-po)*(lan-lada))-lae;
fb=lan*(ladb+2e0*lan-2e0*(1e0-po)*(lan-ladb))/(ladb+2e0*lan+(1e0-po)*(lan-ladb))-lae;
fc=lan*(ladc+2e0*lan-2e0*(1e0-po)*(lan-ladc))/(ladc+2e0*lan+(1e0-po)*(lan-ladc))-lae;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Maxw = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplBurgFrik(double lae, double lan, double po)
{ int h=0, k; double ladb=pgl, lada=lgl, ladc=1e0, ra=1e2, ep=tocras, fa, fb, fc, f[3], fba, fbb, fbc; 
if (!f) { cout << "No memory" << endl; getchar(); exit(1); } f[0]=1e0/8e0; f[1]=f[0]; f[2]=1e0-f[0]-f[1]; //5
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0;
fba=0.0; fbb=0.0; fbc=0.0;
for (k=0; k<3; k++) {
     fba=fba+pow(1e0+((lada/lan-1e0)*f[k]),(-1e0));
     fbb=fbb+pow(1e0+((ladb/lan-1e0)*f[k]),(-1e0));
     fbc=fbc+pow(1e0+((ladc/lan-1e0)*f[k]),(-1e0)); }
fba=fba/3e0; fbb=fbb/3e0; fbc=fbc/3e0;
fa=lan*(1e0+(1e0-po)*(fba*lada/lan-1e0))/(1e0+(1e0-po)*(fba-1e0))-lae;
fb=lan*(1e0+(1e0-po)*(fbb*ladb/lan-1e0))/(1e0+(1e0-po)*(fbb-1e0))-lae;
fc=lan*(1e0+(1e0-po)*(fbc*ladc/lan-1e0))/(1e0+(1e0-po)*(fbc-1e0))-lae;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Bur_Fri = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplShumanVoss(double lae, double lan, double po)
{ double ladb=pgl, lada=lgl, ladc=1e5, ra=1e5, ep=tocras, fa, fb, fc; //6
double lamaa, lamab, lamac, p; int h=0;
p=PoiskPShumanVoss(po); //cout << "Shum Voss p = " << p << endl;
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2.0;
lamaa=lan*lada/(lan+p*(lan-lada));
lamaa=lamaa*(1.0+p*(1.0+p)*(lan-lada)*log(lan*(1.0+p)/p/lada)/(lan+p*(lan-lada)));
fa=lan*pow(po,3.0)+(1-pow(po,3.0))*lamaa-lae;
lamab=lan*ladb/(lan+p*(lan-ladb));
lamab=lamab*(1.0+p*(1.0+p)*(lan-ladb)*log(lan*(1.0+p)/p/ladb)/(lan+p*(lan-ladb)));
fb=lan*pow(po,3.0)+(1.0-pow(po,3.0))*lamab-lae;
lamac=lan*ladc/(lan+p*(lan-ladc));
lamac=lamac*(1.0+p*(1.0+p)*(lan-ladc)*log(lan*(1.0+p)/p/ladc)/(lan+p*(lan-ladc)));
fc=lan*pow(po,3.0)+(1.0-pow(po,3.0))*lamac-lae;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Shu_Vos = " << ladc << endl;
return ladc; }
double PoiskPShumanVoss(double po)
{ double pb=pgl, pa=lgl, pc, fa, fb, fc, ra=1e5, ep=tocras; int h=0;
while ((ra>ep) && (h<hko)) {
pc=(pa+pb)/2.0;
fa=(pow(pa,2.0)+pa)*log(1.0+1.0/pa)-pa-po;
fb=(pow(pb,2.0)+pb)*log(1.0+1.0/pb)-pb-po;
fc=(pow(pc,2.0)+pc)*log(1.0+1.0/pc)-pc-po;
if ((fc*fb>0.0) && (fa*fc<0.0)) pb=pc; if ((fc*fa>0.0) && (fb*fc<0.0)) pa=pc;
ra=fabs(pa-pb); h++; } return pc; }
double ReshNeyaUravKoTeplGoringChirchillSpher(double lae, double lan, double po)
{ double ladb=pgl, lada=lgl, ra=1e4, ep=tocras, vd=1e0-po, ladc, fa, fb, fc; int h=0; //7
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0; 
fa=lae/lan-(2e0+lada/lan-2e0*vd*(1e0-lada/lan))/(2e0+lada/lan+vd*(1e0-lada/lan));
fb=lae/lan-(2e0+ladb/lan-2e0*vd*(1e0-ladb/lan))/(2e0+ladb/lan+vd*(1e0-ladb/lan));
fc=lae/lan-(2e0+ladc/lan-2e0*vd*(1e0-ladc/lan))/(2e0+ladc/lan+vd*(1e0-ladc/lan));
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Gor_Chi_Sph = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplGoringChirchillCyl(double lae, double lan, double po)
{ double ladb=pgl, ladc, lada=lgl, ra=1e2, ep=tocras, vd=1e0-po, fa, fb, fc; int h=0; //8
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0; 
fa=((1e0+lada/lan)/(1e0-lada/lan)-vd)/((1e0+lada/lan)/(1e0-lada/lan)+vd)-lae/lan;
fb=((1e0+ladb/lan)/(1e0-ladb/lan)-vd)/((1e0+ladb/lan)/(1e0-ladb/lan)+vd)-lae/lan;
fc=((1e0+ladc/lan)/(1e0-ladc/lan)-vd)/((1e0+ladc/lan)/(1e0-ladc/lan)+vd)-lae/lan;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Gor_Chi_Cyl = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplBruggman(double lae, double lan, double po)//Модель Бруггмана-Ханаи - только для круглых частиц
{ double ladb=pgl, lada=lgl, ladc, ra=1e4, ep=tocras, fa, fb, fc; int h=0; //9
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2.0;
fa=po-(lada-lae)*pow((lan/lae),(1.0/3.0))/(lada-lan);
fb=po-(ladb-lae)*pow((lan/lae),(1.0/3.0))/(ladb-lan);
fc=po-(ladc-lae)*pow((lan/lae),(1.0/3.0))/(ladc-lan);
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Brug = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplMaxwAcken(double lae, double lan, double po)
{ double ladb=pgl, lada=lgl, ra=1e5, ep=tocras, vd=1.0-po, fa, fb, fc, ladc; int h=0; //10
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2.0;
fa=(1.0-lan/lada)/(2.0*lan/lada+1.0);
fa=(1.0+2.0*vd*fa)/(1.0-vd*fa)-lae/lan;
fb=(1.0-lan/ladb)/(2.0*lan/ladb+1.0);
fb=(1.0+2.0*vd*fb)/(1.0-vd*fb)-lae/lan;
fc=(1.0-lan/ladc)/(2.0*lan/ladc+1.0);
fc=(1.0+2.0*vd*fc)/(1.0-vd*fc)-lae/lan;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Maxw_Aken = " << ladc << endl;
return ladc; }
double KoTeplFrickEllipsoidChasti(double lae, double lan, double po)
{ double ladb=pgl, lada=lgl, ra=1e4, ep=tocras, s=0.0; //11
double AsRa=99e-2, b=2.0, c=b, a=AsRa*b, sc, sa, sb, M, ladc, fa, fb, fc;
double p0, e, lb, xa[3], labc[3], pea[3], peb[3], pec[3];
int h=0, nko=hko, k;
if (a>b) { e=pow(a,2.0)-pow(b,2.0); p0=pow(e,0.5)/a; M=-2.0*p0+log(fabs((1.0+p0)/(1.0-p0))); M=M/pow(e,(3.0/2.0)); }
else { e=pow(b,2.0)-pow(a,2.0); p0=pow(e,0.5)/a; M=-p0+atan(p0); M=M*2.0/pow(e,(3.0/2.0)); }
lb=2.0/a/pow(b,2.0)-M; lb=lb/2.0;
labc[0]=M; labc[1]=lb; labc[2]=lb;
for (k=0; k<3; k++)
    xa[k]=2.0/a/b/c/labc[k]-1.0;
while ((ra>ep) && (h<nko)) {
ladc=(lada+ladb)/2.0;
sa=0.0; sb=0.0; sc=0.0;
for (k=0; k<3; k++) {
    pea[k]=(1.0+xa[k])/(xa[k]+lan/lada); sa=sa+pea[k];
    peb[k]=(1.0+xa[k])/(xa[k]+lan/ladb); sb=sb+peb[k];
    pec[k]=(1.0+xa[k])/(xa[k]+lan/ladc); sc=sc+pec[k]; }
fa=lada+po*(lan-lada)*sa/3.0-lae; //объемная доля растворенного вещества (эллипсоидов)
fb=ladb+po*(lan-ladb)*sb/3.0-lae;
fc=ladc+po*(lan-ladc)*sc/3.0-lae;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Fri_Ellips = " << ladc << endl;
return ladc; }
double opredRdeltaSrVzaPro(double por, double laef, double lavo, int vyb)
{ double L=1e0, Pa=1e1, Pdelt=1e2, delb=PoiskDelta(por,L), a=delb/Pa; //12, 13 - грубая и точная оценки соот-но
double la1, delm, lm, fb, la1u, epf=1e-2, fbn=0.0, fbnn=0.0, ep=1e-3, x, y;
la1=raschTeplVzaimoPronKomp(delb,L,laef,lavo); //поиск КТП  твердого каркаса - грубая оценка
delm=L/Pdelt; //дельта малая - толщина трещин
lm=2.0*L-delm; //эль малая
x=2.0*delb/pow(pi,1.0/2.0)/lm; y=a/delb;
fbn=opredFiBolN(x,y);
if (la1>ep) { fbnn=lavo/la1; fbnn=opredFiBolNN(fbnn, por); } else fbnn=0.0;
if (fbn>1.0) if (fbnn<1.0) fb=fbnn; else fb=0.0; else if (fbnn>1.0) fb=fbn; else 
if (fabs(fbn-fbnn)<epf) fb=(fbn+fbnn)/2.0; else if (fbn>fbnn) fb=fbnn; else fb=fbn; //lam2 - пора, lam1 - твердый каркас //cout << "fb = " << fb << "\tx = " << x << "\ty = " << y << endl;
la1u=utochLamTverKark(delm,a,fb,lm,la1,delb,lavo,laef,L); //уточнение КТП твердого каркаса - тонкая оценка
if (vyb<13) { //cout << "lam_grub_otsen = " << la1 << endl; 
return la1; } else { //cout << "lam_tonk_ots = " << la1u << endl; 
return la1u; } }
double PoiskDelta(double m2, double L) //поиск дельты большой
{ double ca=0.0, cb=1e0, cc, ra=fabs(ca-cb), ep=tocras, fa, fb, fc; int h=0, nit=hko; //cc - относительный размер бруса
while ((ra>ep) && (h<nit)) {
cc=(ca+cb)/2.0;
fa=2.0*pow(ca,3.0)-3.0*pow(ca,2.0)+1.0-m2; 
fb=2.0*pow(cb,3.0)-3.0*pow(cb,2.0)+1.0-m2; 
fc=2.0*pow(cc,3.0)-3.0*pow(cc,2.0)+1.0-m2; 
if ((fc*fb>0.0) && (fa*fc<0.0)) cb=cc; if ((fc*fa>0.0) && (fb*fc<0.0)) ca=cc; 
ra=fabs(ca-cb); h++; } return (cc*L); }
double raschTeplVzaimoPronKomp(double del, double L, double laef, double la_v) //взаимпопроникающие компоненты
{ double R1a, R1b, R1c, R2a, R2b, R2c, fa, fb, fc, R=1.0/laef/L, R3=1.0/la_v/del;
double R4=L/la_v/pow((L-del),2.0), la1b=pgl, la1a=lgl, la1c, ra=1e4, ep=tocras; int h=0;
while ((ra>ep) && (h<hko)) {
la1c=(la1a+la1b)/2.0;    
R1a=L/la1a/pow(del,2.0);
R1b=L/la1b/pow(del,2.0);
R1c=L/la1c/pow(del,2.0);
R2a=1.0/la1a/(L-del);
R2b=1.0/la1b/(L-del);
R2c=1.0/la1c/(L-del);
fa=1.0/R-1.0/R1a-2.0/(R2a+R3)-1.0/R4;
fb=1.0/R-1.0/R1b-2.0/(R2b+R3)-1.0/R4;
fc=1.0/R-1.0/R1c-2.0/(R2c+R3)-1.0/R4;
if ((fc*fb>0.0) && (fa*fc<0.0)) la1b=la1c; 
if ((fc*fa>0.0) && (fb*fc<0.0)) la1a=la1c; 
ra=fabs(la1a-la1b); h++; }
return la1c; }
double utochLamTverKark(double delm, double a, double phib, double lm, double la1, double delb, double lavo, double laef, double L) //структура с переменным сечением компонент
{ double R=1.0/laef/L, b=2.0*a/pow(pi,0.5), delsh=delm/L, r=2.0*delb/pow(pi,0.5), la2=lavo, y=b/r, laz=lavo;
double R6=delm/2.0/laz/(pow(delb,2.0)-pow(a,2.0)), R4, R3, phim, la1a, la1b, ra, ep, la1c, fa, fb, fc;
double lamsha, lamshb, lamshc, R1a, R1b, R1c, Rp0a, Rp0b, Rp0c, R7a, R7b, R7c, R2a, R2b, R2c; 
double R5a, R5b, R5c, Rdsa, Rdsb, Rdsc, Rdelmina, Rdelmaxa, Rdelminb, Rdelmaxb, Rdelminc, Rdelmaxc; //laz - КТП компоненты, заполняющей трещины (воздух)
R4=L/la2/(pow((L-delb),2.0)); 
R3=1.0/delb/la2;
phim=phib*pow(y,2e0)/(phib-pow(y,2e0));
la1a=lgl; la1b=la1+pgl; ra=1e5; ep=tocras; 
int h=0, nit=hko;
while ((ra>ep) && (h<nit)) {
la1c=(la1a+la1b)/2e0;
lamsha=la1a; //термосопротивление шейки
R1a=delm/2.0/pow(a,2e0)/lamsha; 
Rp0a=lm/la1a/pi/pow(b,2e0); //полное сопротивление квадратного бруса
R7a=2e0*Rp0a*phim;
lamshb=la1b; //термосопротивление шейки
R1b=delm/2e0/pow(a,2e0)/lamshb;
Rp0b=lm/la1b/pi/pow(b,2e0); //полное сопротивление квадратного бруса
R7b=2e0*Rp0b*phim;
lamshc=la1c; //термосопротивление шейки
R1c=delm/2e0/pow(a,2e0)/lamshc;
Rp0c=lm/la1c/pi/pow(b,2e0); //полное сопротивление квадратного бруса
R7c=2e0*Rp0c*phim;
R2a=2e0*Rp0a*phib; //тепловое сопротивление восьмой части центрального бруса
R2b=2e0*Rp0b*phib;
R2c=2e0*Rp0c*phib;
R5a=1e0/(L-delb)/la1a;
R5b=1e0/(L-delb)/la1b;
R5c=1e0/(L-delb)/la1c;
Rdsa=1e0/(R1a+R2a)+1e0/(R6+R7a); Rdsa=1e0/Rdsa;
Rdsb=1e0/(R1b+R2b)+1e0/(R6+R7b); Rdsb=1e0/Rdsb;
Rdsc=1e0/(R1c+R2c)+1e0/(R6+R7c); Rdsc=1e0/Rdsc;
Rdelmina=RdeltaMin(R1a,R6,la1a,delb,L,delsh);
Rdelmaxa=RdeltaMax(la1a,a,delb,L,delsh,R1a,R6);
Rdelminb=RdeltaMin(R1b,R6,la1b,delb,L,delsh);
Rdelmaxb=RdeltaMax(la1b,a,delb,L,delsh,R1b,R6);
Rdelminc=RdeltaMin(R1c,R6,la1c,delb,L,delsh);
Rdelmaxc=RdeltaMax(la1c,a,delb,L,delsh,R1c,R6);
fa=1e0/Rdsa+1e0/R4+2e0/(R3+R5a)-1e0/R;
if ((Rdsa>=Rdelmina) && (Rdsa<=Rdelmaxa))
fa=1e0/Rdsa+1e0/R4+2e0/(R3+R5a)-1e0/R;
fb=1e0/Rdsb+1e0/R4+2e0/(R3+R5b)-1e0/R;
if ((Rdsb>=Rdelminb) && (Rdsb<=Rdelmaxb))
fb=1e0/Rdsb+1e0/R4+2e0/(R3+R5b)-1e0/R;
fc=1e0/Rdsc+1e0/R4+2e0/(R3+R5c)-1e0/R; 
if ((Rdsc>=Rdelminc) && (Rdsc<=Rdelmaxc))
fc=1e0/Rdsc+1e0/R4+2e0/(R3+R5c)-1e0/R;  
if ((fc*fb>0.0) && (fa*fc<0.0)) la1b=la1c; 
if ((fc*fa>0.0) && (fb*fc<0.0)) la1a=la1c; 
ra=fabs(la1a-la1b); h++; } 
return la1c; }
double RdeltaMin(double R1, double R6, double la1, double delb, double L, double delsh)
{ double rdma=(1.0/R1+1.0/R6);
return (1.0/rdma+(L-delsh/2.0)/la1/pow(delb,2.0)); }
double RdeltaMax(double la1, double a, double delb, double L, double dsh, double R1, double R6)
{ double R2=L*(1.0-dsh/2.0)/la1/pow(a,2.0), R7=L*(1-dsh/2.0)/la1/(pow(delb,2.0)-pow(a,2.0)), rdma=1.0/(R1+R2)+1.0/(R6+R7);
return (1.0/rdma); }
double ReshNeyaUravKoTeplGoringChirchillCylPlotSl(double laef, double lavo, double po, double srch)
{ double x0=srch/2e0, lan=lavo, ra; int h, k; //cout << "x0 = " << x0 << "\tlaef = " << laef << "\tlavo = " << lan << endl; //x0=(2.+0.7)*1e-3/2./2.;
double kBa, kBb, kBc, CBa, CBb, CBc, fa, fb, fc, kC, ladb, ladc, lada, ep; //14
kC=1e0; ladb=pgl; lada=lgl; ra=1e5; ep=tocras; h=0; 
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0;
kBa=lada/kC/(lan-lada);
kBa=sign(kBa)*pow(fabs(kBa),1e0/3e0); 
kBa=fabs(kBa);
CBa=pow(kBa,2e0)-kBa*x0+pow(x0,2e0);
CBa=fabs(CBa);
CBa=log(pow(CBa,1e0/2e0)/fabs(kBa+x0));
CBa=CBa+pow(3e0,1e0/2e0)*atan((2e0*x0-kBa)/pow(3e0*kBa,1e0/2e0));
CBa=CBa-pow(3e0,1e0/2e0)*atan(1e0/pow(3e0,1e0/2e0));
CBa=CBa*pi*lan*lada/6e0/(lan-lada)/kC/kBa;
CBa=CBa+lada*(1e0-pi*pow(x0,2e0)/4e0);
kBb=ladb/kC/(lan-ladb);
kBb=sign(kBb)*pow(fabs(kBb),1e0/3e0); 
kBb=fabs(kBb);
CBb=pow(kBb,2e0)-kBb*x0+pow(x0,2e0);
CBb=fabs(CBb);
CBb=log(pow(CBb,1e0/2e0)/fabs(kBb+x0));
CBb=CBb+pow(3e0,1e0/2e0)*atan((2e0*x0-kBb)/pow(3e0*kBb,1e0/2e0));
CBb=CBb-pow(3e0,1e0/2e0)*atan(1e0/pow(3e0,1e0/2e0));
CBb=CBb*pi*lan*ladb/6e0/(lan-ladb)/kC/kBb;
CBb=CBb+ladb*(1e0-pi*pow(x0,2e0)/4e0);
kBc=ladc/kC/(lan-ladc);
kBc=sign(kBc)*pow(fabs(kBc),1e0/3e0);
kBc=fabs(kBc);
CBc=pow(kBc,2e0)-kBc*x0+pow(x0,2e0);
CBc=fabs(CBc);
CBc=log(pow(CBc,1e0/2e0)/fabs(kBc+x0));
CBc=CBc+pow(3e0,1e0/2e0)*atan((2e0*x0-kBc)/pow(3e0*kBc,1e0/2e0));
CBc=CBc-pow(3e0,1e0/2e0)*atan(1e0/pow(3e0,1e0/2e0));
CBc=CBc*pi*lan*ladc/6e0/(lan-ladc)/kC/kBc;
CBc=CBc+ladc*(1e0-pi*pow(x0,2e0)/4e0);
fa=CBa-laef;
fb=CBb-laef;
fc=CBc-laef;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "la Gor Cher Cyl Plot Sloi = " << ladc << endl;
return ladc; }
double MetodHashinaShtrikmanaMin(double laef, double lan, double por) //оценка вилки - для изотропных в среднем образцов
{ double lada=lgl, ladb=pgl, ra=1e6, ep=tocras, ta, tb, tc, fa, fb, fc, ladc; int h=0; //15
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2.0;
ta=1/(lada-lan)+(1.0-por)/3.0/lada; ta=lan+por/ta;
tb=1/(ladb-lan)+(1.0-por)/3.0/ladb; tb=lan+por/tb;
tc=1/(ladc-lan)+(1.0-por)/3.0/ladc; tc=lan+por/tc;
fa=ta-laef;
fb=tb-laef;
fc=tc-laef;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "la Hash Shtrik min = " << ladc << endl;
return ladc; }
double MetodHashinaShtrikmanaMax(double laef, double lan, double por)
{ double lada=lgl, ladb=pgl, ra=1e6, ep=tocras, ta, tb, tc, fa, fb, fc, ladc; int h=0; //16
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0;    
ta=1e0/(lada-lan)+por/3e0/lada; ta=lada+(1e0-por)/ta;
tb=1e0/(ladb-lan)+por/3e0/ladb; tb=ladb+(1e0-por)/tb;
tc=1e0/(ladc-lan)+por/3e0/ladc; tc=ladc+(1e0-por)/tc;
fa=ta-laef;
fb=tb-laef;
fc=tc-laef;
if ((fc*fb>0.0) && (fa*fc<0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "la Hash Shtrik max = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplSwift(double lae, double lan)
{ double lada, ladb, ra, ep, ladc, fa, fb, fc; int h=0; //17
lada=lgl; ladb=pgl; ra=1e5; ep=tocras; fa=0.0; fb=0.0; fc=0.0;
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0;
fa=(log(fabs(lan/lada)))*pow((lan/lada-1),(-2e0));
fa=pow((lan/lada-1e0),(-1e0))-fa;
fa=0.577*pi*fa+0.093-lae/lan;
fb=(log(fabs(lan/ladb)))*pow((lan/ladb-1e0),(-2e0));
fb=pow((lan/ladb-1e0),(-1e0))-fb;
fb=0.577*pi*fb+0.093-lae/lan;
fc=(log(fabs(lan/ladc)))*pow((lan/ladc-1e0),(-2e0));
fc=pow((lan/ladc-1e0),(-1e0))-fc;
fc=0.577*pi*fc+0.093-lae/lan;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "la Swift = " << ladc << endl;
return ladc; }
double *opredTvChaKoeTepSrInSpSha(int vyb, double po, double *tem, double *laefm, int lear, double *srla, double *lavoz, double srch) //определение КТП твердой части иными способами
{ double lam, lavo, laef, uo=urovPod(po), up=1e0, ep=tocras, m2=po; int k, f; //if (vyb<2) { for (k=0; k<lear; k++) cout << "tem = " << tem[k] << " lae = " << laefm[k] << " lavo = " << lavoz[k] << endl; }
for (k=0; k<lear; k++) {
lavo=lavoz[k]; laef=laefm[k];
if (vyb == 1) lam=ReshNeyaUravMaxKoTepl(laef, lavo, po); //максимальный КТП - слои идут параллельно тепловому потоку
if (vyb == 2) lam=ReshNeyaUravMinKoTepl(laef, lavo, po); //минимальный КТП - слои идут перпендикулярно тепловому потоку
if (vyb == 3) lam=ReshNeyaUravKoTeplLiechtenecker(laef, lavo, po);
if (vyb == 4) lam=ReshNeyaUravKoTeplMaxwell(laef, lavo, po);
if (vyb == 5) lam=ReshNeyaUravKoTeplBurgFrik(laef, lavo, po);
if (vyb == 6) lam=ReshNeyaUravKoTeplShumanVoss(laef, lavo, po);
if (vyb == 7) lam=ReshNeyaUravKoTeplGoringChirchillSpher(laef, lavo, po);
if (vyb == 8) lam=ReshNeyaUravKoTeplGoringChirchillCyl(laef, lavo, po);
if (vyb == 9) lam=ReshNeyaUravKoTeplBruggman(laef, lavo, po);
if (vyb == 10) lam=ReshNeyaUravKoTeplMaxwAcken(laef, lavo, po);
if (vyb == 11) lam=KoTeplFrickEllipsoidChasti(laef, lavo, po); 
if (vyb == 12) lam=opredRdeltaSrVzaPro(po, laef, lavo, vyb);
if (vyb == 13) lam=opredRdeltaSrVzaPro(po, laef, lavo, vyb);
if (vyb == 14) lam=ReshNeyaUravKoTeplGoringChirchillCylPlotSl(laef, lavo, po, srch);
if (vyb == 15) lam=MetodHashinaShtrikmanaMin(laef, lavo, po);
if (vyb == 16) lam=MetodHashinaShtrikmanaMax(laef, lavo, po); 
if (vyb == 17) lam=ReshNeyaUravKoTeplSwift(laef, lavo);
if (vyb == 18) lam=SysPlastShakhPor(laef, lavo, po);
if (vyb == 19) lam=MetodStarostina(laef, lavo, po);
if (vyb == 20) lam=MetodRusselNeprTverTelo(laef, lavo, po);
if (vyb == 21) lam=MetodRusselNeprVozd(laef, lavo, po);
srla[k]=lam; }
srla=proverkakvi(srla, laefm, lavoz, up, uo, lear); //if (vyb==6) for (k=0; k<lear; k++) cout << vyb << "\tk = " << k << "\tsr_lam_tt = " << srla[k] << "\n";
return srla; }
double sign(double x)
{ if (x>0.0) return 1e0; else if (x<0.0) return (-1e0); else if (x==0.0) return 0.0; }
double *DulnevKoefTepShamot(int vyb, double por, double *tem, int n, double *srk, double *laefm, double *lavo, double srp, int koel, double *srra, double *rpr) //программа определяет КТП твердого каркаса шамота методом Дульнева
{ int k=0; //cout << "\tvyb = " << vyb << "\tkoel = " << koel; 
if (vyb<koel) return BezUchetaRaspedPorSha(por, tem, laefm, lavo, vyb, n, srk, srp); //1 - Кубы, адиаб., 2 - Кубы, изот., 3 - По Одолевскому 1, 4 - По Одолевскому 2, 5 - Цилиндры, адиаб., 6 - Цилиндры, изот.
else return UchetRapsredPorPoRazmSha(por, tem, laefm, lavo, vyb, n, srk, k, srra, rpr); //7 - Кубы, адиаб., 8 - Кубы, изотерм., 9 - Цилиндры, адиаб., 10 - Цилиндры, изот., 11 - По Одолевскому 1, 12 - По Одолевскому 2
} 
double *DulnevKoefTepVerm(int vyb, double por, double *tem, int n, double *srk, double *laefm, double *lavo, double srp, int koel, double *srra, double *raspr, int vybves, int vpv, int dmi, int nac) //определяет КТП твердого каркаса методом Дульнева
{ int k=0; double *ukaz=NULL; //cout << "\tvyb = " << vyb << "\tkoel = " << koel; 
if ((vyb-nac)<koel) ukaz=BezUchetaRaspedPorSha(por, tem, laefm, lavo, vyb-nac, n, srk, srp); //1 - Кубы, адиаб., 2 - Кубы, изот., 3, 4 - По Одолевскому 1, 2, 5 - Цилиндры, адиаб., 6 - Цилиндры, изот.
else ukaz=UchetRapsredPorPoRazmVer(por, tem, laefm, lavo, vyb-nac, n, srk, k, srra, raspr, vybves, vpv, dmi); //7 - Кубы, адиаб., 8 - Кубы, изотерм., 9 - Цилиндры, адиаб., 10 - Цилиндры, изот., 11 - По Одолевскому 1, 12 - По Одолевскому 2
return ukaz; } 
double *BezUchetaRaspedPorSha(double por, double *tem, double *lameff, double *lamvoz, int vy, int n, double *lamvy, double srp) //без учета распределения пор
{ double m2=por, l=srp, lb=l/pow(por, 1e0/3e0), lavo, laef, uo=urovPod(por), up=1e0; int k; //cout << " l = " << l << "\tpor = " << por << "\tlb = " << lb << endl;;
double *lamadi=new double[n], *lamizo=new double[n], *lamodom=new double[n], *lamodos=new double[n];
double *lamcyladi=new double[n], *lamcylizo=new double[n], *lamcylper=new double[n];
if ((!lamadi) || (!lamizo) || (!lamodom) || (!lamodos) || (!lamcyladi) || (!lamcylizo) || (!lamcylper)) 
{ cout << "No memory" << endl; k=getchar(); exit(1); } //cout << "n = " << n << "\tl = " << l << endl; 
for (k=0; k<n; k++) { 
laef=lameff[k]; lavo=lamvoz[k]; lamvy[k]=0.0;
lamadi[k]=DulKoefTep1Adi(lavo, por, laef); //адиабатное разбиение, полость - прямоугольный параллелепипед
lamizo[k]=DulKoefTep1Izoterm(lavo, por, laef); //изотермическое разбиение, полость - прямоугольный параллелепипед
lamodom[k]=DulKoefTep1OdolevMatr(lavo, por, laef); //по Одолевскому (матрицы)
lamodos[k]=DulKoefTep1OdolevStatSm(lavo, por, laef); //по Одолевскому (стат. смесь)
lamcyladi[k]=DulKoefTep1AdiCyl(lavo, por, laef, lb); //адиабатное разбиение, полость - цилиндр
lamcylizo[k]=DulKoefTep1IzotermCyl(lavo, por, laef, lb); //изотермическое разбиение, полость - цилиндр 
lamcylper[k]=DulKoefTep1CylPer(lavo, por, laef, lb); //полость - цилиндр, тепловой поток направлен перпендикулярно оси
} 
for (k=0; k<n; k++) { laef=lameff[k]; lavo=lamvoz[k];
if (!vy) lamvy[k]=lamadi[k]; //куб, адиаб. //22
if (vy == 1) lamvy[k]=lamizo[k]; //куб, изотерм. //23
if (vy == 2) lamvy[k]=lamodom[k]; //по Одолевскому (матрицы) //24
if (vy == 3) lamvy[k]=lamodos[k]; //по Одолевскому (стат. смесь) //25
if (vy == 4) lamvy[k]=lamcyladi[k]; //цилиндр., адиаб. //26
if (vy == 5) lamvy[k]=lamcylizo[k]; //цилиндр., изотерм. //27
if (vy == 6) lamvy[k]=lamcylper[k]; //цилиндр, перпендикулярно тепловому потоку //28
}
lamvy=proverkakvi(lamvy, lameff, lamvoz, up, uo, n);//if (!vy) for (k=0; k<n; k++) cout << "Dul KTP Kub Adiab 22 = " << lamvy[k] << endl; //if (vy==1) for (k=0; k<n; k++) cout << "Dul KTP Kub Izoterm 23 = " << lamvy[k] << endl; //if (vy==2) for (k=0; k<n; k++) cout << "Dul KTP Odolev 24 = " << lamvy[k] << endl; if (vy==3) for (k=0; k<n; k++) cout << "Dul KTP Odolev 25 = " << lamvy[k] << endl; if (vy==4) for (k=0; k<n; k++) cout << "Dul Adiad Cyl 26 = " << lamvy[k] << endl; if (vy==5) for (k=0; k<n; k++) cout << "Dul Izoterm Cyl 27 = " << lamvy[k] << endl; if (vy==6) for (k=0; k<n; k++) cout << "Lam Cyl Per 28 = " << lamvy[k] << endl;
delete []lamadi; delete []lamizo; delete []lamodom; delete []lamodos; 
delete []lamcyladi; delete []lamcylizo; delete []lamcylper;
return lamvy; }
double *UchetRapsredPorPoRazmSha(double po, double *T, double *lameff, double *lamvoz, int vy, int n, double *lamvy, int no, double *srra, double *rpr) //учет распределения пор
{ int k, m=dmsrps, j;
double *vo=new double[m], *pw=new double[m], mx=0.0, lcub, lb, ep=tocras, w=1e2, ot=1e0/3e0; 
if ((!vo) || (!pw)) { cout << "No memory" << endl; k=getchar(); exit(1); } 
if (po<1e0) w=1e0;
for (k=0; k<m; k++) { pw[k]=0.0; vo[k]=0.0; }
j=0; for (k=0; k<m; k++) { if (rpr[k]>ep) {
    pw[j]=rpr[k]*po/w; //объемная доля поры заданного размера в полном (во всем) объеме
    vo[j]=pow(srra[k], 3e0)/pw[j]; //все поры - кубы
	mx=mx+vo[j]; j++; } }
lcub = pow(mx, ot); //оценка максимального значения l - размера образца
lb = lcub / pow(po, ot); //оценка размера поры
double lavo, lame, lamadit, lamizot, lamodot1, lamodot2, lamcyladit, lamcylizot, lamcylpert, uo=urovPod(po), up=1e0;
double *lamcubadia=new double[n], *lamcubizo= new double[n], *lamodo1=new double[n], *lamodo2=new double[n];
double *lamcyladi=new double[n], *lamcylizo=new double[n], *lamcylper=new double[n], m2;
if ((!lamcubadia) || (!lamcubizo) || (!lamodo1) || (!lamodo2) || (!lamcyladi) || (!lamcylizo) || (!lamcylper)) 
{ cout << "No memory" << endl; getchar(); exit(1); }
for (k=0; k<n; k++) {
    lavo=lamvoz[k]; lame=lameff[k]; lamadit=lame; lamizot=lame; lamodot1=lame; 
	lamodot2=lame; lamcyladit=lame; lamcylizot=lame; lamcylpert=lame;
    for (j=0; j<m; j++) { if (pw[j]>0.0) {
lamadit = DulKoefTep1Adi(lavo, pw[j], lamadit); 
lamizot = DulKoefTep1Izoterm(lavo, pw[j], lamizot);
lamodot1 = DulKoefTep1OdolevMatr(lavo, pw[j], lamodot1); 
lamodot2 = DulKoefTep1OdolevStatSm(lavo, pw[j], lamodot2);
lamcyladit = DulKoefTep1AdiCyl(lavo, pw[j], lamcyladit, lb); 
lamcylizot = DulKoefTep1IzotermCyl(lavo, pw[j], lamcylizot, lb); 
lamcylpert = DulKoefTep1CylPer(lavo, pw[j], lamcylpert, lb); } }
lamcubadia[k]=lamadit; lamcubizo[k]=lamizot; lamodo1[k]=lamodot1; 
lamodo2[k]=lamodot2; lamcyladi[k]=lamcyladit; lamcylizo[k]=lamcylizot; 
lamcylper[k]=lamcylpert; }
for (k=0; k<n; k++) { lame=lameff[k]; lavo=lamvoz[k];
if (vy == 7) lamvy[k]=lamcubadia[k]; //куб, адиаб. //29
if (vy == 8) lamvy[k]=lamcubizo[k]; //куб., изотерм. //30
if (vy == 9) lamvy[k]=lamcyladi[k]; //цил., адиаб. //31
if (vy == 10) lamvy[k]=lamcylizo[k]; //цил., изотерм. //32
if (vy == 11) lamvy[k]=lamodo1[k]; //по Одолевскому (матрицы) //33
if (vy == 12) lamvy[k]=lamodo2[k]; //по Одолевскому (стат. смесь) //34
if (vy == 13) lamvy[k]=lamcylper[k]; //цил., перп. тепл. потоку //35
} lamvy=proverkakvi(lamvy, lameff, lamvoz, up, uo, n); //if (vy == 13) for (k=0; k<n; k++) cout << " lam voz = " << lamvoz[k] << endl; if (vy == 7) for (k=0; k<n; k++) cout << "Dul KTP Adi Uchet (" << k << ") = " << lamvy[k] << endl;
delete []lamcubadia; delete []lamcubizo; delete []lamodo1; delete []lamodo2; 
delete []lamcyladi; delete []lamcylizo; delete []lamcylper; delete []pw; delete []vo;
return lamvy; }
double *UchetRapsredPorPoRazmVer(double po, double *T, double *lameff, double *lamvoz, int vy, int n, double *lamvy, int no, double *srra, double *raspr, int vybves, int pkv, int dmi) //учет распределения пор
{ int k, j, m=dmi, vn=8; double t=0.0, ht=1e0, w=1e2, ep=tocras, ot=1e0/3e0; 
t=0.0; for (k=0; k<m; k++) t=t+raspr[k]; 
for (k=0; k<m; k++) raspr[k]=raspr[k]/t; 
t=0.0; for (k=0; k<m; k++) t=t+raspr[k]; //cout << "m = " << m << "\tpo = " << po << "\tt = " << t << endl;
double *vo=new double[m], *pw=new double[m], mx=0.0, lcub=0.0, lb=0.0;
if ((!vo) || (!pw)) { cout << "No memory" << endl; k=getchar(); exit(1); } 
for (k=0; k<m; k++) { pw[k]=0.0; vo[k]=0.0; }
if (po<1e0) w=1e0;
j=0; t=0.0; for (k=0; k<m; k++) { if (raspr[k]>0.0) {
    pw[j]=raspr[k]*po/w; //объемная доля поры заданного размера в полном (во всем) объеме
    vo[j]=pow(srra[k], 3e0)/pw[j]; //все поры - кубы
	mx=mx+vo[j]; j++; t=t+ht; } } //if (vy==vn) { for (k=0; k<m; k++) cout << "raspr = " << raspr[k] << "\tsrra = " << srra[k] << "\t"; cout << endl; for (k=0; k<m; k++) cout << "pw = " << pw[k] << "\tvo = " << vo[k] << "\t"; cout << endl; }
lcub = pow(mx, ot)/t; //оценка максимального значения l - размера образца
lb = lcub / pow(po, ot); //оценка размера поры //cout << "lcub = " << lcub << "\tlb = " << lb << "\tpo = " << po << endl;
double uo=urovPod(po), up=1e0, lavo=0.0, lame=0.0;
double *lamcubadia=new double[n], *lamcubizo=new double[n], *lamodo1=new double[n], *lamodo2=new double[n];
double *lamcyladi=new double[n], *lamcylizo=new double[n], *lamcylper=new double[n];
if ((!lamcubadia) || (!lamcubizo) || (!lamodo1) || (!lamodo2) || (!lamcyladi) || (!lamcylizo) || (!lamcylper)) 
{ cout << "No memory" << endl; k=getchar(); exit(1); } 
for (k=0; k<n; k++) {
	lame=lameff[k]; lavo=lamvoz[k]; lamcubadia[k]=lame; lamcubizo[k]=lame;
	lamodo1[k]=lame; lamodo2[k]=lame; lamcyladi[k]=lame; lamcylizo[k]=lame; lamcylper[k]=lame; //cout << "ktp vozd = " << lavo << "\tte = " << T[k] << endl;
    for (j=0; j<m; j++) {
if (pw[j]>0.0) {
lamcubadia[k] = DulKoefTep1Adi(lavo, pw[j], lamcubadia[k]);
lamcubizo[k] = DulKoefTep1Izoterm(lavo, pw[j], lamcubizo[k]);
lamodo1[k] = DulKoefTep1OdolevMatr(lavo, pw[j], lamodo1[k]); 
lamodo2[k] = DulKoefTep1OdolevStatSm(lavo, pw[j], lamodo2[k]);
lamcyladi[k] = DulKoefTep1AdiCyl(lavo, pw[j], lamcyladi[k], lb); 
lamcylizo[k] = DulKoefTep1IzotermCyl(lavo, pw[j], lamcylizo[k], lb); 
lamcylper[k] = DulKoefTep1CylPer(lavo, pw[j], lamcylper[k], lb); } } }
for (k=0; k<n; k++) { 
if (vy == 7) lamvy[k]=lamcubadia[k]; //куб, адиаб. //29
if (vy == 8) lamvy[k]=lamcubizo[k]; //куб., изотерм. //30
if (vy == 9) lamvy[k]=lamodo1[k]; //по Одолевскому (матрицы) //31
if (vy == 10) lamvy[k]=lamodo2[k]; //по Одолевскому (стат. смесь) //32
if (vy == 11) lamvy[k]=lamcyladi[k]; //цил., адиаб. //33
if (vy == 12) lamvy[k]=lamcylizo[k]; //цил., изотерм. //34
if (vy == 13) lamvy[k]=lamcylper[k]; } //цил., перп. тепл. потоку //35 //cout << "uo = " << uo << "\tup = " << up << "\tpor = " << po << "\tm = " << m << endl; //if (vy == 13) for (k=0; k<n; k++) cout << "Lam po Cyl Perpen 35 = " << lamvy[k] << endl; //cout << endl;
lamvy=proverkakvi(lamvy, lameff, lamvoz, up, uo, n); //if (vy == 8) { for (k=0; k<n; k++) cout << "Lam Kub Izot 30 = " << lamcubizo[k] << "\tlam = " << lamvy[k] << endl; cout << endl; } //if (vy == 13) { for (k=0; k<n; k++) cout << "Lam Cyl Per 35 = " << lamcylper[k] << "\tlam = " << lamvy[k] << endl; cout << endl; } //if (vy == 9) for (k=0; k<n; k++) cout << "Lam po Odolevskomu Matr 31 = " << lamvy[k] << endl; if (vy == 10) for (k=0; k<n; k++) cout << "Lam po Odolevskomu Stat Sm 32 = " << lamvy[k] << endl;//if (vy == 11) for (k=0; k<n; k++) cout << "Lam po Cyl Adiab 33 = " << lamvy[k] << endl; if (vy == 12) for (k=0; k<n; k++) cout << "Lam po Cyl Izoterm 34 = " << lamvy[k] << endl; //if (vy == 13) for (k=0; k<n; k++) cout << "Lam po Cyl Perpen 35 = " << lamvy[k] << endl;
delete []lamcubadia; delete []lamcubizo; delete []lamodo1; delete []lamodo2; 
delete []lamcyladi; delete []lamcylizo; delete []lamcylper; delete []pw; delete []vo;
return lamvy; }
double DulKoefTep1Adi(double lam2, double m2, double lame) //lam1 - КТП твердого каркаса, lam2 - КТП воздуха
{ double lamb=pgl, lama=lgl, ep=tocras, ra=1e0, fa, fb, fc, nua, nub, nuc, lamc; //22
int k=0, nit=hko;
while ((ra>ep) && (k<nit)) {
lamc=(lama+lamb)/2e0; 
nua=lam2/lama; nub=lam2/lamb; nuc=lam2/lamc;
fa=(nua-(nua-1e0)*pow(m2,1e0/3e0)*(1e0-pow(m2,2e0/3e0)))/(nua-pow(m2,1e0/3e0)*(nua-1e0))-lame/lama;
fb=(nub-(nub-1e0)*pow(m2,1e0/3e0)*(1e0-pow(m2,2e0/3e0)))/(nub-pow(m2,1e0/3e0)*(nub-1e0))-lame/lamb;
fc=(nuc-(nuc-1e0)*pow(m2,1e0/3e0)*(1e0-pow(m2,2e0/3e0)))/(nuc-pow(m2,1e0/3e0)*(nuc-1e0))-lame/lamc;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } //cout << "lamc = " << lamc << "\tfc = " << fc << "\tra = " << ra << "\tk = " << k << "\n"; 
return lamc; }
double DulKoefTep1Izoterm(double lam2, double m2, double lame)
{ double lamb=pgl, lama=lgl, ep=tocras, ra=1e0, fa, fb, fc, lamc, nua, nub, nuc; //23 //30
int k=0, nit=hko;
while ((ra>ep) && (k<nit)) {
lamc=(lama+lamb)/2e0;
nua=lam2/lama; nub=lam2/lamb; nuc=lam2/lamc;
fa=(1e0+(nua-1e0)*(pow(m2,2e0/3e0)))/(1e0+pow(m2,2e0/3e0)*(nua-1e0))*(1e0-pow(m2,1e0/3e0))-lame/lama;
fb=(1e0+(nub-1e0)*(pow(m2,2e0/3e0)))/(1e0+pow(m2,2e0/3e0)*(nub-1e0))*(1e0-pow(m2,1e0/3e0))-lame/lamb;
fc=(1e0+(nuc-1e0)*(pow(m2,2e0/3e0)))/(1e0+pow(m2,2e0/3e0)*(nuc-1e0))*(1e0-pow(m2,1e0/3e0))-lame/lamc;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } //cout << "lamc = " << lamc << "\tfc = " << fc << "\tra = " << ra << "\tk = " << k << "\n"; 
return lamc; }
double DulKoefTep1OdolevMatr(double lam2, double m2, double lame) //матричная гетерогенная система - одна фаза образует связную матрицу при любой объемной концентрации этой фазы, система имеет включения в виде кубов, центры которых образуют простую кубическую решетку, ребра параллельны
{ double lamb=pgl, lama=lgl, ep=tocras, ra=1e4, lamc; int k=0, nit=hko; //lam2 - КТП воздуха
double nua, nub, nuc, fa, fb, fc;
while ((ra>ep) && (k<nit)) {
lamc=(lama+lamb)/2e0;
nua=1.0-lama/lam2;
nub=1.0-lamb/lam2;
nuc=1.0-lamc/lam2;
fa=1.0-(1.0-m2)/(1.0/nua-m2/3.0)-lame/lam2;
fb=1.0-(1.0-m2)/(1.0/nub-m2/3.0)-lame/lam2;
fc=1.0-(1.0-m2)/(1.0/nuc-m2/3.0)-lame/lam2;
if ((fc*fb>0) && (fa*fc<0)) lamb=lamc; if ((fc*fa>0) && (fb*fc<0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } 
return lamc; }
double DulKoefTep1OdolevStatSm(double lam2, double m2, double lame)
{ double lamb=pgl, lama=lgl, ep=tocras, ra=1e4, lamc; int k=0, nit=hko; //нет регулярной структуры, статистическая смесь, частицы распределены хаотически
double v1, v2, nua, nub, nuc, fa, fb, fc; //lam2 - КТП воздуха, ищем КТП твердого скелета
while ((ra>ep) && (k<nit)) {
lamc=(lama + lamb)/2.0;
v1=m2; v2=1.0-m2;
nua=((3.0*v1-1.0)*lama+(3.0*v2-1.0)*lam2)/4.0;
nub=((3.0*v1-1.0)*lamb+(3.0*v2-1.0)*lam2)/4.0;
nuc=((3.0*v1-1.0)*lamc+(3.0*v2-1.0)*lam2)/4.0;
fa=nua+pow(pow(nua,2.0)+lama*lam2/2.0,0.5)-lame; //lam2 - КТП в порах
fb=nub+pow(pow(nub,2.0)+lamb*lam2/2.0,0.5)-lame; //lame - ЭКТП
fc=nuc+pow(pow(nuc,2.0)+lamc*lam2/2.0,0.5)-lame;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } 
return lamc; }
double DulKoefTep1AdiCyl(double lam2, double m2, double lame, double d)
{ double fra=0.9, lamb=pgl, lama=lgl, ep=tocras, h=d, ra=1e0, lamc, h1;
double d1, F1, R1a, R1b, R1c, R2, F12, R3a, R3b, R3c, F, R, fa, fb, fc;
int k=0, nit=hko;
while ((ra>ep) && (k<nit)) {
    lamc=(lama+lamb)/2.0;
    h1=fra*d;
    d1=d*pow(4.0*m2/pi/fra, 0.5);
    F1=(pi/4.0)*pow(d1, 2.0); 
    R1a=(h-h1)/2.0/F1/lama;
    R1b=(h-h1)/2.0/F1/lamb;
    R1c=(h-h1)/2.0/F1/lamc;
    R2=(h-h1)/2.0/F1/lam2;
    F12=(pi/4.0)*(pow(d,2.0)-pow(d1,2.0)); 
    R3a=h/F12/lama;
    R3b=h/F12/lamb;
    R3c=h/F12/lamc;
    F=(pi/4.0)*pow(d,2.0); 
    R=h/lame/F;
    fa=R3a*(2.0*R1a+R2)/(2.0*R1a+R2+R3a)-R;
    fb=R3b*(2.0*R1b+R2)/(2.0*R1b+R2+R3b)-R;
    fc=R3c*(2.0*R1c+R2)/(2.0*R1c+R2+R3c)-R;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } 
return lamc; }
double DulKoefTep1IzotermCyl(double lam2, double m2, double lame, double d)
{ double fra = 0.9, lamb = pgl, lama = lgl, ep = tocras, h = d, ra=1e0, lamc;
double h1, d1, F1, R1a, R1b, R1c, R2, F12, R3a, R3b, R3c, F, R, fa, fb, fc;
int k=0, nit=hko;
while ((ra>ep) && (k<nit)) {
    lamc=(lama+lamb)/2e0;
    h1=fra*d;
    d1=d*pow(4.0*m2/pi/fra,0.5);
    F=(pi/4.0)*pow(d,2.0); 
    R1a=(h-h1)/2e0/F/lama;
    R1b=(h-h1)/2e0/F/lamb;
    R1c=(h-h1)/2e0/F/lamc;
    F1=(pi/4e0)*pow(d1,2e0); 
    R2=h1/F1/lam2;
    F12=(pi/4e0)*(pow(d,2e0)-pow(d1,2e0)); 
    R3a=h1/F12/lama;
    R3b=h1/F12/lamb;
    R3c=h1/F12/lamc;
    F=(pi/4e0)*pow(d,2e0); 
    R=h/lame/F;
    fa=2e0*R1a+R3a*R2/(R2+R3a)-R;
    fb=2e0*R1b+R3b*R2/(R2+R3b)-R;
    fc=2e0*R1c+R3c*R2/(R2+R3c)-R;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } 
return lamc; }
double DulKoefTep1CylPer(double lavo, double por, double laef, double lb) 
{ double ra=1e0, lamb=pgl, lam1=lavo, lama=lgl, dt=1e0, h=1e0, Qo=laef*h*dt, rb=pow(por/pi,1e0/2e0)*lb; //lam1 - КТП воздуха, lam2 - КТП ТТ
double lamc, bb, bm, int1, ab, ep=tocras, am, Q1, Q2, fa, fb, fc; //тепловой поток перпендикулярен оси цилиндра
int k=0, nit=hko; 
while ((ra>ep) && (k<nit)) {
    lamc=(lama+lamb)/2e0;
	bb=dt*lam1*lama*h*rb; 
	bm=2e0*(lama-lam1)*rb; am=lam1*lb;
    if (fabs(am)>=fabs(bm)) { 
		ab=pow(pow(am,2e0)-pow(bm,2e0),0.5); 
		int1=2e0*(atan((am+bm)/ab)-atan(bm/ab))/ab;}
	else { 
		ab=pow(pow(bm,2e0)-pow(am,2e0),0.5); 
		int1=(log(fabs(am+bm-ab)/fabs(am+bm+ab))-log(fabs(bm-ab)/fabs(bm+ab)))/ab; }
    Q1=bb*(pi/2e0)/bm-am*bb/bm*int1; Q1=2e0*Q1; 
	Q2=lama*h*(lb-2e0*rb)/lb; 
	fa=Qo-Q1-Q2;
    bb=dt*lam1*lamb*h*rb; 
	bm=2e0*(lamb-lam1)*rb;
	if (fabs(am)>=fabs(bm)) { 
		ab=pow(pow(am,2e0)-pow(bm,2e0),0.5); 
		int1=2e0*(atan((am+bm)/ab)-atan(bm/ab))/ab; }
    else { ab=pow(pow(bm,2e0)-pow(am,2e0),0.5); 
	int1=(log(fabs(am+bm-ab)/fabs(am+bm+ab))-log(fabs(bm-ab)/fabs(bm+ab)))/ab; }
    Q1=bb*(pi/2e0)/bm-am*bb/bm*int1; Q1=2e0*Q1; 
	Q2=lamb*h*(lb-2e0*rb)/lb; fb=Qo-Q1-Q2;
    bb=dt*lam1*lamc*h*rb; 
	bm=2e0*(lamc-lam1)*rb;
	if (fabs(am)>=fabs(bm)) { 
		ab=pow(pow(am,2e0)-pow(bm,2e0),0.5); 
		int1=2e0*(atan((am+bm)/ab)-atan(bm/ab))/ab; }
	else { ab=pow(pow(bm,2e0)-pow(am,2e0),0.5); 
	int1=(log(fabs(am+bm-ab)/fabs(am+bm+ab))-log(fabs(bm-ab)/fabs(bm+ab)))/ab; }
    Q1=bb*(pi/2e0)/bm-am*bb/bm*int1; Q1=2e0*Q1; 
	Q2=lamc*h*(lb-2e0*rb)/lb; fc=Qo-Q1-Q2;
    if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
	k++; ra=fabs(lama-lamb); } 
return lamc; }
double *DulnevZern(double por, double *te, double *srk, int na, double *laefm, double *lavo, double wsio, double walo, double wmgo, double srch, int n, double tena, double dtosc, int dm, double *kuscv, double *tkuscv, double *stchv)
{ int k, f; double uo, up=1e0, epsil, ts, la_voz, la_e, lam, m2, koscve=1e0, stchve;
uo=urovPod(por);
for (k=0; k<n; k++) {
ts=te[k];
stchve=opredKTPTKTochSha(stchv, te, ts, n);
epsil=stchve*koscve; 
stch[k]=epsil; 
la_voz=opredKTPTKTochSha(lavo, te, ts, n);
la_e=opredKTPTKTochSha(laefm, te, ts, n);
lam=opredDulnLam1(por, ts, epsil, la_voz, la_e, srch); 
srk[k]=lam; }
srk=proverkakvi(srk, laefm, lavo, up, uo, n); //for (k=0; k<n; k++) cout << "laefm = " << laefm[k] << "\tlavo = " << lavo[k] << "\tlam = " << srk[k] << endl; 
return srk; }
double opredDulnLam1(double po, double T, double eps, double lavo, double lae, double srch)
{ double ladb=1e5, lada=tocras, ra=1e2, ep=tocras, ladc, fa, fb, fc; int h=0, kit=hko;
while ((ra>ep) && (h<kit)) {
ladc=(lada+ladb)/2e0;
fa=DulnevSloForVer(po, T, eps, lavo, lada, srch); fa=fa-lae;
fb=DulnevSloForVer(po, T, eps, lavo, ladb, srch); fb=fb-lae;
fc=DulnevSloForVer(po, T, eps, lavo, ladc, srch); fc=fc-lae;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(fa-fb); h++; } return ladc; }
double DulnevSloForVer(double m2, double T, double eps, double lavo, double lam1, double srch)
{ double r=srch/2.0, d=srch, epf=1e-2; 
double Nk=pow(pow(m2,2.0)-10.0*m2+9.0,0.5); Nk=(m2+3.0+Nk)/2.0/m2; /*y2=3.3e-3*pow(1-m2,-2./9.); pud=0; hsl=30e-3; rona=0.2; ro1=rona/(1-m2); y2=y2*pow(pud+9.8*ro1*(1-m2)*hsl,1./3.); y1=1e-20; y2=y1;*/
double eta=1e-4, y1=1e-3, y2=y1/pow(eta,0.5); /*y1=(10+50)*1e-4/2;*/
double y3=2.0*pow(Nk-1.0,0.5)/Nk, y4=y3/pow(1.0-m2,1.0/3.0); /*hsh=2e-3/2;*/
double hsh=1e-3, fb, fbn=0.0, fbnn=0.0;
fbn=opredFiBolN(y2,y2/y3); fbnn=opredFiBolNN(lavo/lam1,m2);
if (fbn>1.0) if (fbnn<1.0) fb=fbnn; else fb=0.0; else if (fbnn>1.0) fb=fbn; else 
if (fabs(fbn-fbnn)<epf) fb=(fbn+fbnn)/2.0; else if (fbn>fbnn) fb=fbnn; else fb=fbn; 
double A=pow(y2,2.0)-pow(y1,2.0), F=pow(1.0-pow(y2,2.0),0.5), D=pow(1.0-pow(y3,2.0),0.5), E=pow(y4,2.0)-pow(y3,2.0);
double delsrsz=d*(hsh+1.0/Nk), lamszm=MoleSostTeplVozd(T,delsrsz,lavo), epspr=eps/(2.0-eps), sig=5.668e-8;
double lamszl=4.0*epspr*sig*delsrsz*pow(T,3.0), lamsz=lamszm+lamszl, nusz=lamsz/lam1, w=(lavo/lamszm-nusz*D)/(1.0-nusz);
double y12=pow(y1,2.0)/(hsh/2.0+(1.0-hsh/2.0)*fb), DF=fabs(w-D)/fabs(w-F), nug=nusz, numz=MoleSostTeplVozd(T,hsh*r,lavo)/lam1;
double nu2sp=nug; DF=(D-F+w*log(DF))*2.0*nug/(1.0-nug); double AF=A/(1.0-hsh/2.0-F+hsh/2.0/numz), ADF=1.0/(AF+DF);
ADF=1.0/(D/pow(y3,2.0)+ADF); ADF=ADF+nu2sp*E+y12; return (ADF*lam1/pow(y4,2.0)); }
double *DopFunctOpredKTPTverChas(double *srk, double *tem, int n, double por, double *laefm, int p, int vfv, double *lavo, double r)
{ if (!p) srk=MetodDulnevaSigalovoyPolyDispVer(tem, n, por, laefm, vfv, r);
else if (p==1) srk=MetodDulnevaSigalovoyDopVer(tem, n, por, laefm, vfv, r);
else if (p==2) srk=MetodDulnevaSigalovoyBezIspVer(tem, n, por, laefm, vfv, r); 
else if (p==3) srk=VasFrayObsh(laefm, lavo, por, n);
else { cout << "Irregular number!"; p=getchar(); exit(1); }
double up=1e0, uo=urovPod(por); srk=proverkakvi(srk, laefm, lavo, up, uo, n);
return srk; }
double *MetodDulnevaSigalovoyPolyDispVer(double *tem, int lete, double m22, double *laefm, int vfv, double r)
{ double mg=26e-2, m20=porex /*общая пористость*/, m2=m20-m22; //m2 - межзерновая пористость (внешняя), m22 - внутренняя пористость
if (m2<mg) m2=0.365; //межзерновая пористость по Нижегородову
double A=pow(74e-2/(1e0-m2),1e0/3e0); 
double gu=9.81, hv=3e1*1e-3, rov=125.0, sigb=2.156e5 /*предел прочности при сжатии материала частиц, Па*/, ka=7e0/5e0, T, ladc, ladb, lada; 
double H0=101325.0, H=1e5, kB=1.380658e-23, d=(6e-1*21e-2+65e-2*79e-2)*1e-10, sigse=pi*pow(d,2e0)/4e0, E=8.428e6 /*модуль Юнга вермикулита*/, laef; 
double akk=9e-1 /*коэффициент аккомодации*/, ko1=0.2595/(1-0.2595), ko2=0.2595-ko1, fp=0.2595/m2-(ko1*m2+ko2), ep, uo, up=1e0, ra;
int k, h=0, kit=hko; double *ladi=new double[lete], n, dliSvoPro, cPr, N, dem1, dem2, dem4, B, lamg1, lamg2, lamg4, lamg5, izco;
double sigm1, sigm2, sigm3, sigm3a, sigm3b, sigm3c, sigm4, sigm5, dem5, lam0, sigmr, debo;
double rr=0.0, fa=0.0, fb=0.0, fc=0.0, ch1, ch2, ch3, ch123, ch4a, ch4b, ch4c, gi, ch5, ch6;
if (!ladi) { cout << "No memory" << endl; k=getchar(); exit(1); } 
uo=urovPod(m22);
for (k=0; k<lete; k++) {
ladb=1e5; lada=1e-5; ra=1e2; ep=1e-7; h=0;
while ((ra>ep) && (h<kit)) {
ladc=(lada+ladb)/2e0; T=tem[k];
laef=opredKTPTKTochSha(laefm, tem, T, lete);
izco=2e0*(5.67e-2)*pow(stch[k],2e0)*3e0/4e0;
    lam0=opredTeploprovVozd(T); /*температура в К*/
    n=H/kB/T;
    dliSvoPro=1e0/pow(2e0,1e0/2e0)/n/sigse; //длина свободного пробега молекулы воздуха
    cPr=opredPrVozd(T-tem0); //число Прандтля
    N=pow(pow(m2,2e0)-1e1*m2+9.0,0.5); N=(m2+3e0+N)/2e0/m2;
    debo=rov*gu*hv; /*удельная нагрузка*/
    dem1=r*(A-pi/4e0); /*dem1=2*r*(A-pi/4);*/ if (dem1<0.0) dem1=1e0;
    dem2=7e0*r*pow(debo,7e0/18e0)*pow(sigb,2e0)/pow(N,2e0)/pow(E,8e0/9e0);
    dem4=2e0*r*(pow(2e0,1e0/2e0)*A-1e0); if (dem4<0.0) dem4=1e0; dem5=2e0*pow(2e0,1e0/2e0)*r*A;
    akk=opredKoefAkkomodLandau(T); 
    B=4e0*ka/(ka+1e0)/cPr*((2e0-akk)/akk)*H0*dliSvoPro;
    lamg1=lam0/(1e0+B/H/dem1); lamg2=lam0/(1e0+B/H/dem2); 
	lamg4=lam0/(1e0+B/H/dem4); lamg5=lam0/(1e0+B/H/dem5);
    sigm1=lamg1*pi*r*(A*log(A/(A-1e0))-1e0);
    sigm2=lamg2*pow(N,1e0/2e0)*pow(E,4e0/9e0)*pow(debo,1e0/18e0)*r/3.2/pow(sigb,1e0/2e0);
    sigm2=sigm2*fp; sigm3=(2e-2)*pow(N,1e0/2e0)*pow(E,8e0/9e0)*pow(debo,11e0/18e0)*r/pow(sigb,3e0/2e0);
    sigm3a=lada*sigm3; sigm3b=ladb*sigm3; sigm3c=ladc*sigm3;
    sigm4=lamg4*pi*r*(pow(2e0,1e0/2e0)*A-1e0)/2e0; //большая пора - цилиндр
    if (m2>=74e-2) { sigm5=lamg5*2e0*r*(2e0*pow(A,2e0)-pi)/pow(2e0,1e0/2e0)/A; rr=1e0+pow(A,2e0)/pi; }
    else { sigm5=0.0; rr=1e0+pow(pow(2e0,1e0/2e0)*A-1e0,2e0)/2e0; }
    gi=izco*pow(T/1e2,3e0); 
    sigmr=pi*pow(r,2e0)*gi*rr;
    fa=(pow(2e0,1e0/2e0)*(sigm1+sigm2+sigm4/2+sigm5+sigmr)+sigm3a)/A/r-laef;
    fb=(pow(2e0,1e0/2e0)*(sigm1+sigm2+sigm4/2+sigm5+sigmr)+sigm3b)/A/r-laef;
    fc=(pow(2e0,1e0/2e0)*(sigm1+sigm2+sigm4/2+sigm5+sigmr)+sigm3c)/A/r-laef;
    if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
	if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc;
	ra=fabs(fa-fb); h++; }
ladi[k]=ladc; } return ladi; }
double *MetodDulnevaSigalovoyDopVer(double *tem, int lete, double m22, double *laefm, int vfv, double r) /*m22 - пористость вермикулита фракции 2-0,7 мм, из статьи "Теплопроводность моно- и полидисперсных материалов"*/
{ double m20=porex /*общая пористость*/, m2=m20-m22, mg=26e-2; if (m2<mg) m2=0.365; //межзерновая пористость по Нижегородову
double A=pow(74e-2/(1e0-m2),1e0/3e0) /*alpha - размер воздушной оболочки (ореола)*/, gu=9.8, hv=30e-3, rov=125., sigb=2.156e5; //предел прочности при сжатии материала частиц, Па
double ka=7e0/5e0, H0=101325.0, H=1e5, kB=1.380658e-23, d=((6e-1)*(21e-2)+(65e-2)*(79e-2))*(1e-10), sigse=pi*pow(d,2.0)/4.0;
double E=8.428e6 /*модуль Юнга вермикулита*/, akk=0.9 /*коэффициент аккомодации*/, fa=0.0, fb=0.0, fc=0.0, up=1.0, uo, izco;
double ko1=0.2595/(1e0-0.2595), ko2=0.2595-ko1, fp=0.2595/m2-(ko1*m2+ko2), lada, ladb, ladc, ra, ep; int k, h, kit=hko;
double *ladi=new double[lete], sigm5, n, dliSvoPro, cPr, N, dem1, dem2, dem3, dem4, dem5, B, lamg1, sigm1, sigm2, sigm3a, sigm3b, sigm3c, sigm4; 
double gi, rr=0.0, sigmr, ch1, ch2, ch3, ch123, ch4a, ch4b, ch4c, ch5, ch6, Zs, Zss, V, G, lamrs, lamrss, X, T, laef, lam0;
uo=urovPod(m22);
for (k=0; k<lete; k++) {
ladb=1e5; lada=1e-5; ra=1e2; ep=1e-7; h=0; 
while ((ra>ep) && (h<kit)) {
ladc=(lada+ladb)/2e0; T=tem[k];
	laef=opredKTPTKTochSha(laefm, tem, T, lete); 
	izco=2e0*(5.67e-2)*pow(stch[k],2e0)*3e0/4e0;
	lam0=opredTeploprovVozd(T); /*температура в К*/
    n=H/kB/T; dliSvoPro=1e0/pow(2.0,0.5)/n/sigse; /*длина свободного пробега молекулы воздуха*/
    cPr=opredPrVozd(T-tem0); /*число Прандтля*/ N=pow(pow(m2,2e0)-10.0*m2+9e0,1e0/2e0); N=(m2+3e0+N)/2e0/m2;
    dem1=2e0*r*(A-pi/4e0); if (dem1<0.0) dem1=2e0*r*(A-2e0/3e0); if (dem1<0.0) dem1=1e0;
	dem3=2e0*r*(pow(2e0,1e0/2e0)*A-1e0); if (dem3<0.0) dem3=1e0; dem5=2.0*pow(2.0,1e0/2e0)*r*A;
    akk=opredKoefAkkomodLandau(T); 
    B=4e0*ka/(ka+1e0)*((2e0-akk)/akk)*H0*dliSvoPro/cPr;
    X=4.45*(A*log(A/fabs(A-1e0))-1e0)/(1e0+B/H/dem1); //X - молекулярный перенос тепла через ореолы вокруг частицы
    Zs=2.23*(pow(2e0,1e0/2e0)*A-1e0)/(1e0+B/H/dem3); //Z, Z' - перенос тепла через поры
    Zss=2.23/(1e0+B/H/dem3)/(pow(2e0,1e0/2e0)*A-1e0);
    V=(pow(A,2e0)-pi/2e0)/(1e0+B/H/dem5)/A; //V - через дополнительные пути (столбы) при пористости более 75 %
    G=izco*pow(T/1e2,3e0);
    lamrs=4.45*G*r*(1+0.5*pow(pow(2e0,1e0/2e0)*A-1e0,2e0))/A;
    lamrss=4.45*G*r*(1+pow(A,2e0)/pi);
    if (m2>=74e-2) { fa=lam0*(X+Zss+V)/A+lamrss+lada-laef; 
	fb=lam0*(X+Zss+V)/A+lamrss+ladb-laef; fc=lam0*(X+Zss+V)/A+lamrss+ladc-laef; }
    else if (m2>=26e-2) { fa=lam0*(X+Zs)/A+lamrs+lada-laef;    
    fb=lam0*(X+Zs)/A+lamrs+ladb-laef; fc=lam0*(X+Zs)/A+lamrs+ladc-laef; }
    if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
	if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc;
	ra=fabs(fa-fb); h++; }
ladi[k]=ladc; } return ladi; }
double *MetodDulnevaSigalovoyBezIspVer(double *tem, int lete, double m22, double *laefm, int vfv, double r) //m22 - пористость вермикулита
{ double m20 = porex /*общая пористость*/, m2=m20-m22, mg=26e-2; if (m2<mg) m2=0.365; //межзерновая пористость по Нижегородову //из статьи "Теплопроводность зернистых систем"
double A=pow(74e-2/(1e0-m2),1e0/3e0) /*alpha - размер воздушной оболочки (ореола)*/, gu=9.8, hv=30e-3, rov=125.0, sigb=2.156e5; //предел прочности при сжатии материала частиц, Па
double ka=7e0/5e0, H0=101325.0, H=1e5, kB=1.380658e-23, d=(6e-1*21e-2+65e-2*79e-2)*1e-10, sigse=pi*pow(d,2e0)/4e0;
double E=8.428e6 /*модуль Юнга вермикулита*/, Gms=20.58e6 /*модуль сдвига*/, akk=9e-1 /*коэффициент аккомодации*/;
double fa=0.0, fb=0.0, fc=0.0, up=1e0, uo;
double ko1=0.2595/(1.0-0.2595), ko2=0.2595-ko1, fp=0.2595/m2-(ko1*m2+ko2), dem1, dem2, dem4, dem5; int h, k, kit=hko;
double ra, ep, lada, ladc, ladb, lam0, n, N, dliSvoPro, B, ch1, ch2, ch3, ch123, ch4a, ch4b, ch4c, gi, ch5, ch6,izco;
double *ladi=new double[lete], cPr, lamg1, sigm1, sigm2, sigm3, sigm4, sigm5, T, debo; 
double sigm3a, sigm3b, sigm3c, rr=0.0, sigmr, laef;
uo=urovPod(m22);
for (k=0; k<lete; k++) {
ladb=1e5; lada=1e-5; ra=1e2; ep=1e-7; h=0; kit=1e3;
while ((ra>ep) && (h<kit)) {
ladc=(lada+ladb)/2e0; T=tem[k]; 
laef=opredKTPTKTochSha(laefm, tem, T, lete);
    lam0=opredTeploprovVozd(T);
	izco=2e0*(5.67e-2)*pow(stch[k],2e0)*3e0/4e0;
    n=H/kB/T;
    dliSvoPro=1e0/pow(2e0,1e0/2e0)/n/sigse; //длина свободного пробега молекулы воздуха
    cPr=opredPrVozd(T-tem0); //число Прандтля
    N=pow(pow(m2,2e0)-1e1*m2+9e0,1e0/2e0); N=(m2+3e0+N)/2e0/m2;
    debo=rov*gu*hv; //удельная нагрузка
    dem1=r*(A-pi/4e0); if (dem1<0.0) dem1=1e0;
    dem2=7e0*r*pow(debo,7e0/18e0)*pow(sigb,1e0/2e0)/pow(N,1e0/2e0)/pow(E,8e0/9e0);
    dem4=2e0*r*(pow(2e0,1e0/2e0)*A-1e0); if (dem4<0.0) dem4=1e0;
    dem5=2e0*pow(2e0,1e0/2e0)*r*A;
    akk=opredKoefAkkomodLandau(T); 
    B=4e0*ka/(ka+1e0)/cPr*((2e0-akk)/akk)*H0*dliSvoPro;
    ch1=4.45*(A*log(A/fabs(A-1e0))-1e0)/(1e0+B/H/dem1);
    ch2=fp*0.44*pow(N,1e0/2e0)*pow(E,4e0/9e0)*pow(debo,1e0/18e0)/pow(sigb,1e0/2e0)/(1e0+B/H/dem2);
    ch3=2.23*(pow(2e0,1e0/2e0)*A-1e0)/(1e0+B/H/dem4);
    ch4a=(2e-2)*fp*lada*pow(N,1e0/2e0)*pow(E,8e0/9e0)*pow(debo,11e0/18e0)/pow(sigb,3e0/2e0)/A;
    ch4b=(2e-2)*fp*ladb*pow(N,1e0/2e0)*pow(E,8e0/9e0)*pow(debo,11e0/18e0)/pow(sigb,3e0/2e0)/A;
    ch4c=(2e-2)*fp*ladc*pow(N,1e0/2e0)*pow(E,8e0/9e0)*pow(debo,11e0/18e0)/pow(sigb,3e0/2e0)/A;
    gi=izco*pow(T/1e2,3e0);
    if ((m2>=2595e-4) && (m2<=74e-2)) {
        ch5=4.45*gi*r*(1.0+pow(pow(2e0,1e0/2e0)*A-1e0,2e0)/2e0)/A;
		ch123=lam0*(ch1+ch2+ch3)/A;
		fa=ch123+ch4a+ch5-laef; fb=ch123+ch4b+ch5-laef; fc=ch123+ch4c+ch5-laef; }
    if (m2>74e-2) {
        ch3=2.23/(1e0+B/H/dem4)/(pow(2e0,1e0/2e0)*A-1e0);
        ch6=(pow(A,2e0)-1.58)/A/(1e0+B/H/dem5);
        ch5=4.45*gi*r*(1e0+pow(A,2e0)/pi)/A;
        fa=(ch1+ch2+ch3+ch6)*lam0/A+ch4a+ch5-laef;
        fb=(ch1+ch2+ch3+ch6)*lam0/A+ch4b+ch5-laef;
		fc=(ch1+ch2+ch3+ch6)*lam0/A+ch4c+ch5-laef; }
    if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
	if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc;
	ra=fabs(fa-fb); h++; }
ladc=urovOtsechen(ladc, laef, uo); 
ladc=urovPodder(ladc, laef, up);
ladi[k]=ladc; } return ladi; }
double opredKoefAkkomodLandau(double T)
{ double PP=6.6260755e-34/2.0/pi, kB=1.380658e-23, NA=6.0221409e23, R=kB*NA, mu=29e-3;
double a=3e-9, gamv=7.0/5.0, m=mu/NA, ro=opredPlotnVozd(T-tem0);
double c=pow(gamv*R*T/mu,0.5), alpha=pow(a*c,2.0); alpha=(kB*T)/alpha; alpha=pow(alpha,3.0/2.0);
alpha=alpha/6.0/ro; alpha=alpha/pow(2.0*pi*m,0.5); return alpha; }
double opredPlotnVozd(double T)
{ double *te=arrTempAir(), *ro=arrPlotnAir(), r=opredKTPTKTochSha(ro, te, T, dmvtk);
return r; }
double *arrTempAir()
{ double *te=new double[dmvtk]; int k=0;
if (!te) { cout << "No memory!"; k=getchar(); exit(1); } else {
te[k]=0.0; k++; te[k]=1e1; k++; te[k]=2e1; k++; te[k]=3e1; k++; te[k]=4e1; k++; te[k]=5e1; k++;
te[k]=6e1; k++; te[k]=7e1; k++; te[k]=8e1; k++; te[k]=9e1; k++; te[k]=1e2; k++; te[k]=12e1; k++;
te[k]=14e1; k++; te[k]=16e1; k++; te[k]=18e1; k++; te[k]=2e2; k++; te[k]=25e1; k++; te[k]=3e2; k++;
te[k]=35e1; k++; te[k]=4e2; k++; te[k]=5e2; k++; te[k]=6e2; k++; te[k]=7e2; k++; te[k]=8e2; k++;
te[k]=9e2; k++; te[k]=1e3; k++; te[k]=11e2; k++; te[k]=12e2; k++; } return te; }
double *arrPlotnAir()
{ double *ar=new double[dmvtk]; int k=0;
if (!ar) { cout << "No memory!"; k=getchar(); exit(1); }
ar[k]=1.293; k++; ar[k]=1.247; k++; ar[k]=1.205; k++; ar[k]=1.165; k++; ar[k]=1.128; k++; ar[k]=1.093; k++;
ar[k]=1.06; k++; ar[k]=1.029; k++; ar[k]=1.0; k++; ar[k]=0.972; k++; ar[k]=0.946; k++; ar[k]=0.898; k++;
ar[k]=0.854; k++; ar[k]=0.815; k++; ar[k]=0.779; k++; ar[k]=0.746; k++; ar[k]=0.674; k++; ar[k]=0.615; k++;
ar[k]=0.566; k++; ar[k]=0.524; k++; ar[k]=0.456; k++; ar[k]=0.404; k++; ar[k]=0.362; k++; ar[k]=0.329; k++;
ar[k]=0.301; k++; ar[k]=0.277; k++; ar[k]=0.257; k++; ar[k]=0.239; k++; return ar; }
double opredUrovPodderM03(double por)
{ int n=3, k; double p, s, *pam=new double[n];
if (!pam) { cout << "No memory!" << endl; k=getchar(); exit(1); }
pam[2]=-8.93227577791347; pam[1]=8.89444783404518; pam[0]=1.06833435021354;
p=0.0; s=0.0; for (k=0; k<n; k++) { s=s+pam[k]*pow(por,p); p=p+1.0; }
delete []pam; return 1.0/(1.0-s*por); }
double *VasFraySha(double *tem, int n)
{ int ke=9, k=0; double *te=new double[ke], *lsk=new double[ke], *vfs=new double[n];
if ((!te) || (!lsk) || (!vfs)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
k=0;	te[k]=1e2; k++; te[k]=18e1; k++; te[k]=22e1; k++; te[k]=3e2; k++; 
		te[k]=4e2; k++; te[k]=6e2;  k++; te[k]=8e2;  k++; te[k]=1e3; k++; te[k]=12e2; k++;
k=0;	lsk[k]=1.52; k++; lsk[k]=2.25; k++; lsk[k]=2.50; k++; lsk[k]=2.85; k++;
		lsk[k]=3.24; k++; lsk[k]=3.28; k++; lsk[k]=3.46; k++; lsk[k]=4.06; k++; lsk[k]=4.98; k++;
for (k=0; k<n; k++) { vfs[k]=0.0; vfs[k]=opredKTPTKTochSha(lsk, te, tem[k], ke); } delete []lsk; delete []te; return vfs; }
double SysPlastShakhPor(double laef, double lam2, double po)
{ double lamb=pgl, lama=lgl, ep=tocras, ra=1e0, fa, fb, fc, lamc; //18
int k=0, nit=hko;
while ((ra>ep) && (k<nit)) {
lamc=(lama+lamb)/2e0;
if (po<=5e-1) {
fa=lam2*(4.0*po/(1.0+lam2/lama)+lama*(1.0-2.0*po)/lam2)-laef;
fb=lam2*(4.0*po/(1.0+lam2/lamb)+lamb*(1.0-2.0*po)/lam2)-laef;
fc=lam2*(4.0*po/(1.0+lam2/lamc)+lamc*(1.0-2.0*po)/lam2)-laef; }
else if (po<1e0) {
fa=lam2*(4.0*(1.0-po)/(1+lam2/lama)+2.0*(po-1.0))-laef;
fb=lam2*(4.0*(1.0-po)/(1+lam2/lamb)+2.0*(po-1.0))-laef;
fc=lam2*(4.0*(1.0-po)/(1+lam2/lamc)+2.0*(po-1.0))-laef; }
if ((fc*fb>0) && (fa*fc<0)) lamb=lamc; if ((fc*fa>0) && (fb*fc<0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } //cout << "la Shakh poryadok = " << lamc << endl;
return lamc; }
double MetodStarostina(double laef, double lam2, double po)
{ double lamb=pgl, lama=lgl, ep=tocras, ra=1e0, fa, fb, fc, lamc; //19
int k=0, nit=hko;
while ((ra>ep) && (k<nit)) {
lamc=(lama+lamb)/2e0;
fa=pow(lama,2.0)*pow(po,2.0/3.0)+(lam2-lama)*lama;
fa=fa/(lama+(pow(po,2.0/3.0)-po)*(lam2-lama))-laef;
fb=pow(lamb,2.0)*pow(po,2.0/3.0)+(lam2-lamb)*lamb;
fb=fb/(lamb+(pow(po,2.0/3.0)-po)*(lam2-lamb))-laef;
fc=pow(lamc,2.0)*pow(po,2.0/3.0)+(lam2-lamc)*lamc;
fc=fc/(lamc+(pow(po,2.0/3.0)-po)*(lam2-lamc))-laef;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } //cout << "la metod Starostina = " << lamc << endl;
return lamc; }
double MetodRusselNeprTverTelo(double laef, double lam2, double po)
{ double lamb=pgl, lama=lgl, ep=tocras, ra=1e0, fa, fb, fc, lamc; //20
int k=0, nit=hko;
while ((ra>ep) && (k<nit)) {
lamc=(lama+lamb)/2e0;
fa=lama*po+(lama/lam2)*(1e0-pow(po,2e0/3e0));
fa=fa/(po-pow(po,2e0/3e0)+(lam2/lama)*(1e0-pow(po,2e0/3e0)+po))-laef;
fb=lamb*po+(lamb/lam2)*(1e0-pow(po,2e0/3e0));
fb=fb/(po-pow(po,2e0/3e0)+(lam2/lamb)*(1e0-pow(po,2e0/3e0)+po))-laef;
fc=lamc*po+(lamc/lam2)*(1e0-pow(po,2e0/3e0));
fc=fc/(po-pow(po,2e0/3e0)+(lam2/lamc)*(1e0-pow(po,2e0/3e0)+po))-laef;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } //cout << "la metod Russelya nepr tver telo = " << lamc << endl;
return lamc; }
double MetodRusselNeprVozd(double laef, double lam2, double po)
{ double lamb=pgl, lama=lgl, ep=tocras, ra=1e0, fa, fb, fc, lamc; //21
int k=0, nit=hko;
while ((ra>ep) && (k<nit)) {
lamc=(lama+lamb)/2e0;
fa=lama*pow(1e0-po,2e0/3e0)+1e0-pow(1e0-po,2e0/3e0);
fa=fa/((lama/lam2)*(pow(1e0-po,2e0/3e0)-1e0+po)+(2e0-pow(1e0-po,2e0/3e0)-po))-laef;
fb=lamb*pow(1e0-po,2e0/3e0)+1e0-pow(1e0-po,2e0/3e0);
fb=fb/((lamb/lam2)*(pow(1e0-po,2e0/3e0)-1e0+po)+(2e0-pow(1e0-po,2e0/3e0)-po))-laef;
fc=lamc*pow(1e0-po,2e0/3e0)+1e0-pow(1e0-po,2e0/3e0);
fc=fc/((lamc/lam2)*(pow(1e0-po,2e0/3e0)-1e0+po)+(2e0-pow(1e0-po,2e0/3e0)-po))-laef;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } //cout << "la metod Russelya nepr vozd = " << lamc << endl;
return lamc; }
double *VasFrayObsh(double *laefm, double *lavo, double po, int n)
{ int k, nit=hko, h=0; double L=1e0, x=Poiskkx((1.0-po)/4.0), db=L*x, hlm=x/(0.5-x), hlb=hlm/(1.0+hlm);
double lamb=pgl, lama=lgl, ep=tocras, ra=1e0, fa, fb, fc, lamc, mu=1e-2, E=8.428e6;
double Pn=125.0*9.8*30e-3, Aa, Ab, Ac, kk=(1.5+2e0)/2.0, kb=(2.2+2.9)/2e0;
double km=(4e0+5e0)/2e0, nugza, nugzb, nugzc, lamz, nu=2e0*(1e0-pow(mu,2e0))/E, kc=(35e-2+45e-2)/2e0;
double rp=(725e-3)*pow(nu*Pn*L/2e0,1e0/3e0), lamka, lamkb, lamkc, lampna, lampnb, lampnc, Q=pow(74e-2/(1e0-po),1e0/3e0);
double hsh=L*km*1e-3, sfk=pi*pow(rp,2e0), lamksha, lamkshb, lamkshc, nuga, nugb, nugc, m20=porex, m2=po, uo=urovPod(po);
double laef, lavoz, *lto=new double[n], up=1e0; if (!lto) { cout << "No memory!" << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) {
	laef=laefm[k]; lavoz=lavo[k]; h=0;
while ((ra>ep) && (h<nit)) {
lamc=(lama+lamb)/2e0;
lamz=lavoz; nugza=lamz/lama; nugzb=lamz/lamb; nugzc=lamz/lamc;
if (Pn<3e5) { lampna=lama*pow(Pn,2e0/3e0)*kc/75.0/Q; 
lampnb=lamb*pow(Pn,2e0/3e0)*kc/75.0/Q; 
lampnc=lamc*pow(Pn,2e0/3e0)*kc/75.0/Q; }
else { lampna=lama*pow(Pn/E,4e0/9e0)*kb/Q; 
lampnb=lamb*pow(Pn/E,4e0/9e0)*kb/Q;
lampnc=lamc*pow(Pn/E,4e0/9e0)*kb/Q; }
lamksha=2e0*pow(2e0,1e0/2e0)*sfk*lama/L/Q/hsh/kk; 
lamkshb=2e0*pow(2e0,1e0/2e0)*sfk*lamb/L/Q/hsh/kk; 
lamkshc=2e0*pow(2e0,1e0/2e0)*sfk*lamc/L/Q/hsh/kk;
lamka=lamksha+lampna;
lamkb=lamkshb+lampnb;
lamkc=lamkshc+lampnc;
Aa=pow(hlb,2e0)*1e3*nugza/4e0/kk/km+lamka/lama;
Ab=pow(hlb,2e0)*1e3*nugzb/4e0/kk/km+lamkb/lamb;
Ac=pow(hlb,2e0)*1e3*nugzc/4e0/kk/km+lamkc/lamc;
nuga=nugza; nugb=nugzb; nugc=nugzc;
fa=1e0/(pow(1e0/hlb,2e0)+Aa)+nuga*pow(1e0-hlb,2e0)+2e0/(1e0+hlm+1e0/(nuga*hlb))-laef/lama;
fb=1e0/(pow(1e0/hlb,2e0)+Ab)+nugb*pow(1e0-hlb,2e0)+2e0/(1e0+hlm+1e0/(nugb*hlb))-laef-lamb;
fc=1e0/(pow(1e0/hlb,2e0)+Ac)+nugc*pow(1e0-hlb,2e0)+2.0/(1e0+hlm+1e0/(nugc*hlb))-laef/lamc;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
h++; ra=fabs(lama-lamb); } 
lto[k]=lamc; }
return lto; }
double Poiskkx(double kx)
{double xa=0.0, xb=5e-1, xc, ra=fabs(xa-xb), ep=tocras, fa, fb, fc; 
int h=0, nit=hko; //xc - относительный размер бруса
while ((ra>ep) && (h<nit)) {
xc=(xa+xb)/2e0;
fa=4.0*pow(xa,3.0)-3.0*pow(xa,2.0)+kx; 
fb=4.0*pow(xb,3.0)-3.0*pow(xb,2.0)+kx; 
fc=4.0*pow(xc,3.0)-3.0*pow(xc,2.0)+kx; 
if ((fc*fb>0.0) && (fa*fc<0.0)) xb=xc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) xa=xc; 
ra=fabs(xa-xb); h++; } return (xc); }
double urovPod(double po)
{ double m2=porex, uo, pormakvi=53e-2, porg=5e-1, porveit=6e-1, pormikvi=36e-2;
if (po>porveit) m2=porex; else if (po>pormikvi) m2=pormakvi; else m2=po;
if (po<porist) uo=opredUrovPodderM03(po); else if (po<porg) uo=1e0/(1e0-po); else uo=1e0/(1e0-m2);
return uo; }
double *proverkakvi(double *ta, double *laefm, double *lavo, double up, double uo, int n)
{ int k, f; double lvmi=1e3, e=1e-3; 
for (k=0; k<n; k++) if (lvmi>lavo[k]) lvmi=lavo[k];
for (k=0; k<n; k++) {
	ta[k]=urovOtsechen(ta[k], laefm[k], uo);
	ta[k]=urovPodder(ta[k], lavo[k], up); }
f=1; for (k=0; k<n; k++) if ((ta[k]<lvmi) && (f>0)) { f=-1; break; }
if (f<0) for (k=0; k<n; k++) ta[k]=0.0; 
return ta; }
double **opredSredKTPTK()
{
	int k=1, f=cemtk, j=0, q=j, w=j, b=11, l=0, fl=1, na=0, *ceti=new int[f]; 
	int nnvyuv=2, kvf=0, jvsv=0, qvmi=0, qvuv=0, d=0;
	int vyfrve=0, vysove=0, vymivmf=1, vyukve=1, c=0, nnvyfv=2, nnvysv=3, nnvmivmf=3;
	double *ttk=NULL, cf=0.0, hf=1e0, e=1e-6, *cet=NULL, tna=2e2+273.15, dt=1e2, *ktptk=NULL;
	double **muv=new double*[b], *tem=new double[f], *cemt=new double[f], **mu=NULL;
	double *nvyfv=new double[nnvyfv], *nvysv=new double[nnvysv], *nvmivmf=new double[nnvmivmf], **muc=new double*[b];
	double *nvyuv=new double[nnvyuv], **mut=new double*[b], *tepv=new double[cemtk], tf=0.0; 
	if ((!nvyfv) || (!nvysv) || (!nvmivmf) || (!nvyuv)) 
	{ cout << "No memory" << endl; k=getchar(); exit(1); } 
	k=0; tem[k]=tna; for (k=1; k<cemtk; k++) tem[k]=tem[k-1]+dt;
	for (k=0; k<cemtk; k++) { cemt[k]=0.0; tepv[k]=0.0; }
	k=0; nvyfv[k]=0; k++; nvyfv[k]=1; //фракции вермикулита
	k=0; nvysv[k]=0; k++; nvysv[k]=1; k++; nvysv[k]=2; //состояния вермикулита
	k=0; nvmivmf[k]=0; k++; nvmivmf[k]=1; k++; nvmivmf[k]=2; //стационарные методы измерений - 2019 и 2020
	k=0; nvyuv[k]=1; k++; nvyuv[k]=2; b=0; na=0; //укладка вермикулита
	for (kvf=na; kvf<nnvyfv; kvf++) {
		vyfrve=nvyfv[kvf];
		for (jvsv=0; jvsv<nnvysv; jvsv++) {
			vysove=nvysv[jvsv];
			if ((!vyfrve) || (vyfrve==2)) {
				for (qvmi=0; qvmi<nnvmivmf; qvmi++) {
					vymivmf=nvmivmf[qvmi]; 
					mu=initarrver(fl, vyfrve, vyukve, vymivmf, vysove);
					k=0; ktptk=mu[k]; k++; ttk=mu[k]; k++; cet=mu[k]; //cout << "\tce = " << cet[0] << endl;
					muv[b]=ktptk; mut[b]=ttk; muc[b]=cet; b++; 
					if (vysove>0) break; } }
			else if (vyfrve==1) {
				for (qvuv=0; qvuv<nnvyuv; qvuv++) {
					vyukve=nvyuv[qvuv];
					mu=initarrver(fl, vyfrve, vyukve, vymivmf, vysove);
					k=0; ktptk=mu[k]; k++; ttk=mu[k]; k++; cet=mu[k]; //cout << "\tce = " << cet[0] << endl;
					muv[b]=ktptk; mut[b]=ttk; muc[b]=cet; b++;
				} } } } //cout << "\tb = " << b << endl;
	for (j=0; j<b; j++) { 
		ttk=mut[j]; ktptk=muv[j]; cet=muc[j]; 
		l=0; tf=cet[l]; c=0; cf=0.0; while (cf<(tf-e)) { cf=cf+hf; c++; } cout << "\tce = " << cf << endl;
		for (q=0; q<c; q++) {
			cout << "ktp_tk = " << ktptk[q] << "\ttem = " << ttk[q] << endl;
			for (k=0; k<cemtk; k++) {
				if (fabs(tem[k]-ttk[q])<=hf) {
					tepv[k]=tepv[k]+ktptk[q];
					cemt[k]=cemt[k]+hf; 
					break; } } } }
		q=0; for (k=0; k<cemtk; k++) {
			cf=cemt[k]; if (cf>e) tepv[k]=tepv[k]/cf; else tepv[k]=0.0; 
			c=0; tf=0.0; while (tf<(cf-e)) { tf=tf+hf; c++; } ceti[k]=c;
			if (c==b) { q++; } }
		for (k=0; k<cemtk; k++) if (tepv[k]>e) cout << "ktp_tk_sr = " << tepv[k] << "\ttem = " << tem[k] << "\tce = " << cemt[k] << endl;
		f=1; double *ktpvy=new double[q], *tevy=new double[q], *cemvy=new double[f]; 
		tf=0.0; for (k=0; k<q; k++) tf=tf+hf; k=0; cemvy[k]=tf;
		j=0; for (k=0; k<cemtk; k++) { c=ceti[k]; if (b==c) { tevy[j]=tem[k]; ktpvy[j]=tepv[k]; j++; } }
	for (k=0; k<b; k++) { cet=muv[k]; if (cet) delete[]cet; cet=mut[k]; if (cet) delete[]cet; cet=muc[k]; if (cet) delete[]cet; }
	if (nvyfv) delete[]nvyfv; if (nvysv) delete[]nvysv;
	if (nvmivmf) delete[]nvmivmf; if (nvyuv) delete[]nvyuv;
	if (muv) delete[]muv; if (mut) delete[]mut; if (muc) delete[]muc; 
	if (tem) delete[]tem; if (cemt) delete[]cemt; if (ceti) delete[]ceti;
	f=3; muv=new double*[f]; k=0; muv[k]=ktpvy; k++; muv[k]=tevy; k++; muv[k]=cemvy;
	return muv;
}