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
using namespace std;
const double pi = acos(-1e0), enu=1e-3, tocrasov=1e-6;
const int dmkoscrt=14;
//-----
int dlar=1, N=6;
double *koefoslab(double, double, double, double *, int, double *);
double epsisred(double, double *, double *, int, double *, double *, int);
double *izmMasChast(double *, double, int, double *);
double trapz(double *, double *, int);
double *Kramers_Kronig(double *, int, double *, char *, char *, char *, int, int);
double *dliny_voln(double *, int);
double *epsilnu(double *, int, double *);
double *Koef_Pogl(double *);
double *PokazPrelomAl(double *, double *, double *, double, int, double *);
double *PokazPrelomPribl(double *, double *, double *, int, double *, double);
void zapisvfile(double *, int, char *, char *);
double *ronusha(double *, int);
double *ronukvi(double *, int);
double reflsha(double);
double reflkvi(double, int);
double reflitom(double);
double *reshMetKram(double **, double *, int);
double *oslaintestepcher(int, int, double *);
double *RasKorOV(double, double, double, double, double, double);
double UsredMasOV(double **, int, int);
void vyvomatr(double **, int, int);
double KoefPoglRossel(double, double *, double *, double *, int);
double *KoefPoglRosselNac(double *, int, int, double, double, double, double *, double *, int, double, double, double *, double *, int, int, int);
double nsredPlank2(double *, double *, double, int);
void osvpamov(char *, char *, char *, char *);
double usredVelichPlank(double *, double *, double *, double, int, double);
double LuchKTPChudnovsky(double *, double, int, double);
double PoiskReflPhi(double, double, double, double);
double ReflSred(double);
double PoiskTetaShtrNach(double, double, double, double);
double PraCha10(double, double, double, double);
double *PoiSpeKoeOtr(double *, int);
double provUgla(double);
double *GRAYDIFFSPEC(int, double *, double *, double *, double *, double **, int *, double *, double *, int);
double *GRAYDIFF(int, double *, double *, double *, double **, int *, double *, double *, int);
double VIEW(int, int, double *, int);
double perpplates(double, double, double, double, double, double, double);
double parlplates(double, double, double, double, double, double, double);
double *GAUSS(double **, double *, int, double *);
double SeryeStenkiRasIzl(double, double, double, double *, double *, double *, double *, double *, int *, int);
double SeryeStenkiRasIzlDifPov(double, double, double, double *, double *, double *, double *, double *, int *, int);
double *SEMIGRAY(int, double *, double **, double **, double **, double ***, int *, double *, double *, double *, double *, int, int);
double opredKTPTKTochSha(double *, double *, double, int);
double bbfn(double);
double *MetodGaussa(double **, double *, int, double *);
double *reshMetObrMatr(double **, double *, int, double *);
char **napolStrok(int, int);
char **napNazvFile(int, int, char *, char);
//-------------------
char **napolStrok()
{	int k=0, dsov=60, n=2; 
char *snmov=new char[dsov], *sfnoov=new char[dsov], **ms=new char*[n];
if ((!snmov) || (!sfnoov) || (!ms)) { cout << "No memory!" << endl; k=getchar(); exit(1); } 
for (k=0; k<dsov; k++) { snmov[k]='\0'; sfnoov[k]='\0'; }
k=0;
snmov[k]='N'; k++; snmov[k]='o'; k++; snmov[k]='_'; k++; snmov[k]='m'; k++; snmov[k]='e'; k++; snmov[k]='m'; k++; 
snmov[k]='o'; k++; snmov[k]='r'; k++; snmov[k]='y'; k++; snmov[k]='!'; k++; snmov[k]='\0';
k=0;
sfnoov[k]='F';  k++; sfnoov[k]='i';  k++; sfnoov[k]='l';  k++; sfnoov[k]='e';  k++; sfnoov[k]='_';  k++; sfnoov[k]='i'; k++;
sfnoov[k]='s';  k++; sfnoov[k]='_';  k++; sfnoov[k]='n';  k++; sfnoov[k]='o';  k++; sfnoov[k]='t';  k++; sfnoov[k]='_'; k++;
sfnoov[k]='o';  k++; sfnoov[k]='p';  k++; sfnoov[k]='e';  k++; sfnoov[k]='n';  k++; sfnoov[k]='!';  k++; sfnoov[k]='\0'; 
k=0; ms[k]=snmov; k++; ms[k]=sfnoov;
return ms;
}
char **napNazvFile(int vyve, int vm, char *snm, char pbf)
{
int k=0, dsov=60; 
char *sndov=new char[dsov], *ssdov=new char[dsov], *szfovk=new char[dsov], *szfov=new char[2*dsov];
char *skpovk=NULL, *sdvovk=NULL, *sppovk=NULL, *skpov=NULL, *sppov=NULL, *sdvov=NULL;
for (k=0; k<dsov; k++) { sndov[k]='\0'; ssdov[k]='\0'; szfovk[k]='\0'; }
k=0; 
sndov[k]='D'; k++; sndov[k]=':'; k++; sndov[k]='\\';k++; sndov[k]='\\'; k++; sndov[k]='_';  k++; sndov[k]='А';  k++;
sndov[k]='с'; k++; sndov[k]='п'; k++; sndov[k]='и'; k++; sndov[k]='р';  k++; sndov[k]='а';  k++; sndov[k]='н';  k++;
sndov[k]='т'; k++; sndov[k]='у'; k++; sndov[k]='р'; k++; sndov[k]='а';  k++; sndov[k]='\\'; k++; sndov[k]='\\'; k++;
sndov[k]='t'; k++; sndov[k]='m'; k++; sndov[k]='p'; k++; sndov[k]='\\'; k++; sndov[k]='\\'; k++; sndov[k]='\0';
k=0;
ssdov[k]='C';  k++; ssdov[k]=':'; k++; ssdov[k]='\\'; k++; ssdov[k]='\\'; k++; ssdov[k]='U';  k++; ssdov[k]='s';  k++;
ssdov[k]='e';  k++; ssdov[k]='r'; k++; ssdov[k]='s';  k++; ssdov[k]='\\'; k++; ssdov[k]='\\'; k++; ssdov[k]='А';  k++;
ssdov[k]='н';  k++; ssdov[k]='д'; k++; ssdov[k]='р';  k++; ssdov[k]='е';  k++; ssdov[k]='й';  k++; ssdov[k]='\\'; k++;
ssdov[k]='\\'; k++; ssdov[k]='D'; k++; ssdov[k]='o';  k++; ssdov[k]='c';  k++; ssdov[k]='u';  k++; ssdov[k]='m';  k++;
ssdov[k]='e';  k++; ssdov[k]='n'; k++; ssdov[k]='t';  k++; ssdov[k]='s';  k++; ssdov[k]='\\'; k++; ssdov[k]='\\'; k++;
ssdov[k]='_';  k++; ssdov[k]='А'; k++; ssdov[k]='с';  k++; ssdov[k]='п';  k++; ssdov[k]='и';  k++; ssdov[k]='р';  k++;
ssdov[k]='а';  k++; ssdov[k]='н'; k++; ssdov[k]='т';  k++; ssdov[k]='у';  k++; ssdov[k]='р';  k++; ssdov[k]='а';  k++;
ssdov[k]='\\'; k++; ssdov[k]='\\';k++; ssdov[k]='t';  k++; ssdov[k]='m';  k++; ssdov[k]='p';  k++; ssdov[k]='\\'; k++;
ssdov[k]='\\'; k++; ssdov[k]='\0';
k=0;
szfovk[k]='V'; k++; szfovk[k]='y'; k++; szfovk[k]='v'; k++; szfovk[k]='o'; k++; szfovk[k]='d'; k++; szfovk[k]='v'; k++; 
szfovk[k]='F'; k++; szfovk[k]='i'; k++; szfovk[k]='l'; k++; szfovk[k]='e'; k++; szfovk[k]='-'; k++; szfovk[k]=pbf; k++; 
szfovk[k]='.'; k++;  szfovk[k]='t'; k++; szfovk[k]='x'; k++; szfovk[k]='t'; k++; szfovk[k]='\0';
skpovk=new char[dsov]; sdvovk=new char[dsov]; sppovk=new char[dsov]; 
skpov=new char[2*dsov]; sdvov=new char[2*dsov]; sppov=new char[2*dsov];
if ((!skpovk) || (!sdvovk) || (!sppovk) || (!skpov) || (!sppov) || (!sdvov)) { cout << snm << endl; k=getchar(); exit(1); } 
for (k=0; k<(2*dsov); k++) { skpov[k]='\0'; sdvov[k]='\0'; sppov[k]='\0'; szfov[k]='\0'; }
for (k=0; k<dsov; k++) { skpovk[k]='\0'; sdvovk[k]='\0'; sppovk[k]='\0'; }
k=0;
skpovk[k]='K'; k++; skpovk[k]='o'; k++; skpovk[k]='e'; k++; skpovk[k]='f'; k++; skpovk[k]='f'; k++; skpovk[k]='i'; k++;
skpovk[k]='c'; k++; skpovk[k]='i'; k++; skpovk[k]='e'; k++; skpovk[k]='n'; k++; skpovk[k]='t'; k++; skpovk[k]='_'; k++;
skpovk[k]='p'; k++; skpovk[k]='o'; k++; skpovk[k]='g'; k++; skpovk[k]='l'; k++; skpovk[k]='o'; k++; skpovk[k]='s'; k++;
skpovk[k]='c'; k++; skpovk[k]='h'; k++; skpovk[k]='e'; k++; skpovk[k]='n'; k++; skpovk[k]='i'; k++; skpovk[k]='y'; k++;
skpovk[k]='a'; k++; skpovk[k]='_'; k++; 
if (!vyve) { skpovk[k]='s'; k++; skpovk[k]='h'; k++; skpovk[k]='a'; k++; }
else if (vyve==1) { skpovk[k]='v'; k++; skpovk[k]='e'; k++; skpovk[k]='r'; k++; }
else if (vyve==2) { skpovk[k]='i'; k++; skpovk[k]='t'; k++; skpovk[k]='o'; k++; skpovk[k]='m'; k++; 
if (!vm)   { skpovk[k]='4'; k++; skpovk[k]='4'; k++; }
if (vm==1) { skpovk[k]='6'; k++; skpovk[k]='2'; k++; }
if (vm==2) { skpovk[k]='8'; k++; skpovk[k]='6'; k++; }
if (vm==2) { skpovk[k]='1'; k++; skpovk[k]='0'; k++; }
skpovk[k]='0'; k++;
}
else if (vyve==3) { skpovk[k]='k'; k++; skpovk[k]='v'; k++; skpovk[k]='i'; k++; 
	 if (vm==4) { skpovk[k]='4'; k++; }
else if (vm==5) { skpovk[k]='5';  k++; }
else if (vm==6) { skpovk[k]='6';  k++; }
else if (vm==7) { skpovk[k]='7';  k++; }
else if (vm==8) { skpovk[k]='8';  k++; }
else if (vm==9) { skpovk[k]='9';  k++; }
else if (vm==10) { skpovk[k]='1';  k++; skpovk[k]='0';  k++; }
else { cout << "Net takoy marki KVI!" << endl; k=getchar(); exit(1); }
skpovk[k]='0'; k++; skpovk[k]='0'; k++; }
skpovk[k]='.'; k++; skpovk[k]='t'; k++; skpovk[k]='x'; k++; skpovk[k]='t'; k++; skpovk[k]='\0';
k=0;
sdvovk[k]='D'; k++; sdvovk[k]='l'; k++; sdvovk[k]='i'; k++; sdvovk[k]='n'; k++; sdvovk[k]='a'; k++; sdvovk[k]='_'; k++; 
sdvovk[k]='v'; k++; sdvovk[k]='o'; k++; sdvovk[k]='l'; k++; sdvovk[k]='n'; k++; sdvovk[k]='y'; k++; sdvovk[k]='_'; k++; 
if (!vyve) { sdvovk[k]='s'; k++; sdvovk[k]='h'; k++; sdvovk[k]='a'; k++; }
else if (vyve==1) { sdvovk[k]='v'; k++; sdvovk[k]='e'; k++; sdvovk[k]='r'; k++; }
else if (vyve==2) { sdvovk[k]='i'; k++; sdvovk[k]='t'; k++; sdvovk[k]='o'; k++; sdvovk[k]='m'; k++; }
else if (vyve==3) { sdvovk[k]='k'; k++; sdvovk[k]='v'; k++; sdvovk[k]='i'; k++; }
sdvovk[k]='.'; k++; sdvovk[k]='t'; k++; sdvovk[k]='x'; k++; sdvovk[k]='t'; k++; sdvovk[k]='\0';
k=0;
sppovk[k]='P'; k++; sppovk[k]='o'; k++; sppovk[k]='k'; k++; sppovk[k]='a'; k++; sppovk[k]='z'; k++; sppovk[k]='a'; k++; 
sppovk[k]='t'; k++; sppovk[k]='_'; k++; sppovk[k]='p'; k++; sppovk[k]='r'; k++; sppovk[k]='e'; k++; sppovk[k]='l'; k++;
sppovk[k]='o'; k++; sppovk[k]='m'; k++; sppovk[k]='l'; k++; sppovk[k]='e'; k++; sppovk[k]='n'; k++; sppovk[k]='_'; k++;
if (!vyve) { sppovk[k]='s'; k++; sppovk[k]='h'; k++; sppovk[k]='a'; k++; }
else if (vyve==1) { sppovk[k]='v'; k++; sppovk[k]='e'; k++; sppovk[k]='r'; k++; }
else if (vyve==2) { sppovk[k]='i'; k++; sppovk[k]='t'; k++; sppovk[k]='o'; k++; sppovk[k]='m'; k++; }
else if (vyve==3) { sppovk[k]='k'; k++; sppovk[k]='v'; k++; sppovk[k]='i'; k++; }
sppovk[k]='.'; k++; sppovk[k]='t'; k++; sppovk[k]='x'; k++; sppovk[k]='t'; k++; sppovk[k]='\0'; //0 - шамот, 1 - вермикулит, 2 - ИТОМ, 3 - КВИ
strcpy(szfov, ssdov); strcat(szfov, szfovk); k=strlen(szfov)+1; szfov[k]='\0';
strcpy(sppov,ssdov); strcat(sppov,sppovk); k=strlen(sppov)+1; sppov[k]='\0';
strcpy(skpov,ssdov); strcat(skpov,skpovk); k=strlen(skpov)+1; skpov[k]='\0'; 
strcpy(sdvov,ssdov); strcat(sdvov,sdvovk); k=strlen(sdvov)+1; sdvov[k]='\0'; 
k=4; char **unau=new char*[k];
if (szfovk) delete []szfovk; if (skpovk) delete[]skpovk;
if (sdvovk) delete []sdvovk; if (sppovk) delete []sppovk; 
if (sndov) delete[]sndov; if (ssdov) delete[]ssdov;
k=0; unau[k]=szfov; k++; unau[k]=sppov; k++; unau[k]=skpov; k++; unau[k]=sdvov;
return unau;
}
void osvpamov(char *szfov, char *sppov, char *skpov, char *sdvov)
{ if (szfov) delete[]szfov; if (sppov) delete[]sppov; if (skpov) delete []skpov; if (sdvov) delete []sdvov; }
double *koefoslab(double wmg, double wsi, double wal, double *tere, int n, double *kuo)
{ int lt=n, k; 
double wo=wmg+wsi+wal, kmg=0.0, kal=0.0, ksi=0.0, ht=1e0, eps=1e-6;
double *mgo=new double[lt], *alo=new double[lt], *sio=new double[lt];
if ((!mgo) || (!alo) || (!sio)) { cout << "No memory" << endl; k=getchar(); exit(1); }
sio=oslaintestepcher(lt,0,sio);
alo=oslaintestepcher(lt,1,alo);
mgo=oslaintestepcher(lt,2,mgo); //for (k=0; k<n; k++) cout << "sio = " << sio[k] << "\talo = " << alo[k] << "\tmgo = " << mgo[k] << endl; //cout << "wmg = " << wmg << "\twsio = " << wsi << "\twal = " << wal << "\two = " << wo << endl;
for (k=0; k<n; k++) {
kmg=mgo[k]; kal=alo[k]; ksi=sio[k];
if (kmg<0.0) kmg=0.0; if (kmg>ht) kmg=ht; 
if (kal<0.0) kal=0.0; if (kal>ht) kal=ht; 
if (ksi<0.0) ksi=0.0; if (ksi>ht) ksi=ht;
if (fabs(wo)>eps) 
kuo[k]=(kmg*wmg+kal*wal+ksi*wsi)/wo; 
else kuo[k]=0.0; } //for (k=0; k<n; k++) cout << "kuo = " << kuo[k] << "\ttere = " << tere[k] << endl;
delete []mgo; delete []alo; delete []sio; return kuo; }
double epsisred(double T, double *tdkusctm, double *dkusctm, int dkusctl, double *dkoscet, double *dkoscem, int dkoscel, int vybves, int vybmar)
{ int k=1, j=0;
char **ms=napolStrok(), *snm=ms[j], *sfno=ms[k], *szfov=NULL, *sppov=NULL, *skpov=NULL, *sdvov=NULL; 
if (ms) delete[]ms; ms=napNazvFile(vybves, vybmar, snm);
k=0; szfov=ms[k]; k++; sppov=ms[k]; k++; skpov=ms[k]; k++; sdvov=ms[k]; if (ms) delete[]ms;
double *dv=dliny_voln(dv, snm, sdvov, sfno);
double *npp=new double[dlar], hf=1e0; if (!npp) { cout << snm << endl; k=getchar(); exit(1); }
for (k=0; k<dlar; k++) npp[k]=0.0; //cout << "ce = " << dlar << endl;
npp=Kramers_Kronig(dv, dlar, npp, sppov, sfno, snm, vybves, vybmar); //cout << "k = " << k << endl;
double dkusct=opredKTPTKTochSha(dkusctm, tdkusctm, T, dkusctl); 
if (dkusct>hf) dkusct=hf; if (dkusct<0.0) dkusct=0.0; 
double dkosce=opredKTPTKTochSha(dkoscem, dkoscet, T, dkoscel); 
if (dkosce>hf) dkosce=hf; if (dkosce<0.0) dkosce=0.0; 
double ume=dkusct*dkosce; 
double *epsil=new double[dlar];
if (!epsil) { cout << snm << endl; k=getchar(); exit(1); } 
epsil=epsilnu(npp, dlar, epsil); //cout << "ume = " << ume << "\t";
double epssr=usredVelichPlank(dv, epsil, npp, T, dlar, ume); //cout << "T = " << T << "\teps_s = " << epssr << "\tume = " << ume << "\n";
if (dv) delete []dv; if (npp) delete []npp; if (epsil) delete []epsil;
osvpamov(szfov, sppov, skpov, sdvov);
if (snm) delete[]snm; if (sfno) delete[]sfno;
return epssr; }
double reflVes(double T, int vybves, int vybmar)
{ int k=1, j=0;
char **ms=napolStrok(), *snm=ms[j], *sfno=ms[k], *szfov=NULL, *sppov=NULL, *skpov=NULL, *sdvov=NULL; 
if (ms) delete[]ms;
ms=napNazvFile(vybves, vybmar, snm);
k=0; szfov=ms[k]; k++; sppov=ms[k]; k++; skpov=ms[k]; k++; sdvov=ms[k]; if (ms) delete[]ms;
double *dv=NULL;
dv=dliny_voln(dv, 0);
double *npp=new double[dlar];
if (!npp) { cout << snm << endl; k=getchar(); exit(1); } 
for (k=0; k<dlar; k++) npp[k]=0.0;
npp=Kramers_Kronig(dv, dlar, npp, sppov, sfno, snm, vybves, vybmar);
double *ronu=ronusha(npp, dlar), t=0.0, ume=1e0;
t=usredVelichPlank(dv, ronu, npp, T, dlar, ume);
if (dv) delete []dv; if (npp) delete []npp; if (ronu) delete []ronu; 
osvpamov(szfov, sppov, skpov, sdvov);
return t; }
double *Kramers_Kronig(double *dv, int lear, double *nnu, char *sppov, char *sfno, char *snm, int vv, int vm)
{	int k=0; ofstream fo;
fo.open(sppov, ios_base::out | ios_base::trunc); 
if (!fo.is_open()) { cout << sfno << endl; k=getchar(); exit(1); }
	double c0=299792458.0, *nu=new double[lear], *nus=new double[lear], *nnus=new double[lear], c=0.0;
	double *alsr=new double[lear], e=1e-20, ht=1e0, ko=1e-2; //for (k=0; k<lear; k++) alsr[k]=0.0;
	if ((!nu) || (!nnus) || (!nus) || (!alsr)) { cout << snm << endl; k=getchar(); exit(1); } 
	alsr=Koef_Pogl(alsr); 
	for (k=0; k<lear; k++) if (dv[k]>e) nu[k]=2e0*pi*c0/dv[k]; else nu[k]=0.0; //for (k=0; k<lear; k++) cout << "k = " << k << "\tnu = " << nu[k] << endl; cout << "lear = " << lear << endl; k=getchar();
nus=izmMasChast(nu, enu, lear, nus); //for (k=0; k<lear; k++) cout << "k = " << k << "\tnu = " << nus[k] << endl; cout << "lear = " << lear << endl; k=getchar();
nnus=PokazPrelomAl(nu, nus, alsr, c0, lear, nnus); //for (k=0; k<lear; k++) cout << "k = " << k << "\tnu = " << nnus[k] << endl; cout << "lear = " << lear << endl; k=getchar();
nnu=PokazPrelomPribl(nu, nus, nnus, lear, nnu, enu); //for (k=0; k<lear; k++) cout << "k = " << k << "\tnu = " << nnu[k] << endl; cout << "lear = " << lear << endl; k=getchar();
if (nu) delete[]nu; if (nnus) delete[]nnus; if (nus) delete[]nus; if (alsr) delete[]alsr;
if (!vv) c=1.56-ht;
else if (vv==1) c=1.543-ht;
else if (vv==2) { 
double w1=0.0, w2=0.0, rosha=219.0*1e1, rover=25.0*1e1, phi1=0.0, phi2=0.0, v1=0.0, v2=0.0;
if (!vm) { w1=6e1*ko; w2=ht-w1; }
else if (vm==1) { w1=5e1*ko; w2=ht-w1; }
else if (vm==2) { w1=25.0*ko; w2=ht-w1; }
else if (vm==3) { w1=2e1*ko; w2=ht-w1; }
else { cout << "Oshibka v vybore marki ITOM!"; k=getchar(); exit(1); }
v1=w1/rover; v2=w2/rosha; 
phi1=v1/(v1+v2); phi2=v2/(v1+v2); 
c=(1.543-ht)*phi1+(1.56-ht)*phi2; }
else if (vv==3) { double wv=0.0; 
if (vm==4) wv=49e-2; else if (vm==5) wv=37e-2; else if (vm==6) wv=31e-2; else if (vm==7) wv=27e-2;
else if (vm==8) wv=23e-2; else if (vm==9) wv=22e-2; else if (vm==10) wv=2e-1; 
else { cout << "Oshibka v vybore marki KVI!"; k=getchar(); exit(1); }
c=(1.543-ht)*wv+(ht-wv)*(1.56-ht); }
else { cout << "Oshibka v vybore veschestva!"; k=getchar(); exit(1); }
for (k=0; k<lear; k++) { nnu[k]=nnu[k]+c; fo << nnu[k] << endl; } fo.close(); 
return nnu; }
double *ronu(double *npp, int lear, char *snm)
{ int k=0; double *ronu=new double[lear], n=0.0, hf=1e0; if (!ronu) { cout << snm << endl; k=getchar(); exit(1); } 
for (k=0; k<lear; k++) {
	n=fabs(npp[k]);
    ronu[k]=0.5+((n-hf)*(3.0*n+hf))/(6.0*pow(n+hf,2.0))-(2.0*pow(n,3.0)*(pow(n,2.0)+2.0*n-hf))/((pow(n,2.0)+hf)*(pow(n,4.0)-hf)); 
    ronu[k]=ronu[k]+(8.0*pow(n,4.0)*(pow(n,4.0)+hf)*log(n))/((pow(n,2.0)+hf)*pow((pow(n,4.0)-hf),2.0));
    ronu[k]=ronu[k]+pow(n,2.0)*pow((pow(n,2.0)-hf),2.0)*log(fabs((n-hf)/(n+hf)))/pow((pow(n,2.0)+hf),3.0); }
return ronu; }
double *epsilnu(double *npp, int lear, double *epsarr)
{ int k=0; double eps=0.0, n=0.0, hf=1e0;
for (k=0; k<lear; k++) {
	n=fabs(npp[k]);
	eps=(4.0*n+2.0)/3.0/pow((n+hf),2.0);
	eps=eps+2.0*pow(n,3.0)*(pow(n,2.0)+2.0*n-hf)/(pow(n,2.0)+hf)/(pow(n,4.0)-hf);
	eps=eps-8.0*pow(n,4.0)*(pow(n,4.0)+hf)*log(n)/(pow(n,2.0)+hf)/pow((pow(n,4.0)-hf),2.0);
	eps=eps-pow(n,2.0)*log((n-hf)/(n+hf))*pow((pow(n,2.0)-hf),2.0)/pow((pow(n,2.0)+hf),3.0); 
	epsarr[k]=eps; }
return epsarr; }
double *PokazPrelomAl(double *nu, double *nus, double *ka, double vl, int arle, double *np, char *snm)
{ int k, j, p=arle-1, q; double dkpdo, podln, *fn=new double[arle], hf=1e0; 
if (!fn) { cout << snm << endl; k=getchar(); exit(1); } 
for (k=0; k<arle; k++) {
    q=0; for (j=0; j<p; j++) {
        dkpdo=(ka[j+1]-ka[j])/(nu[j+1]-nu[j]);
        podln=((nu[j]+nus[k])/(nu[j]-nus[k]));
        podln=fabs(podln);
        fn[q]=dkpdo*log(podln);
        q++; }
    fn[p]=fn[p-1];
    np[k]=hf+(vl/pi)*trapz(nu,fn,arle)/2.0/nu[k]; }
delete []fn; return np; }
double *izmMasChast(double *nu, double epnu, int chel, double *nus)
{ int k, p=chel-1;
for (k=0; k<p; k++) nus[k]=epnu*(nu[k+1]-nu[k])+nu[k];
nus[p]=epnu*(nu[p-1]-nu[p-2])+nu[p-1];
return nus; }
double *PokazPrelomPribl(double *x, double *xs, double *ys, int arle, double *y, double e)
{ int k, p=arle-1; double ko=0.0;
for (k=0; k<arle; k++) y[k]=0.0;
for (k=0; k<(p-1); k++) {
	ko=(ys[k+1]-ys[k])/(xs[k+1]-xs[k]); 
	y[k]=ko*(x[k]-xs[k])+ys[k]; }
k=p-1;
ko=(ys[k]-ys[k-1])/(xs[k]-xs[k-1]);
y[k]=(x[k]-xs[k-1])*ko+ys[k];
y[p]=(ys[p]-ys[k])*e+ys[k];
return y; }
double *Koef_Pogl(double *kp, char *snm, char *skpov, char *sfno)
{	ifstream fin; double p; int k=0, q=100; char *s=new char[q]; 
if (!s) { cout << snm << endl; k=getchar(); exit(1); }
fin.open(skpov, ios_base::in); for (k=0; k<q; k++) s[k]='\0'; //cout << "k = " << k << endl;
if (!fin.is_open()) { cout << sfno << endl; k=getchar(); exit(1); }
k=0; while ((!fin.eof()) && (k<dlar)) { fin.getline(s,q,'\n'); p=atof(s); kp[k]=p; k++; }
fin.close(); if (s) delete []s; 
return kp; }
double *dliny_voln(double *dv, char *snm, char *sdvov, char *sfno)
{ ifstream fin; double p, ht=1e-2; int k, q=100, lear; char *s=new char[q]; 
if (!s) { cout << snm << endl; k=getchar(); exit(1); } for (k=0; k<q; k++) s[k]='\0';
fin.open(sdvov, ios_base::in); 
if (!fin.is_open()) { cout << sfno << endl; k=getchar(); exit(1); }
k=0; while (!fin.eof()) { fin.getline(s,q,'\n'); p=atof(s); k++; } 
fin.clear(); fin.seekg(0);
lear=k; dlar=lear; dv=new double[dlar]; if (!dv) { cout << snm << endl; k=getchar(); exit(1); } 
for (k=0; k<dlar; k++) dv[k]=0.0;
k=0; while ((!fin.eof()) && (k<dlar)) { fin.getline(s,q,'\n'); p=atof(s); dv[k]=p;  k++; }
fin.close();
for (k=0; k<lear; k++) dv[k]=ht/dv[k];
if (s) delete []s; 
return dv; } 
void zapisvfile(double *vyvod, int dlina, char *nazf, char *sfno)
{ 	if (dlina>0) { int k; FILE *fo = fopen(nazf, "a+"); 
	if (!fo) { cout << sfno << endl; k=getchar(); exit(1); } cout << "dl_zap" << dlina << endl;
	for (k=0; k<dlina; k++) fprintf(fo,"%0.15le\n",vyvod[k]); //fprintf(fo,"%\n\n\n"); 
	fclose(fo); } }
double trapz(double *arx, double *ary, int lenarr)
{ int k; double s=0.0;
for (k=1; k<lenarr; k++) s=s+(ary[k]+ary[k-1])*(arx[k]-arx[k-1])/2e0;
return s; }
double *oslaintestepcher(int dm, int vy, double *sc) //рассчитывает ослабление интегральной степени черноты для SiO2, Al2O3, MgO
{	int n=2, k, j, w=0, p=6, q=10, m, l; double **scv=new double*[p], **hsv=new double*[p];
	double *hs=new double[n], **kor=new double*[q], t0, t1, *vscs=new double[dm], *vscm=new double[dm], *vsca=new double[dm], *po=NULL;
if ((!scv) || (!sc) || (!hsv) || (!hs) || (!kor)) { cout << "No memory" << endl; k=getchar(); exit(1);}
	k=0;
	sc[k]=1e0;   k++; sc[k]=0.94;  k++; sc[k]=0.87;  k++; sc[k]=0.801; k++; sc[k]=0.736; k++; 
	sc[k]=0.676; k++; sc[k]=0.635; k++; sc[k]=0.590; k++; sc[k]=0.567; k++; sc[k]=0.543; k++; 
	sc[k]=0.53;  k++; sc[k]=0.525; k++; sc[k]=0.515; k++; sc[k]=0.507; //степень черноты магнезита
	w=0; scv[w]=sc; sc=new double[dm];
	k=0;
	hs[k]=98e-2; k++; hs[k]=2e-2;
	k=0;
	t0=hs[k]; k++; t1=hs[k];
	k=0;
	hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1);
	hsv[w]=hs; w++; hs=new double[n]; //Магнезит: 98 % - MgO, 2 % - SiO2
	if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1); }
k=0;
sc[k]=1e0;   k++; sc[k]=0.976; k++; sc[k]=0.949; k++; sc[k]=0.905; k++; sc[k]=0.859; k++; 
sc[k]=0.812; k++; sc[k]=0.774; k++; sc[k]=0.737; k++; sc[k]=0.709; k++; sc[k]=0.681; k++; 
sc[k]=0.661; k++; sc[k]=0.639; k++; sc[k]=0.626; k++; sc[k]=0.620; //степень черноты шамота
scv[w]=sc; sc=new double[dm]; 
k=0; 
hs[k]=56e-2; k++; hs[k]=396e-3; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Шамот: 56 % - SiO2, 39,6 % - Al2O3 - шамот, считаем, что магния нет в шамоте
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1);}
k=0;
sc[k]=1e0;    k++; sc[k]=98e-2;  k++; sc[k]=951e-3; k++; sc[k]=92e-2;  k++; sc[k]=883e-3; k++; 
sc[k]=853e-3; k++; sc[k]=821e-3; k++; sc[k]=79e-2;  k++; sc[k]=767e-3; k++; sc[k]=746e-3; k++;
sc[k]=73e-2;  k++; sc[k]=715e-3; k++; sc[k]=705e-3; k++; sc[k]=692e-3; //степень черноты корундошамота
scv[w]=sc; sc=new double[dm]; 
k=0; 
hs[k]=28e-2; k++; hs[k]=7e-1; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Корундошамот: 28 % - SiO2, 70 % - Al2O3 - корундошамот
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1);}
k=0;
sc[k]=1.0000; k++; sc[k]=983e-3; k++; sc[k]=936e-3; k++; sc[k]=867e-3; k++; sc[k]=819e-3; k++; 
sc[k]=721e-3; k++; sc[k]=659e-3; k++; sc[k]=593e-3; k++; sc[k]=541e-3; k++; sc[k]=49e-2;  k++; 
sc[k]=453e-3; k++; sc[k]=429e-3; k++; sc[k]=403e-3; k++; sc[k]=384e-3; //степень черноты каолинового теплоизоляционного кирпича (КТК)
scv[w]=sc; sc=new double[dm];
k=0; 
hs[k]=57e-2; k++; hs[k]=4e-1; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //КТК: 57 % - SiO2, 40 % - Al2O3
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1); }
k=0;
sc[k]=1.000;  k++; sc[k]=984e-3; k++; sc[k]=941e-3; k++; sc[k]=882e-3; k++; sc[k]=813e-3; k++; 
sc[k]=751e-3; k++; sc[k]=695e-3; k++; sc[k]=641e-3; k++; sc[k]=594e-3; k++; sc[k]=558e-3; k++; 
sc[k]=53e-2;  k++; sc[k]=499e-3; k++; sc[k]=479e-3; k++; sc[k]=462e-3; //степень черноты муллита
scv[w]=sc; sc=new double[dm];
k=0; 
hs[k]=28e-2; k++; hs[k]=72e-2; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Муллит: 28 % - SiO2, 72 % - Al2O3
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1);}
k=0;
sc[k]=1e0;    k++; sc[k]=984e-3; k++; sc[k]=953e-3; k++; sc[k]=917e-3; k++; sc[k]=854e-3; k++; 
sc[k]=808e-3; k++; sc[k]=756e-3; k++; sc[k]=711e-3; k++; sc[k]=578e-3; k++; sc[k]=523e-3; k++; 
sc[k]=495e-3; k++; sc[k]=468e-3; k++; sc[k]=448e-3; k++; sc[k]=429e-3; //степень черноты кремнезема
scv[w]=sc; 
k=0; 
hs[k]=985e-3; k++; hs[k]=1e-2; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Кремнезем: 98,5 % - SiO2, 1 % - Al2O3
hsv[w]=hs; w++;
for (j=0; j<dm; j++) { k=0; 
for (m=1; m<p; m++)
	for (l=m+1; l<p; l++) {
		kor[k]=RasKorOV(hsv[m][0],hsv[m][1],hsv[l][0],hsv[l][1],scv[m][j],scv[l][j]); k++; }
		vscs[j]=UsredMasOV(kor,q,0); vsca[j]=UsredMasOV(kor,q,1); } 
if (vy==2) { hs=hsv[0]; sc=scv[0]; 
for (k=0; k<dm; k++) vscm[k]=(sc[k]-vscs[k]*hs[1])/hs[0]; } 
for (k=0; k<q; k++) { sc=kor[k]; delete []sc; }
for (k=0; k<p; k++) { hs=hsv[k]; delete []hs; } delete []hsv; //0 - SiO2, 1 - Al2O3, 2 - MgO
for (k=0; k<p; k++) { po=scv[k]; delete []po; } delete []scv;
if (!vy) { delete []vsca; delete []vscm; po=vscs; }
if (vy==1) { delete []vscm; delete []vscs; po=vsca; }
if (vy==2) { delete []vsca; delete []vscs; po=vscm; } 
return po;
}
void vyvomatr(double **matr, int m, int n)
{ int k, j; for (k=0; k<m; k++) { for (j=0; j<n; j++) {	cout << matr[k][j] << "\t"; } cout << endl; } }
double *RasKorOV(double a11, double a12, double a21, double a22, double b1, double b2)
{ int l=2, k; double **mas=new double*[l], *st=new double[l], *bs=new double[l], *x, *kor=new double[l];
if ((!mas) || (!st) || (!bs)) { cout << snmov << endl; k=getchar(); exit(1);}
st[0]=a11; st[1]=a12; mas[0]=st; st=new double[l];
if (st) { st[0]=a21; st[1]=a22; mas[1]=st; } else { cout << snmov << endl; k=getchar(); exit(1);}
bs[0]=b1; bs[1]=b2;
x=reshMetKram(mas,bs,l);
if ((x[0]>=0.0) && (x[1]>=0.0) && (x[0]<=1.0) && (x[1]<=1.0)) {
	kor[0]=x[0]; kor[1]=x[1]; }
else { kor[0]=0.0; kor[1]=0.0; }
for (k=0; k<l; k++) { st=mas[k]; delete []st; } delete []bs;
return kor; }
double UsredMasOV(double **kor, int q, int vy)
{ double s1=0.0, s2=0.0, p=0.0, *s; int j;
for (j=0; j<q; j++) { s=kor[j];
    if ((s[0]>0.0) && (s[1]>0.0)) {
        s1=s1+s[0]; s2=s2+s[1]; p=p+1e0; } }
if (!vy) return s1/p; else return s2/p; }
double KoefPoglRossel(double tem, double *dl, double *alfs, double *npp, int n)
{ double PP=6.6260755e-34, PB=1.380658e-23, c0=299792458.0, C1=PP*pow(c0,2.0), C2=PP*c0/PB;
double sig=2e0*C1*pow(pi,5e0)/15.0/pow(C2,4e0), np2=nsredPlank2(dl,npp,tem,n), t, chi, zna, chasc, chasz, rs;
double *Ibc=new double[n], *Ibz=new double[n], *dv=new double[n], dlv, c1m, c2m; int k;
if ((!Ibc) || (!Ibz) || (!dv)) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<n; k++) {
c1m=C1/pow(npp[k],2e0);
c2m=C2/npp[k];
t=(pi/2.0)*(c1m*c2m)/sig;
dlv=dl[k]/npp[k];
chi=exp(c2m/tem/dlv);
zna=pow((chi-1.0),2.0);
Ibz[k]=t*chi/zna/pow(dlv,6e0)/pow(tem,5e0)/np2; 
Ibc[k]=Ibz[k]/alfs[k];
dv[k]=dlv; }
chasc=trapz(dv,Ibc,n);
chasz=trapz(dv,Ibz,n); 
rs=chasc/chasz;
rs=1e0/rs;
delete []Ibc; delete []Ibz; delete []dv;
return rs; }
double nsredPlank2(double *dv, double *npp, double tem, int n)
{ int k; double *npp2=new double[n], t, ume=1e0; 
if (!npp2) {cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) npp2[k]=pow(npp[k],2e0); 
t=usredVelichPlank(dv, npp2, npp, tem, n, ume);
delete []npp2; return t; }
double usredVelichPlank(double *dv, double *uv, double *npp, double tem, int n, double umensh)
{ double PP=6.6260755e-34, PB=1.380658e-23, c0=299792458.0, vl, c1, c2, lambda; 
int k; double *Ib, *Ibn, *dl, nc, nz, e=1e-20;
Ib=new double[n]; Ibn=new double[n]; dl=new double[n];
if ((!Ib) || (!Ibn) || (!dl)) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<n; k++) {
    vl=c0/npp[k];
    c1=PP*pow(vl,2e0);
    c2=PP*vl/PB;
    lambda=dv[k]/npp[k];
    Ib[k]=2e0*pi*c1/(pow(lambda,5e0)*(exp(c2/lambda/tem)-1e0));
Ibn[k]=uv[k]*Ib[k];
dl[k]=lambda; }
nc=trapz(dl,Ibn,n);
nz=trapz(dl,Ib,n);
if (fabs(nz)>e) nc=umensh*nc/nz; else nc=0.0;
delete []Ibn; delete []Ib; delete []dl;
return nc; }
double *KoefPoglRosselNac(double *tem, int ident, int arle, double wmgo, double wsio, double walo, double *dkoscet, double *dkoscem, int dkoscel, double dkosups, double dkokpiko, double *tdkusctm, double *dkusctm, int dkusctl, int vi, int vyve)
{	int k=0, j=0, dm=0, vm=vi; 
	napstrdir(vyve, vm); //cout << "vyve = " << vyve << "vm = " << vm << endl; 
	double *dl=NULL, *alfs=NULL, *npp=NULL, *kpr=new double[arle], krpk=0.0, temk=0.0, hf=1e0;
	double *ktr=new double[arle], dkusct=1e0, dkosce=1e0; //КТП по Росселанду
	double dko2=dkosups, dko4=dkokpiko, dko0=1e0; //учет пористой структуры и соотношения КП/КО
	dl=dliny_voln(dl, snm, sdvov, sfno);
	alfs=new double[dlar]; if (!alfs) { cout << snm << endl; k=getchar(); exit(1); } for (k=0; k<dlar; k++) alfs[k]=0.0;
	alfs=Koef_Pogl(alfs); npp=new double[dlar]; if (!npp) { cout << snm << endl; k=getchar(); exit(1); } 
	for (k=0; k<dlar; k++) npp[k]=0.0; dm=dlar; 
	npp=Kramers_Kronig(dl, dlar, npp, sppov, sfno, snm, vybves, vybmar);
	for (k=0; k<dm; k++) alfs[k]=alfs[k]*dko2*dko4; //cout << "dlark = " << dlark << endl; //for (k=0; k<(dm/50); k++) cout << alfs[k] << "\t"; for (k=0; k<(dm/50); k++) cout << npp[k] << "\t"; 
	double *GnpT=new double[arle], sigma=5.67e-8, *npT=new double[arle]; //производная dn/dT, функция <n>(T)
	if ((!kpr) || (!GnpT) || (!ktr) || (!npT)) { cout << snm << endl; k=getchar(); exit(1); }
	for (k=0; k<arle; k++) {
	dkusct=1e0; dkosce=1e0; dko0=dkusct*dkosce; temk=tem[k];
	npT[k]=usredVelichPlank(dl, npp, npp, temk, dm, dko0); 
	cout << "\tal_sred = " << usredVelichPlank(dl, alfs, npp, temk, dm, dko0); //dkusct=opredKTPTKTochSha(dkusctm, tdkusctm, temk, dkusctl); dkosce=opredKTPTKTochSha(dkoscem, dkoscet, temk, dkoscel);
	if ((dkosce<=0.0) || (dkosce>1e0)) dkosce=1e0; if ((dkusct<=0.0) || (dkusct>1e0)) dkusct=1e0; 
	for (j=0; j<dm; j++) alfs[j]=alfs[j]*dko0;
	krpk=KoefPoglRossel(temk, dl, alfs, npp, dm); kpr[k]=krpk; 
	for (j=0; j<dm; j++) alfs[j]=alfs[j]/dko0; 
	cout << "\ttemp = " << temk << "\tKoef_pogl_Ross = " << krpk << "\tdko2 = " << dko2 << "\tdko4 = " << dko4 << endl; 
	} for (k=1; k<arle; k++) GnpT[k-1]=(npT[k]-npT[k-1])/(tem[k]-tem[k-1]); k=arle-1; GnpT[k]=2e0*GnpT[k-1]-GnpT[k-2];
	for (k=0; k<arle; k++) {
	temk=tem[k]; krpk=kpr[k];
    krpk=fabs(8e0*npT[k]*sigma*pow(temk,3e0)/(3e0*krpk));
    krpk=krpk*fabs(2e0*npT[k]+temk*GnpT[k]);
	ktr[k]=krpk; }
	delete []dl; delete []npp; delete []alfs; delete []kpr; delete []GnpT; delete []npT;
	osvpamov(vyve); return ktr; }
double LuchKTPChudnovsky(double *Ab, double tem, int kost, double razm)
{ double s=0.0, t=0.0; int k; for (k=0; k<kost; k++) { s=s+Ab[k]; t=t+1e0; } s=s/t;
return ((3e0/4e0)*2e0*(5.67e-8)*pow(s,2e0)*pow(tem,3e0)*razm); }
double PoiskReflPhi(double phi, double phis, double n1, double n2)
{ double rpa=(n2*cos(phi)-n1*cos(phis))/(n2*cos(phi)+n1*cos(phis)); //коэффициент отражения параллельный плоскости падения
double rpe=(n1*cos(phi)-n2*cos(phis))/(n1*cos(phi)+n2*cos(phis)); //коэффициент отражения перпендикулярный к плоскости падения
return (fabs(pow(rpa,2e0))+fabs(pow(rpe,2e0)))/2e0; }
double ReflSred(double T)
{ int ide=0, k=1; 
napstrdir(k,k);
double *dv=NULL; dv=dliny_voln(dv, snm, sdvov, sfno);
double *npp=new double[dlar]; 
if (!npp) { cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<dlar; k++) npp[k]=0.0;
npp=Kramers_Kronig(dv, dlar, npp, sppov, sfno, snm, vybves, vybmar);
double *ronu=PoiSpeKoeOtr(npp,dlar), ume=1e0;
double t=usredVelichPlank(dv, ronu, npp, T, dlar, ume);
delete []dv; delete []npp; delete []ronu; k=1; osvpamov(k);
return t; }
double *PoiSpeKoeOtr(double *npp, int lear)
{ double maxphi=pi/2e0, minphi=0.0, phi, hphi=enu, n2, n1=1e0, *rs=new double[lear], p, r, psi; int k;
for (k=0; k<lear; k++) { phi=minphi; n2=npp[k]; p=0.0; r=0.0;
while (phi<maxphi) { psi=PoiskTetaShtrNach(n2-n1, phi, n1, n2); 
phi=phi+hphi; p=p+1e0; r=r+PoiskReflPhi(phi, psi, n1, n2); }
rs[k]=r/p; } return rs; }
double PoiskTetaShtrNach(double dnr, double phipad, double n1, double n2)
{ double phipre=0.0, a=0.0, b=pi/2e0, c, ep=tocrasov, ra=fabs(a-b);
double fa, fb, fc, phiprel;
int Nit=10000, k=0;
if (dnr<0.0) { phipre=asin(n2/n1); phipre=provUgla(phipre); }
if (((phipad<phipre) && (dnr<0.0)) || (dnr>0.0)) {
while ((ra>ep) && (k<Nit)) {
    c=(a+b)/2e0;
    fa=PraCha10(phipad,a,n1,n2);
    fb=PraCha10(phipad,b,n1,n2);
    fc=PraCha10(phipad,c,n1,n2);
    if ((fc*fb>0.0) && (fa*fc<0.0)) b=c; 
	if ((fc*fa>0.0) && (fb*fc<0.0)) a=c; 
	k++; ra=fabs(a-b); }
	phiprel=provUgla(c); }
else phiprel=pi/2e0;
return phiprel; }
double PraCha10(double phi, double phis, double n1, double n2)
{ return n1*sin(phi)-n2*sin(phis); }
double provUgla(double ugol)
{ double ugo=ugol;
if (ugo>pi/2e0) ugo=pi-ugo;
if (ugo<0.0) ugo=fabs(ugo);
return ugo; }
double SeryeStenkiRasIzl(double ww, double hh, double ll, double *T, double *EPS, double *HOs, double *RHOs, double *A, int *ns, int chst) //N серых диффузно излучающих поверхностей с диффузным и зеркальным коэффициентами отражения
{ int N=chst, fl=0; double b=0.0, t=0.0, epsi=1e-10; //t=T[ns[0]]; b=T[ns[1]]; if (fabs(t-b)<epsi) return 0.0;
double **Fs=new double*[N], z=0.0, r=0.0, y=0.0; int i=0, j=0, f=3, k=0; 
double *PIN=new double[N], *POUT=new double[N], *arg=new double[f], *p=NULL;
for (i=0; i<N; i++) { p=new double[N]; for (j=0; j<N; j++) p[j]=0.0; Fs[i]=p; }
	double sumq=0.0, sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	int *id=new int[N], iclsd=0; for (i=0; i<N; i++) id[i]=1; // Поверхности 1 и 3 (bottom and top)
    	for (i=0; i<N; i++) if (!(id[i])) PIN[i]=q[i]; else PIN[i]=sigma*pow(T[i], 4e0); //Заполнение массива PIN ПТП q и температурой T - перевод температур в ППИ, Fs(1,2)=F1-2, 
	iclsd=1; //угловые коэффициенты; конфигурация замкнута (iclsd=1), диагональные элементы не нужны //for (j=0; j<N; j++) cout << "RHOs ( " << j << " ) = " << RHOs[j] << endl;  for (j=0; j<N; j++) cout << "EPS ( " << j << " ) = " << EPS[j] << endl;
	k=0; arg[k]=fabs(hh); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(ww); //h,l,w
	k=39; i=0; j=1; Fs[i][j]=VIEW(k, f, arg, fl); //cout << "Fs(1,2) = " << Fs[0][1] << "\th = " << hh << "\tl = " << ll << "\tw = " << ww << "\tf = " << f << endl; //Fs(1,2)=F12
    k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(hh);     // a,b,c
	k=38; r=VIEW(k, f, arg, fl); 
	z=0.0; t=parlplates(ww, ww, 2.0*ww, ll, z, ll, hh); k=1; b=RHOs[k];
    i=0; j=2; Fs[i][j]=r+b*t; //Fs13=F13+rhos2*F1(2)3	
	t=perpplates(ww, 2.0*ww, z, hh, ll, z, ll); k=1; b=RHOs[k];
    i=0; j=3; k=1; Fs[i][j]=Fs[i][k]+b*t; //cout << "t = " << t << "\tb = " << b << endl; //Fs14=F14+rhos2*F1(2)4
    i=0; j=4; k=1; Fs[i][j]=0.5*(1.0-(1.0-RHOs[k])*Fs[i][k]-Fs[i][2]-Fs[i][3]); //смещены два прямоугольника на ww по оси X
    i=0; j=5; k=4; Fs[i][j]=Fs[i][k]; 
    k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(hh); //h,l,w
    k=39; r=VIEW(k, f, arg, fl); 
	k=0; b=RHOs[k];
	t=perpplates(hh, 2.0*hh, z, ww, ll, z, ll);
	i=1; j=2; Fs[i][j]=r+b*t; //Fs23=F23+rhos1*F1(1)3 //смещены два прямоугольника на hh по оси X, т.к. у них разные размеры
	k=0; arg[k]=fabs(ll); k++; arg[k]=fabs(hh); k++; arg[k]=fabs(ww); //a,b,c
    k=38; r=VIEW(k, f, arg, fl);
	k=0; b=RHOs[k];
	t=parlplates(hh, hh, 2.0*hh, ll, z, ll, ww);
	i=1; j=3; Fs[i][j]=r+b*t;  //Fs24=F24+rhos1*F2(1)4
    i=1; j=4; k=0; r=A[k]*Fs[k][i]/A[i];
	Fs[i][j]=0.5*(1.0-(1.0-RHOs[k])*r-Fs[i][2]-Fs[i][3]); //Fs25
    i=1; j=5; k=4; Fs[i][j]=Fs[i][k]; //Fs26
    k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); k++; arg[k]=2.0*fabs(hh);  // a,b,c
	k=38; r=VIEW(k, f, arg, fl);
	k=1; b=RHOs[k];
	t=parlplates(ww, ww, 2.0*ww, ll, z, ll, 2.0*hh);
    i=2; j=0; Fs[i][i]=RHOs[j]*(r+b*t); //Fs33=rhos1*F3(1)3+rhos1*rhos2*F3(12+21)3 //излучение стенки саму на себя
	r=perpplates(z, ww, hh, 2.0*hh, ll, z, ll);
	k=1; b=RHOs[k];
	t=perpplates(ww, 2.0*ww, hh, 2.0*hh, ll, z, ll);
    i=2; j=3; k=0; Fs[i][j]=Fs[k][j]+RHOs[k]*(r+b*t); // Fs34=Fs1-4+rhos1*F3(1)-3+rhos1*rhos2*F3(12+21)-4
    i=2; j=4; k=0; r=A[k]*Fs[k][i]/A[i]; t=A[1]*Fs[1][i]/A[i];
	Fs[i][j]=0.5*(1.0-(1.0-RHOs[k])*r-(1.0-RHOs[1])*t-Fs[i][i]-Fs[i][3]);
    i=2; j=5; k=4; Fs[i][j]=Fs[i][k];
    k=0; arg[k]=fabs(hh); k++; arg[k]=fabs(ll); k++; arg[k]=2.0*fabs(ww); //a,b,c
	k=38; r=VIEW(k, f, arg, fl);
	k=0; b=RHOs[k];
	t=parlplates(hh, hh, 2.0*hh, ll, z, ll, 2.0*ww);
    i=3; j=1; Fs[i][i]=RHOs[j]*(r+b*t); //Fs44=rhos2*F4(2)-4+rhos1*rhos2*F4(12+21)-4
    i=3; j=4; k=0; r=A[k]*Fs[k][i]; t=A[1]*Fs[1][i]; b=A[2]*Fs[2][i];
	Fs[i][j]=0.5*(1.0-((1.0-RHOs[k])*r+(1.0-RHOs[1])*t+b)/A[i]-Fs[i][i]);
    i=3; j=5; k=4; Fs[i][j]=Fs[i][k];
    i=4; j=5; k=0; r=A[k]*Fs[k][i]; t=A[1]*Fs[1][4]; b=A[2]*Fs[2][4]+A[3]*Fs[3][4];
	Fs[i][j]=0.5*(1.0-((1.0-RHOs[k])*r+(1.0-RHOs[1])*t+b)/A[4]); //Fs56 //for (j=0; j<N; j++) cout << "A ( " << j << " ) = " << A[j] << endl; //for (i=0; i<N; i++) { y=RHOs[i]+EPS[i]; cout << "Edinitsa = " << y << endl; } //for (i=0; i<N; i++) cout << "Temp ( " << i << " ) = " << T[i] << endl; //POUT=GRAYDIFFSPEC(iclsd, A, EPS, RHOs, HOs, Fs, id, PIN, POUT, N); //Решение системы уравнений методом Гаусса //for (i=0; i<N; i++) cout << POUT[i] << endl; 
	POUT=GRAYDIFFSPEC(iclsd, A, EPS, RHOs, HOs, Fs, id, PIN, POUT, N);
	sumq=0.0; //проверка полной ППИ
    	for (i=0; i<N; i++) //Вывод - перевод в температуры
		{ if (!(id[i])) T[i]=pow((POUT[i]/sigma), 1e0/4e0); else q[i]=POUT[i];
        	QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; 
			cout << "Surface " << i << "\tT [K] = " << T[i] << "\tq [W/m2] = " << q[i] << "\tQA [W] = " << QA[i]/ll << endl;
		} r=0.0; z=1.0; for (i=0; i<N; i++) if (ns[i]==i) { r=r+z*q[i]; z=-z; } 
		for (i=0; i<N; i++) { p=Fs[i]; delete []p; } delete []Fs; delete []PIN; 
		delete []POUT; delete []arg; delete []q; delete []QA; delete []id; 
		return r; }
double *GRAYDIFF(int iclsd, double *A, double *EPS, double *HOs, double **Fs, int *ID, double *PIN, double *POUT, int chst)
{
int N=chst, i, j;
double **qm=new double*[N], **em=new double*[N], **pm=new double*[N], *B=new double[N], *p, *idf=new double[N], hf=1e0, s, ikr, epsi=1e-10;
for (i=0; i<N; i++) { p=new double[N]; qm[i]=p; p=new double[N]; em[i]=p; p=new double[N]; pm[i]=p; }
for (i=0; i<N; i++) { s=0.0; for (j=0; j<ID[i]; j++) s=s+hf; idf[i]=s; }
for (i=1; i<N; i++)
	 for (j=0; j<i; j++)
		  if (fabs(A[i])>epsi) 
			  Fs[i][j]=(A[j]/A[i])*Fs[j][i]; else Fs[i][j]=0.0; //Для закрытой конфигурации, необходимо рассчитать диагональные элементы из правила суммирования
	if (iclsd==1) { for (i=0; i<N; i++) {
		Fs[i][i]=1e0;
		for (j=0; j<N; j++) {
			if (j!=i) {
				Fs[i][i]=Fs[i][i]-Fs[i][j]; } } } }
	 for (i=0; i<N; i++) { // Заполнение матриц коэффициентов q и e
		for (j=0; j<N; j++) {
			if (i==j) ikr=1e0; else ikr=0.0; // Символ Кронекера delta_ij //cout << "i = " << i << "\tj = " << j << "\tk = " << ikr << endl;
			if (fabs(EPS[j])>epsi) 
				qm[i][j]=ikr/EPS[j]-(1e0/EPS[j]-1e0)*Fs[i][j];
			em[i][j]=ikr-Fs[i][j]; } }
	 for (i=0; i<N; i++) { // Заполнение матрицы выходных коэффициентов POUT и RHS
			B[i]=-HOs[i];
			for (j=0; j<N; j++) {
					pm[i][j]=qm[i][j]*idf[j]-em[i][j]*(1e0-idf[j]);
					B[i]=B[i]+(em[i][j]*idf[j]-qm[i][j]*(1e0-idf[j]))*PIN[j]; } }
	 POUT=MetodGaussa(pm, B, N, POUT); //for (i=0; i<N; i++) cout << "pout = " << POUT[i] << endl; cout << endl; POUT=GAUSS(pm, B, N, POUT); //for (i=0; i<N; i++) cout << "pout = " << POUT[i] << "\t"; cout << endl;
	 double sumq = 0.0, *T=new double[N], sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	 for (i=0; i<N; i++) { if (!(ID[i])) T[i]=pow((POUT[i]/sigma),1e0/4e0); else q[i]=POUT[i];
     QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; } //printf("Summa vsekh potokov = %0.10lf\n", sumq); 
	 delete []T; delete []q; delete []QA;
	 for (i=0; i<N; i++) { for (j=0; j<N; j++) printf("%0.6lf\t", Fs[i][j]); cout << endl; }
	 for (i=0; i<N; i++) { p=qm[i]; delete []p; p=em[i]; delete []p; p=pm[i]; delete []p; } 
	 delete []qm; delete []em; delete []pm; delete []B; delete []idf; return POUT;
}
double *GRAYDIFFSPEC(int iclsd, double *A, double *EPS, double *RHOs, double *HOs, double **Fs, int *ID, double *PIN, double *POUT, int chst)
{ int N=chst, i=0, j=0; 
double **qm=new double*[N], **em=new double*[N], **pm=new double*[N], *B=new double[N], *p, *idf=new double[N], hf=1e0, s, ikr, y, epsi=1e-10; // Вычисление недостающих угловых коэффициентов - левой нижней части матрицы из свойства взаимности
for (i=0; i<N; i++) { p=new double[N]; qm[i]=p; p=new double[N]; em[i]=p; p=new double[N]; pm[i]=p; }
for (i=0; i<N; i++) { s=0.0; for (j=0; j<ID[i]; j++) s=s+hf; idf[i]=s; }
for (i=1; i<N; i++)
	 for (j=0; j<i; j++)
		  if (fabs(A[i])>epsi) 
			  Fs[i][j]=(A[j]/A[i])*Fs[j][i]; else Fs[i][j]=0.0; //Для закрытой конфигурации, необходимо рассчитать диагональные элементы из правила суммирования
	if (iclsd==1) { for (i=0; i<N; i++) {
		Fs[i][i]=1e0;
		for (j=0; j<N; j++) {
			if (j!=i) {
				Fs[i][i]=Fs[i][i]-(1e0-RHOs[j])*Fs[i][j]; } }
		if (RHOs[i]<1e0) Fs[i][i]=Fs[i][i]/(1e0-RHOs[i]); else RHOs[i]=0.0; } }
	 for (i=0; i<N; i++) { // Заполнение матриц коэффициентов q и e
		for (j=0; j<N; j++) {
			if (i==j) ikr=1e0; else ikr=0.0; // Символ Кронекера delta_ij //cout << "i = " << i << "\tj = " << j << "\tk = " << ikr << endl;
			if (fabs(EPS[j])>epsi) 
				qm[i][j]=ikr/EPS[j]-((1.0-RHOs[j])/EPS[j]-1.0)*Fs[i][j];
			em[i][j]=ikr-(1.0-RHOs[j])*Fs[i][j]; } }
	 for (i=0; i<N; i++) { // Заполнение матрицы выходных коэффициентов POUT и RHS
			B[i]=-HOs[i];
			for (j=0; j<N; j++) {
					pm[i][j]=qm[i][j]*idf[j]-em[i][j]*(1e0-idf[j]);
					B[i]=B[i]+(em[i][j]*idf[j]-qm[i][j]*(1e0-idf[j]))*PIN[j]; } } //cout << "Fs" << endl; for (i=0; i<N; i++) { for (j=0; j<N; j++) printf("%0.4lf\t",Fs[i][j]); cout << endl; } //for (i=0; i<N; i++) { y=0.0; for (j=0; j<N; j++) y=y+(1e0-RHOs[j])*Fs[i][j]; cout << "stroka ( " << i << " ) = " << y << endl; } //j=getchar();
	 POUT=MetodGaussa(pm, B, N, POUT); //for (i=0; i<N; i++) cout << "pout = " << POUT[i] << endl; cout << endl; POUT=GAUSS(pm, B, N, POUT); for (i=0; i<N; i++) cout << "pout = " << POUT[i] << endl; cout << endl; POUT=reshMetObrMatr(pm, B, N, POUT); for (i=0; i<N; i++) cout << "pout = " << POUT[i] << endl; cout << endl;
	 for (i=0; i<N; i++) { for (j=0; j<N; j++) printf("%0.8lf\t", Fs[i][j]); cout << endl; }
	 double sumq = 0.0, *T=new double[N], sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	for (i=0; i<N; i++) { if (!(ID[i])) T[i]=pow((POUT[i]/sigma),1e0/4e0); else q[i]=POUT[i];
    QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; } //printf("Summa = %0.10lf\n", sumq); 
	delete []T; delete []q; delete []QA;
	for (i=0; i<N; i++) { p=qm[i]; delete []p; p=em[i]; delete []p; p=pm[i]; delete []p; } 
	 delete []qm; delete []em; delete []pm; delete []B; delete []idf; 
	 return POUT; } //Инвертирование матрицы коэффициентов POUT и умножение на матрицу RHS, чтобы получить POUT
double *SEMIGRAY(int iclsd, double *A, double **EPS, double **RHOs, double **HOs, double ***Fs, int *ID, double *q, double *T, double *PIN, double *POUT, int chst, int chel)
{ int N=chst, i, j, k, *id1=new int[N];
double ikr, sigma=5.670E-8, *epsl=new double[N], *rhosl=new double[N], **FSl=new double*[N], *q1=new double[N], *HOs2=new double[N], *p; // Расчет плотности теплового потока  для внешнего излучения (для кажого диапазона)
for (j=0; j<N; j++) { p=new double[N]; FSl[j]=p; }
    for (i=0; i<N; i++) { id1[i]=1; PIN[i]=0.0; } // Установить значения eps, rhos и Fs для k-ого диапазона
	k=0; for (i=0; i<N; i++) { epsl[i]=EPS[k][i]; rhosl[i]=RHOs[k][i];
        for (j=i; j<N; j++) FSl[i][j]=Fs[k][i][j]; }
	GRAYDIFFSPEC(iclsd, A, epsl, rhosl, HOs[k], FSl, id1, PIN, q1, N); //Вычисление неизвестных плотностей потоков тепловой энергии и температур при q2=q-q1
    k=1; for (i=0; i<N; i++) {
        epsl[i]=EPS[k][i]; rhosl[i]=RHOs[k][i];
        for (j=i; j<N; j++) { 
			FSl[i][j]=Fs[k][i][j];
			FSl[j][i]=A[i]/A[j]*FSl[i][j]; }
        if (iclsd==1) {
			FSl[i][i]=1e0;
            for (j=0; j<N; j++) 
				if (j!=i)
                    FSl[i][i]=FSl[i][i]-(1e0-rhosl[j])*FSl[i][j]; }
        HOs2[i]=-q1[i]/EPS[k][i];
        for (j=0; j<N; j++) HOs2[i]=HOs2[i]+((1e0-rhosl[j])/EPS[2][j]-1e0)*FSl[i][j]*q1[j];
        if (!(ID[i])) PIN[i]=q[i];
		else PIN[i]=sigma*pow(T[i],4e0); }
    GRAYDIFFSPEC(iclsd, A, epsl, rhosl, HOs2, FSl, ID, PIN, POUT, N);
    for (i=0; i<N; i++) {
        if (!(ID[i])) T[i]=pow((POUT[i]/sigma),1e0/4e0);      
        else q[i]=POUT[i]; } 
	delete []id1; delete []epsl; delete []rhosl; delete []q1; delete []HOs2;
	for (j=0; j<N; j++) { p=FSl[j]; delete []p; } delete []FSl;
	return POUT;
}
double VIEW(int NO, int NARG, double *ARG, int fl)
{ double VIEW=0.0, epsi=1e-15, e=epsi, ov=1e0/2e0, A=0.0, B=0.0, C=0.0; int k=0;
if ((NO<38) || (NO>39)) { cout << "Illegal value for number (NO =" << NO << ") " << endl; k=getchar(); exit(1); }
if (NO==38) { if (NARG!=3) {cout << "Wrong number of input parameters (NARG =" << NARG << ") for NO = " << NO << endl; k=getchar(); exit(1);}
double X=0.0, Y=0.0, RTX=0.0, RTY=0.0, RT=0.0;
      k=0; A=ARG[k]; k++; B=ARG[k]; k++; C=ARG[k];
	  if (fabs(C)>e) { X=A/C; Y=B/C; } else { X=0.0; Y=0.0; }
	  RTX=pow(1.0+X*X, ov); RTY=pow(1.0+Y*Y, ov); RT=pow(1.0+X*X+Y*Y, ov);
      if ((fabs(X)>e) && (fabs(Y)>e))     
	  VIEW=(log(RTX*RTY/RT)+X*RTY*atan(X/RTY)+Y*RTX*atan(Y/RTX)-X*atan(X)-Y*atan(Y))*2.0/(pi*X*Y);
	  else VIEW=0.0; }
else if (NO==39) { if(NARG!=3) { cout << "Wrong number of input parameters (NARG =" << NARG << ") for NO = " << NO << endl; k=getchar(); exit(1); }
	double H=0.0, L=0.0, W=0.0, HH=0.0, WW=0.0, W2=0.0, H2=0.0, HW2=0.0, HW=0.0, H12=0.0, W12=0.0, C1=0.0, C2=0.0, C3=0.0;
      k=0; H=ARG[k]; k++; L=ARG[k]; k++; W=ARG[k];
	  if (fabs(L)>e) { HH=H/L; WW=W/L; } else { HH=0.0; WW=0.0; } 
	  W2=WW*WW; H2=HH*HH; HW2=H2+W2; HW=sqrt(H2+W2); H12=H2+1.0; W12=W2+1.0;
	  if (HW>e) { C2=W2*(H12+W2)/W12/HW2; C3=H2*(H12+W2)/H12/HW2; } else { C2=0.0; C3=0.0; }
      C1=W12*H12/(H12+W2);  
	  double atw=0.0, ath=0.0, athw=0.0;
	  if (fabs(WW)>e) atw=atan(1e0/WW); else atw=pi/2e0; 
	  if (fabs(HH)>e) ath=atan(1e0/HH); else ath=pi/2e0; 
	  if (fabs(HW)>e) athw=atan(1e0/HW); else athw=pi/2e0;
      if (fabs(WW)>e) VIEW=(WW*atw+HH*ath-HW*athw+0.25*(log(C1)+W2*log(C2)+H2*log(C3)))/(pi*WW); else VIEW=0.0; //cout << "H = " << HH << "\tW = " << WW << "\tHW = " << HW << "\tpi = " << pi << "\tC1 = " << C1 << "\tC2 = " << C2 << "\tC3 = " << C3 << "atw = " << atw << "\tath = " << ath << "\tathw = " << athw << "\tV = " << VIEW << endl; 
} //if ((fl==13) || (fl==15) || (fl==4)) { cout << "UK = " << VIEW << "\t"; cout << "A = " << A << "\tB = " << B << "\tC = " << C << "\t"; }
return VIEW; }
double perpplates(double X1, double X2, double Y1, double Y2, double Z1, double Z2, double Z3)
{ int NARG=3, k=3, n=39, fl=0; double perppltf=0.0, A=0.0, F=0.0, *ARG=new double[k], epsi=1e-15, e=epsi;
      k=0; ARG[k]=Y2; k++; ARG[k]=Z3; k++; ARG[k]=X2;
      A=X2*Z3;
      k=39; fl++; if (fabs(A*Y2)>e) F=A*VIEW(k, NARG, ARG, fl); //cout << "F1 per = " << F << "\t";
	  k=0; ARG[k]=Y1; A=X2*Z3; fl++;
      if (fabs(A*Y1)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F2 per = " << F << "\t";
      k=0; ARG[k]=Y2; k=2; ARG[k]=X1; A=X1*Z3; fl++;
      if (fabs(A*Y2)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F3 per = " << F << "\t";
	  k=0; ARG[k]=Y1; A=X1*Z3; fl++;
      if (fabs(A*Y1)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F4 per = " << F << "\t";
      k=0; ARG[k]=Y2; k++; ARG[k]=Z2; A=X1*Z2; fl++;
      if (fabs(A*Y2)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F5 per = " << F << "\t";
      k=0; ARG[k]=Y1; A=X1*Z2; fl++;
      if (fabs(A*Y1)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F6 per = " << F << "\t";
      k=0; ARG[k]=Y2; k=2; ARG[k]=X2; A=X2*Z2; fl++;
      if (fabs(A*Y2)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F7 per = " << F << "\t";
      k=0; ARG[k]=Y1; A=X2*Z2; fl++;
      if (fabs(A*Y1)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F8 per = " << F << "\t";
      k=0; ARG[k]=Y2; k++; ARG[k]=(Z3-Z1); A=X2*(Z3-Z1); fl++;
      if (fabs(A*Y2)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F9 per = " << F << "\t";
      k=0; ARG[k]=Y1; A=X2*(Z3-Z1); fl++;
      if (fabs(A*Y1)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F10 per = " << F << "\t";
      k=0; ARG[k]=Y2; k=2; ARG[k]=X1; A=X1*(Z3-Z1); fl++;
      if (fabs(A*Y2)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F11 per = " << F << "\t";
      k=0; ARG[k]=Y1; A=X1*(Z3-Z1); fl++;
      if (fabs(A*Y1)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F12 per = " << F << "\t";
      k=0; ARG[k]=Y2; k++; ARG[k]=(Z2-Z1); k++; ARG[k]=X2; A=X2*(Z2-Z1); fl++;
      if (fabs(A*Y2)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F13 per = " << F << "\t";
      k=0; ARG[k]=Y1; A=X2*(Z2-Z1); fl++;
	  if (fabs(A*Y1)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F14 per = " << F << "\t";
      k=0; ARG[k]=Y2; k=2; ARG[k]=X1; A=X1*(Z2-Z1); fl++;
      if (fabs(A*Y2)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F15 per = " << F << "\t";
      k=0; ARG[k]=Y1; A=X1*(Z2-Z1); fl++;
      if(fabs(A*Y1)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F16 per = " << F << "\n";
      if ((fabs(Z1)>e) && (fabs(X2-X1)>e)) perppltf=F/(2.0*(X2-X1)*Z1); else perppltf=0.0;
	  if (ARG) delete []ARG; 
	  return perppltf; }
double parlplates(double X1, double X2, double X3, double Y1, double Y2, double Y3, double C)
{ int NARG=3, j=0, k=0, n=38, fl=0; double parlpltf, A=0.0, F=0.0, *ARG=new double[NARG], epsi=1e-10, e=epsi;
      k=0; ARG[k]=X3; k++; ARG[k]=Y3; k++; ARG[k]=C; k=0; j=1; A=ARG[k]*ARG[j]; fl++;
      if (fabs(A)>e) F=A*VIEW(n, NARG, ARG, fl); //cout << "F1 par = " << F << "\t";
      k=1; ARG[k]=Y2; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F2 par = " << F << "\t";
	  k=1; ARG[k]=Y3-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F3 par = " << F << "\t";
	  k=1; ARG[k]=Y2-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F4 par = " << F << "\t";
	  k=0; ARG[k]=X2; j=1; ARG[j]=Y3; A=ARG[k]*ARG[j]; fl++;
	  if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F5 par = " << F << "\t";
      k=1; ARG[k]=Y2; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F6 par = " << F << "\t";
	  k=1; ARG[k]=Y3-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F7 par = " << F << "\t";
      k=1; ARG[k]=Y2-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F8 par = " << F << "\t";
	  k=0; ARG[k]=X3-X1; j=1; ARG[j]=Y3; A=ARG[k]*ARG[j]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F9 par = " << F << "\t";
      k=1; ARG[k]=Y2; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F10 par = " << F << "\t";
	  k=1; ARG[k]=Y3-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F11 par = " << F << "\t";
      k=1; ARG[k]=Y2-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F12 par = " << F << "\t";
	  k=0; ARG[k]=X2-X1; j=1; ARG[j]=Y3; A=ARG[k]*ARG[j]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F13 par = " << F << "\t";
      k=1; ARG[k]=Y2; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F14 par = " << F << "\t";
	  k=1; ARG[k]=Y3-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F15 par = " << F << "\t";
      k=1; ARG[k]=Y2-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F16 par = " << F << "\n";
      if ((fabs(X1)>e) && (fabs(Y1)>e)) parlpltf=F/(4.0*(X1*Y1));
      if (ARG) delete []ARG; 
	  return parlpltf; }
double *GAUSS(double **A, double *B, int chst, double *X)
{ int N=chst, I, *L=new int[N], K, J, LK;
double *S=new double[N], SMAX, RMAX, R, XMULT, SUM, epsi=1e-10;
if ((!S) || (!L)) { cout << "No memory!" << endl; K=getchar(); exit(1); }
      for (I=0; I<N; I++)
			{ L[I]=I; SMAX=0.0; //L[I] - номер строки
			for (J=0; J<N; J++)
				if (fabs(A[I][J])>SMAX) 
					SMAX=fabs(A[I][J]);
			S[I]=SMAX; } //S[I] - массив максимальных элементов в строке
      for (K=0; K<N-1; K++) 
        { RMAX=0.0; J=K;
        for (I=K; I<N; I++) 
		{ if (fabs(S[L[I]])>epsi) 
				R=fabs(A[L[I]][K])/S[L[I]]; else R=0.0; //делим первые элементы на S[I]
		if (R>RMAX) { J=I; RMAX=R; } } //находим 'элемент с максимальным отношением K-го элемента к максимальному в K-ой строке, записываем его номер
		LK=L[J]; L[J]=L[K]; L[K]=LK; //меняем местами L[K] - текущий и L[J] - с максимальным
		for (I=K+1; I<N; I++)
			{ XMULT = A[L[I]][K]/A[LK][K]; //делим каждый элемент K-ого столбца на максимальный, начиная с (K+1)-го элемента строки 
					if (A[LK][K]==0.0)
						XMULT=0.0;
				for (J=K+1; J<N; J++)
					A[L[I]][J]=A[L[I]][J]-XMULT*A[LK][J]; //вычитаем, чтобы K-ый элемент в K-ом столбце обнулился
				A[L[I]][K] = XMULT; } }
      for (K=0; K<N-1; K++)
        for (I=K+1; I<N; I++)      
          B[L[I]] = B[L[I]] - A[L[I]][K]*B[L[K]];
      if (fabs(A[L[N-1]][N-1])>epsi) 
		  X[N-1] = B[L[N-1]]/A[L[N-1]][N-1]; else X[N-1]=0.0;
      for (I=N-1; I>=0; I--) {
		  SUM = B[L[I]];
		  for (J=I+1; J<N; J++)
          SUM = SUM - A[L[I]][J]*X[J];
        if (fabs(A[L[I]][I])>epsi) 
			X[I] = SUM/A[L[I]][I];
		else X[I]=0.0; }
	  delete []L; delete []S; return X; }
double SeryeStenkiRasIzlDifPov(double ww, double hh, double ll, double *T, double *EPS, double *HOs, double *RHOs, double *A, int *ns, int chst) //N серых диффузных поверхностей
{ int N=chst, fl=0; double b=0.0, t=0.0, epsi=1e-10; //t=T[ns[0]]; b=T[ns[1]]; if (fabs(t-b)<epsi) return 0.0;
double **Fs=new double*[N], z, r, y; int i=0, j=0, f=3, k=0, d=0; 
double *PIN=new double[N], *POUT=new double[N], *arg=new double[f], *p=NULL;
for (i=0; i<N; i++) { p=new double[N]; for (j=0; j<N; j++) p[j]=0.0; Fs[i]=p; }
	double sumq=0.0, sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	int *id=new int[N], iclsd=0; for (i=0; i<N; i++) id[i]=1; // Поверхности 1 и 3 (bottom and top)
    	for (i=0; i<N; i++) if (!(id[i])) PIN[i]=q[i]; else PIN[i]=sigma*pow(T[i], 4e0); //Заполнение массива PIN ПТП q и температурой T - перевод температур в ППИ, Fs(1,2)=F1-2, 
	iclsd=1; //угловые коэффициенты; конфигурация замкнута (iclsd=1), диагональные элементы не нужны //for (j=0; j<N; j++) cout << "RHOs ( " << j << " ) = " << RHOs[j] << endl;  for (j=0; j<N; j++) cout << "EPS ( " << j << " ) = " << EPS[j] << endl;
	k=0; arg[k]=fabs(hh); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(ww); //h,l,w
	k=39; i=0; j=1; Fs[i][j]=VIEW(k, f, arg, fl); //Fs12 //cout << "Fs(0,1) = " << Fs[0][1] << "\th = " << hh << "\tl = " << ll << "\tw = " << ww << "\tf = " << f << endl; //Fs(1,2)=F1-2, F12
    	k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(hh);
	k=38; i=0; j=2; Fs[i][j]=VIEW(k, f, arg, fl); //Fs13
    	k=3; d=1; Fs[i][k]=Fs[i][d]; //Fs14
    	k=0; arg[k]=fabs(hh); k++; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); 
    k=39; i=0; j=4; Fs[i][j]=VIEW(k, f, arg, fl); //Fs15
	k=5; Fs[i][k]=Fs[i][j]; //Fs16
	   	k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(hh); 
    k=39; i=1; j=2; Fs[i][j]=VIEW(k, f, arg, fl); //Fs23
    	k=0; arg[k]=fabs(ll); k++; arg[k]=fabs(hh); k++; arg[k]=fabs(ww); 
    k=38; i=1; j=3; Fs[i][j]=VIEW(k, f, arg, fl); //Fs24
    	k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(hh); k++; arg[k]=fabs(ll);
    k=39; i=1; j=4; Fs[i][j]=VIEW(k, f, arg, fl); //F25 
    i=1; j=5; k=4; Fs[i][j]=Fs[i][k]; //Fs26
	i=2; j=3; k=0; Fs[i][j]=Fs[k][j]; //Fs34
	i=2; j=4; k=0; Fs[i][j]=Fs[k][j]; //Fs35
	i=2; j=5; k=0; Fs[i][j]=Fs[k][j]; //Fs36
	i=3; j=4; k=1; Fs[i][j]=Fs[k][j]; //Fs45
	i=3; j=5; k=1; Fs[i][j]=Fs[k][j]; //Fs46
	k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(hh); k++; arg[k]=fabs(ll);
    k=38; i=4; j=5; Fs[i][j]=VIEW(k, f, arg, fl); //F56 
	POUT=GRAYDIFF(iclsd, A, EPS, HOs, Fs, id, PIN, POUT, N);
	sumq=0.0; //проверка полной ППИ
    	for (i=0; i<N; i++) //Вывод - перевод в температуры
		{ if (!(id[i])) T[i]=pow((POUT[i]/sigma), 1e0/4e0); else q[i]=POUT[i];
        	QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; 
			cout << "Surface " << i << "\tT [K] = " << T[i] << "\tq [W/m2] = " << q[i] << "\tQA [W] = " << QA[i]/ll << endl;
		} r=0.0; z=1.0; for (i=0; i<N; i++) if (ns[i]==i) { r=r+z*q[i]; z=-z; } 
		for (i=0; i<N; i++) { p=Fs[i]; delete []p; } delete []Fs; delete []PIN; 
		delete []POUT; delete []arg; delete []q; delete []QA; delete []id; 
		cout << "Summa vsekh teplovykh potokov = " << sumq << endl;
		return r; }
double F0_lamT(double lT0)
{ int k=20, j=0; double *F0=new double[k], *lamT=new double[k], dolya=0.0, hf=1e0;
if ((!F0) || (!lamT)) { cout << "No memory!" << endl; j=getchar(); exit(1); }
lamT[j]=0.18889e-2; F0[j]=0.05059; j++; lamT[j]=0.22222e-2; F0[j]=0.10503; j++; lamT[j]=0.24444e-2; F0[j]=0.14953; j++;
lamT[j]=0.27222e-2; F0[j]=0.21033; j++; lamT[j]=0.28889e-2; F0[j]=0.24803; j++; lamT[j]=0.31667e-2; F0[j]=0.31067; j++;
lamT[j]=0.33889e-2; F0[j]=0.35933; j++; lamT[j]=0.36111e-2; F0[j]=0.40585; j++; lamT[j]=0.38889e-2; F0[j]=0.46031; j++;
lamT[j]=0.41111e-2; F0[j]=0.50066; j++; lamT[j]=0.44444e-2; F0[j]=0.55573; j++; lamT[j]=0.47778e-2; F0[j]=0.60449; j++;
lamT[j]=0.51667e-2; F0[j]=0.65402; j++; lamT[j]=0.56111e-2; F0[j]=0.70211; j++; lamT[j]=0.61667e-2; F0[j]=0.75146; j++;
lamT[j]=0.68889e-2; F0[j]=0.80152; j++; lamT[j]=0.78889e-2; F0[j]=0.85171; j++; lamT[j]=0.93889e-2; F0[j]=0.90031; j++;
lamT[j]=1.25556e-2; F0[j]=0.95094; j++; lamT[j]=2.33333e-2; F0[j]=0.99051; j++;
dolya=opredKTPTKTochSha(F0, lamT, lT0, k); if (dolya<0.0) dolya=0.0; if (dolya>=hf) dolya=hf; 
delete []F0; delete []lamT; return dolya; }
double bbfn(double X)
{   double CC=1.5e1/pow(pi,4e0), C2=1.4388e4;
	double EPS=1e-16, V=C2/X, EX=exp(V), M=0.0, BBFN=0.0, VM, BM, EM=1e0, hf=EM; 
	int j=0, jk=1000;
do {
    M=M+hf;
    VM=M*V;
    BM=(6e0+VM*(6e0+VM*(3e0+VM)))/pow(M, 4e0);
    EM=EM/EX;
    BBFN=BBFN+BM*EM; j++; }
while (((pow(VM,3e0)*EM)>EPS) && (j<jk));
    BBFN=CC*BBFN;
    return BBFN; }