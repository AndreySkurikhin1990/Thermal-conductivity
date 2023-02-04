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
const double pi = acos(-1e0), enu=1e-3, tocrasov=1e-4;
const int dsov=60, dmkoscrt=14;
//-----
int dlar=1, dlarsh=1, dlari=1, dlark=1, N=6;
char *sndov=NULL, *ssdov=NULL, *skpovsk=NULL, *sdvovsk=NULL, *szfovk=NULL, *snmov=NULL;
char *sfnoov=NULL, *sppovv=NULL, *sppovvk=NULL, *sppovs=NULL, *sppovsk=NULL;
char *skpovs=NULL, *sdvovs=NULL, *szfov=NULL, *sdvovv=NULL, *sdvovvk=NULL, *skpovv=NULL;
char *skpovvk=NULL, *sdvovi=NULL, *sppovi=NULL, *skpovi=NULL, *sdvovik=NULL;
char *sppovik=NULL, *skpovik=NULL, *skpovk=NULL, *sdvovk=NULL, *sppovk=NULL;
char *sdvovkk=NULL, *sppovkk=NULL, *skpovkk=NULL;
double *koefoslab(double, double, double, double *, int, double *);
double epsisredver(double, double *, double *, int, double *, double *, int);
double epsisreditom(double, double *, double *, int, double *, double *, int, int);
double epsisredkvi(double, double *, double *, int, double *, double *, int, int);
double epsisredsha(double, double *, double *, int, double *, double *, int);
double *izmMasChast(double *, double, int, double *);
double trapz(double *, double *, int);
double *Kramers_Kronig_ver(double *, int, double *);
double *Kramers_Kronig_sha(double *, int, double *);
double *Kramers_Kronig_itom(double *, int, int, double *);
double *Kramers_Kronig_kvi(double *, int, double *, int);
double *dliny_voln_ver(double *, int);
double *dliny_voln_sham(double *, int);
double *dliny_voln_itom(double*, int, int);
double *dliny_voln_kvi(double *, int);
double *epsilnu(double *, int, double *);
double *Koef_Pogl_ver(double *);
double *Koef_Pogl_sha(double *);
double *Koef_Pogl_itom(double *, int);
double *Koef_Pogl_kvi(double *);
double *PokazPrelomAl(double *, double *, double *, double, int, double *);
double *PokazPrelomPribl(double *, double *, double *, int, double *, double);
void zapisvfile(double *, int, char *);
double *ronusha(double *, int);
double *ronukvi(double *, int);
double reflsha(double);
double reflkvi(double, int);
double reflitom(double);
void napstrdir(int, int);
double *reshMetKram(double **, double *, int);
double *oslaintestepcher(int, int, double *);
double *RasKorOV(double, double, double, double, double, double);
double UsredMasOV(double **, int, int);
void vyvomatr(double **, int, int);
double KoefPoglRossel(double, double *, double *, double *, int);
double *KoefPoglRosselNac(double *, int, int, double, double, double, double *, double *, int, double, double, double *, double *, int, int, int);
double nsredPlank2(double *, double *, double, int);
void osvpamov(int);
double usredVelichPlank(double *, double *, double *, double, int, double);
double LuchKTPChudnovsky(double *, double, int, double);
double PoiskReflPhi(double, double, double, double);
double ReflSredSha(double);
double ReflSredVer(double);
double ReflSreditom(double, int);
double ReflSredkvi(double, int);
double PoiskTetaShtrNach(double, double, double, double);
double PraCha10(double, double, double, double);
double *PoiSpeKoeOtr(double *, int);
double provUgla(double);
double *GRAYDIFFSPEC(int, double *, double *, double *, double *, double **, int *, double *, double *, int);
double *GRAYDIFF(int, double *, double *, double *, double **, int *, double *, double *, int);
double VIEW(int, int, double *);
double perppltf(double, double, double, double, double, double, double);
double parlpltf(double, double, double, double, double, double, double);
double *GAUSS(double **, double *, int, double *);
double SeryeStenkiRasIzl(double, double, double, double *, double *, double *, double *, double *, int *, int);
double opredKTPTKTochSha(double *, double *, double, int);
double bbfn(double);
double *MetodGaussa(double **, double *, int, double *);
double *reshMetObrMatr(double **, double *, int, double *);
double *SEMIGRAY(int, double *, double **, double **, double **, double ***, int *, double *, double *, double *, double *, int, int);
//-------------------
void napstrdir(int vyve, int vm)
{	int k; 
sndov=new char[dsov]; ssdov=new char[dsov]; szfovk=new char[dsov]; 
snmov=new char[dsov]; sfnoov=new char[dsov]; szfov=new char[2*dsov];
for (k=0; k<(2*dsov); k++) { szfov[k]='\0'; } 
for (k=0; k<dsov; k++) { sndov[k]='\0'; ssdov[k]='\0'; szfovk[k]='\0'; snmov[k]='\0'; sfnoov[k]='\0'; }
k=0; 
sndov[k]='D'; k++; sndov[k]=':'; k++; sndov[k]='\\';k++; sndov[k]='\\'; k++; sndov[k]='_';  k++; sndov[k]='À';  k++;
sndov[k]='ñ'; k++; sndov[k]='ï'; k++; sndov[k]='è'; k++; sndov[k]='ð';  k++; sndov[k]='à';  k++; sndov[k]='í';  k++;
sndov[k]='ò'; k++; sndov[k]='ó'; k++; sndov[k]='ð'; k++; sndov[k]='à';  k++; sndov[k]='\\'; k++; sndov[k]='\\'; k++;
sndov[k]='t'; k++; sndov[k]='m'; k++; sndov[k]='p'; k++; sndov[k]='\\'; k++; sndov[k]='\\'; k++; sndov[k]='\0';
k=0;
ssdov[k]='C';  k++; ssdov[k]=':'; k++; ssdov[k]='\\'; k++; ssdov[k]='\\'; k++; ssdov[k]='U';  k++; ssdov[k]='s';  k++;
ssdov[k]='e';  k++; ssdov[k]='r'; k++; ssdov[k]='s';  k++; ssdov[k]='\\'; k++; ssdov[k]='\\'; k++; ssdov[k]='À';  k++;
ssdov[k]='í';  k++; ssdov[k]='ä'; k++; ssdov[k]='ð';  k++; ssdov[k]='å';  k++; ssdov[k]='é';  k++; ssdov[k]='\\'; k++;
ssdov[k]='\\'; k++; ssdov[k]='D'; k++; ssdov[k]='o';  k++; ssdov[k]='c';  k++; ssdov[k]='u';  k++; ssdov[k]='m';  k++;
ssdov[k]='e';  k++; ssdov[k]='n'; k++; ssdov[k]='t';  k++; ssdov[k]='s';  k++; ssdov[k]='\\'; k++; ssdov[k]='\\'; k++;
ssdov[k]='_';  k++; ssdov[k]='À'; k++; ssdov[k]='ñ';  k++; ssdov[k]='ï';  k++; ssdov[k]='è';  k++; ssdov[k]='ð';  k++;
ssdov[k]='à';  k++; ssdov[k]='í'; k++; ssdov[k]='ò';  k++; ssdov[k]='ó';  k++; ssdov[k]='ð';  k++; ssdov[k]='à';  k++;
ssdov[k]='\\'; k++; ssdov[k]='\\';k++; ssdov[k]='t';  k++; ssdov[k]='m';  k++; ssdov[k]='p';  k++; ssdov[k]='\\'; k++;
ssdov[k]='\\'; k++; ssdov[k]='\0';
k=0;
szfovk[k]='V'; k++; szfovk[k]='y'; k++; szfovk[k]='v'; k++; szfovk[k]='o'; k++; szfovk[k]='d'; k++; szfovk[k]='v'; k++; 
szfovk[k]='F'; k++; szfovk[k]='i'; k++; szfovk[k]='l'; k++; szfovk[k]='e'; k++; szfovk[k]='.'; k++; szfovk[k]='t'; k++;
szfovk[k]='x'; k++; szfovk[k]='t'; k++; szfovk[k]='\0';
k=0;
snmov[k]='N'; k++; snmov[k]='o'; k++; snmov[k]='_'; k++; snmov[k]='m'; k++; snmov[k]='e'; k++; snmov[k]='m'; k++; 
snmov[k]='o'; k++; snmov[k]='r'; k++; snmov[k]='y'; k++; snmov[k]='!'; k++; snmov[k]='\0';
k=0;
sfnoov[k]='F';  k++; sfnoov[k]='i';  k++; sfnoov[k]='l';  k++; sfnoov[k]='e';  k++; sfnoov[k]='_';  k++; sfnoov[k]='i'; k++;
sfnoov[k]='s';  k++; sfnoov[k]='_';  k++; sfnoov[k]='n';  k++; sfnoov[k]='o';  k++; sfnoov[k]='t';  k++; sfnoov[k]='_'; k++;
sfnoov[k]='o';  k++; sfnoov[k]='p';  k++; sfnoov[k]='e';  k++; sfnoov[k]='n';  k++; sfnoov[k]='!';  k++; sfnoov[k]='\0'; 
strcpy(szfov,ssdov); strcat(szfov,szfovk); k=strlen(szfov)+1; szfov[k]='\0';
//0 - øàìîò
if (!vyve) { skpovsk=new char[dsov]; sdvovsk=new char[dsov]; sppovsk=new char[dsov]; 
skpovs=new char[2*dsov]; sdvovs=new char[2*dsov]; sppovs=new char[2*dsov]; 
if ((!skpovsk) || (!sdvovsk) || (!sppovsk) || (!skpovs) || (!sppovs) || (!sdvovs)) 
	{ cout << snmov << endl; k=getchar(); exit(1); } 
if ((!sndov) || (!ssdov) || (!szfovk) || (!snmov) || (!sfnoov) || (!szfov)) { cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<(2*dsov); k++) { skpovs[k]='\0'; sdvovs[k]='\0'; sppovs[k]='\0'; } 
for (k=0; k<dsov; k++) { skpovsk[k]='\0'; sdvovsk[k]='\0'; sppovsk[k]='\0'; }
k=0;
skpovsk[k]='K'; k++; skpovsk[k]='o'; k++; skpovsk[k]='e'; k++; skpovsk[k]='f'; k++; skpovsk[k]='f'; k++; skpovsk[k]='i'; k++;
skpovsk[k]='c'; k++; skpovsk[k]='i'; k++; skpovsk[k]='e'; k++; skpovsk[k]='n'; k++; skpovsk[k]='t'; k++; skpovsk[k]='_'; k++;
skpovsk[k]='p'; k++; skpovsk[k]='o'; k++; skpovsk[k]='g'; k++; skpovsk[k]='l'; k++; skpovsk[k]='o'; k++; skpovsk[k]='s'; k++;
skpovsk[k]='c'; k++; skpovsk[k]='h'; k++; skpovsk[k]='e'; k++; skpovsk[k]='n'; k++; skpovsk[k]='i'; k++; skpovsk[k]='y'; k++;
skpovsk[k]='a'; k++; skpovsk[k]='_'; k++; skpovsk[k]='s'; k++; skpovsk[k]='h'; k++; skpovsk[k]='a'; k++; skpovsk[k]='.'; k++;
skpovsk[k]='t'; k++; skpovsk[k]='x'; k++; skpovsk[k]='t'; k++; skpovsk[k]='\0';
k=0;
sdvovsk[k]='D'; k++; sdvovsk[k]='l'; k++; sdvovsk[k]='i'; k++; sdvovsk[k]='n'; k++; sdvovsk[k]='a'; k++; sdvovsk[k]='_'; k++; 
sdvovsk[k]='v'; k++; sdvovsk[k]='o'; k++; sdvovsk[k]='l'; k++; sdvovsk[k]='n'; k++; sdvovsk[k]='y'; k++; sdvovsk[k]='_'; k++; 
sdvovsk[k]='s'; k++; sdvovsk[k]='h'; k++; sdvovsk[k]='a'; k++; sdvovsk[k]='.'; k++; sdvovsk[k]='t'; k++; sdvovsk[k]='x'; k++; 
sdvovsk[k]='t'; k++; sdvovsk[k]='\0';
k=0;
sppovsk[k]='P'; k++; sppovsk[k]='o'; k++; sppovsk[k]='k'; k++; sppovsk[k]='a'; k++; sppovsk[k]='z'; k++; sppovsk[k]='a'; k++; 
sppovsk[k]='t'; k++; sppovsk[k]='_'; k++; sppovsk[k]='p'; k++; sppovsk[k]='r'; k++; sppovsk[k]='e'; k++; sppovsk[k]='l'; k++;
sppovsk[k]='o'; k++; sppovsk[k]='m'; k++; sppovsk[k]='l'; k++; sppovsk[k]='e'; k++; sppovsk[k]='n'; k++; sppovsk[k]='_'; k++;
sppovsk[k]='s'; k++; sppovsk[k]='h'; k++; sppovsk[k]='a'; k++; sppovsk[k]='.'; k++; sppovsk[k]='t'; k++; sppovsk[k]='x'; k++; 
sppovsk[k]='t'; k++; sppovsk[k]='\0';
strcpy(sppovs,ssdov); strcat(sppovs,sppovsk); k=strlen(sppovs)+1; sppovs[k]='\0';
strcpy(skpovs,ssdov); strcat(skpovs,skpovsk); k=strlen(skpovs)+1; skpovs[k]='\0'; 
strcpy(sdvovs,ssdov); strcat(sdvovs,sdvovsk); k=strlen(sdvovs)+1; sdvovs[k]='\0'; 
}
//1 - âåðìèêóëèò
else if (vyve==1) { skpovvk=new char[dsov]; sdvovvk=new char[dsov]; sppovvk=new char[dsov]; 
skpovv=new char[2*dsov]; sdvovv=new char[2*dsov]; sppovv=new char[2*dsov]; 
if ((!skpovvk) || (!sdvovvk) || (!sppovvk) || (!skpovv) || (!sdvovv) || (!sppovv)) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<(2*dsov); k++) { skpovv[k]='\0'; sdvovv[k]='\0'; sppovv[k]='\0'; } 
for (k=0; k<dsov; k++) { skpovvk[k]='\0'; sdvovvk[k]='\0'; sppovvk[k]='\0'; }
k=0;
skpovvk[k]='K'; k++; skpovvk[k]='o'; k++; skpovvk[k]='e'; k++; skpovvk[k]='f'; k++; skpovvk[k]='f'; k++; skpovvk[k]='i'; k++;
skpovvk[k]='c'; k++; skpovvk[k]='i'; k++; skpovvk[k]='e'; k++; skpovvk[k]='n'; k++; skpovvk[k]='t'; k++; skpovvk[k]='_'; k++;
skpovvk[k]='p'; k++; skpovvk[k]='o'; k++; skpovvk[k]='g'; k++; skpovvk[k]='l'; k++; skpovvk[k]='o'; k++; skpovvk[k]='s'; k++;
skpovvk[k]='c'; k++; skpovvk[k]='h'; k++; skpovvk[k]='e'; k++; skpovvk[k]='n'; k++; skpovvk[k]='i'; k++; skpovvk[k]='y'; k++;
skpovvk[k]='a'; k++; skpovvk[k]='_'; k++; skpovvk[k]='v'; k++; skpovvk[k]='e'; k++; skpovvk[k]='r'; k++; skpovvk[k]='.'; k++;
skpovvk[k]='t'; k++; skpovvk[k]='x'; k++; skpovvk[k]='t'; k++; skpovvk[k]='\0';
k=0;
sdvovvk[k]='D'; k++; sdvovvk[k]='l'; k++; sdvovvk[k]='i'; k++; sdvovvk[k]='n'; k++; sdvovvk[k]='y'; k++; sdvovvk[k]='_'; k++; 
sdvovvk[k]='v'; k++; sdvovvk[k]='o'; k++; sdvovvk[k]='l'; k++; sdvovvk[k]='n'; k++; sdvovvk[k]='_'; k++; sdvovvk[k]='v'; k++; 
sdvovvk[k]='e'; k++; sdvovvk[k]='r'; k++; sdvovvk[k]='.'; k++; sdvovvk[k]='t'; k++; sdvovvk[k]='x'; k++; sdvovvk[k]='t'; k++;
sdvovvk[k]='\0'; 
k=0;
sppovvk[k]='P'; k++; sppovvk[k]='o'; k++; sppovvk[k]='k'; k++; sppovvk[k]='a'; k++; sppovvk[k]='z'; k++; sppovvk[k]='a'; k++;
sppovvk[k]='t'; k++; sppovvk[k]='_'; k++; sppovvk[k]='p'; k++; sppovvk[k]='r'; k++; sppovvk[k]='e'; k++; sppovvk[k]='l'; k++;
sppovvk[k]='o'; k++; sppovvk[k]='m'; k++; sppovvk[k]='l'; k++; sppovvk[k]='e'; k++; sppovvk[k]='n'; k++; sppovvk[k]='_'; k++;
sppovvk[k]='v'; k++; sppovvk[k]='e'; k++; sppovvk[k]='r'; k++; sppovvk[k]='.'; k++; sppovvk[k]='t'; k++; sppovvk[k]='x'; k++;
sppovvk[k]='t'; k++; sppovvk[k]='\0';
strcpy(sppovv,ssdov); strcat(sppovv,sppovvk); k=strlen(sppovv)+1; sppovv[k]='\0'; //cout << sppovv << endl;
strcpy(skpovv,ssdov); strcat(skpovv,skpovvk); k=strlen(skpovv)+1; skpovv[k]='\0'; //cout << skpovv << endl;
strcpy(sdvovv,ssdov); strcat(sdvovv,sdvovvk); k=strlen(sdvovv)+1; sdvovv[k]='\0'; //cout << sdvovv << endl;
}
//2 - ÈÒÎÌ
else if (vyve==2) { sdvovik=new char[dsov]; sppovik=new char[dsov]; skpovik=new char[dsov];
	sdvovi=new char[2*dsov]; sppovi=new char[2*dsov]; skpovi=new char[2*dsov]; 
if ((!sdvovi) || (!sppovi) || (!skpovi) || (!sdvovik) || (!sppovik) || (!skpovik)) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<(2*dsov); k++) { sdvovi[k]='\0'; sppovi[k]='\0'; skpovi[k]='\0'; } 
for (k=0; k<dsov; k++) { sdvovik[k]='\0'; sppovik[k]='\0'; skpovik[k]='\0'; }
k=0;
sdvovik[k]='D'; k++; sdvovik[k]='l'; k++; sdvovik[k]='i'; k++; sdvovik[k]='n'; k++; sdvovik[k]='y'; k++; sdvovik[k]='_'; k++;
sdvovik[k]='v'; k++; sdvovik[k]='o'; k++; sdvovik[k]='l'; k++; sdvovik[k]='n'; k++; sdvovik[k]='_'; k++; sdvovik[k]='i'; k++;
sdvovik[k]='t'; k++; sdvovik[k]='o'; k++; sdvovik[k]='m'; k++; sdvovik[k]='.'; k++; sdvovik[k]='t'; k++; sdvovik[k]='x'; k++; 
sdvovik[k]='t'; k++; sdvovik[k]='\0'; k++; 
k=0;
sppovik[k]='P'; k++; sppovik[k]='o';  k++; sppovik[k]='k'; k++; sppovik[k]='a'; k++; sppovik[k]='z'; k++; sppovik[k]='a'; k++; 
sppovik[k]='t'; k++; sppovik[k]='e';  k++; sppovik[k]='l'; k++; sppovik[k]='_'; k++; sppovik[k]='p'; k++; sppovik[k]='r'; k++; 
sppovik[k]='e'; k++; sppovik[k]='l';  k++; sppovik[k]='o'; k++; sppovik[k]='m'; k++; sppovik[k]='l'; k++; sppovik[k]='e'; k++; 
sppovik[k]='n'; k++; sppovik[k]='i';  k++; sppovik[k]='y'; k++; sppovik[k]='a'; k++; sppovik[k]='_'; k++; sppovik[k]='i'; k++; 
sppovik[k]='t'; k++; sppovik[k]='o';  k++; sppovik[k]='m'; k++; sppovik[k]='.'; k++; sppovik[k]='t'; k++; sppovik[k]='x'; k++;
sppovik[k]='t'; k++; sppovik[k]='\0'; k++; 
k=0;
skpovik[k]='K'; k++; skpovik[k]='o'; k++; skpovik[k]='e'; k++; skpovik[k]='f'; k++; skpovik[k]='f'; k++; skpovik[k]='i'; k++;
skpovik[k]='c'; k++; skpovik[k]='i'; k++; skpovik[k]='e'; k++; skpovik[k]='n'; k++; skpovik[k]='t'; k++; skpovik[k]='_'; k++;
skpovik[k]='p'; k++; skpovik[k]='o'; k++; skpovik[k]='g'; k++; skpovik[k]='l'; k++; skpovik[k]='o'; k++; skpovik[k]='s'; k++;
skpovik[k]='c'; k++; skpovik[k]='h'; k++; skpovik[k]='e'; k++; skpovik[k]='n'; k++; skpovik[k]='i'; k++; skpovik[k]='y'; k++;
skpovik[k]='a'; k++; skpovik[k]='_'; k++; skpovik[k]='i'; k++; skpovik[k]='t'; k++; skpovik[k]='o'; k++; skpovik[k]='m'; k++;
if ((vm<0) || (vm>3)) { cout << "Net takoy marki ITOM! (napstr)" << endl; k=getchar(); exit(1); }
else if (!vm) { skpovik[k]='4'; k++; skpovik[k]='4'; k++; skpovik[k]='0'; k++;  }
else if (vm==1) { skpovik[k]='6'; k++; skpovik[k]='2'; k++; skpovik[k]='0'; k++;  }
else if (vm==2) { skpovik[k]='8'; k++; skpovik[k]='6'; k++; skpovik[k]='0'; k++; }
else if (vm==3) { skpovik[k]='1'; k++; skpovik[k]='0'; k++; skpovik[k]='0'; k++; skpovik[k]='0'; k++; }
skpovik[k]='.'; k++; skpovik[k]='t'; k++; skpovik[k]='x'; k++; skpovik[k]='t'; k++; skpovik[k]='\0'; 
strcpy(sdvovi,ssdov); strcat(sdvovi,sdvovik); k=strlen(sdvovi)+1; szfov[k]='\0'; 
strcpy(sppovi,ssdov); strcat(sppovi,sppovik); k=strlen(sppovi)+1; szfov[k]='\0'; 
strcpy(skpovi,ssdov); strcat(skpovi,skpovik); k=strlen(skpovi)+1; szfov[k]='\0'; 
}
//3 - ÊÂÈ
else if (vyve==3) { sdvovkk=new char[dsov]; sppovkk=new char[dsov]; skpovkk=new char[dsov];
sdvovk=new char[2*dsov]; sppovk=new char[2*dsov]; skpovk=new char[2*dsov]; 
if ((!sdvovk) || (!sppovk) || (!skpovk) || (!sdvovkk) || (!sppovkk) || (!skpovkk)) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<(2*dsov); k++) { sdvovk[k]='\0'; sppovk[k]='\0'; skpovk[k]='\0'; } 
for (k=0; k<dsov; k++) { sdvovkk[k]='\0'; sppovkk[k]='\0'; skpovkk[k]='\0'; }
k=0;
sdvovkk[k]='D';  k++; sdvovkk[k]='l'; k++; sdvovkk[k]='i';  k++; sdvovkk[k]='n'; k++; sdvovkk[k]='y'; k++; sdvovkk[k]='_'; k++;
sdvovkk[k]='v';  k++; sdvovkk[k]='o'; k++; sdvovkk[k]='l';  k++; sdvovkk[k]='n'; k++; sdvovkk[k]='_'; k++; sdvovkk[k]='k'; k++;
sdvovkk[k]='v';  k++; sdvovkk[k]='i'; k++; sdvovkk[k]='.';  k++; sdvovkk[k]='t'; k++; sdvovkk[k]='x'; k++; sdvovkk[k]='t'; k++;
sdvovkk[k]='\0'; k++; 
k=0;
sppovkk[k]='P'; k++; sppovkk[k]='o';  k++; sppovkk[k]='k'; k++; sppovkk[k]='a'; k++; sppovkk[k]='z'; k++; sppovkk[k]='a'; k++; 
sppovkk[k]='t'; k++; sppovkk[k]='e';  k++; sppovkk[k]='l'; k++; sppovkk[k]='_'; k++; sppovkk[k]='p'; k++; sppovkk[k]='r'; k++; 
sppovkk[k]='e'; k++; sppovkk[k]='l';  k++; sppovkk[k]='o'; k++; sppovkk[k]='m'; k++; sppovkk[k]='l'; k++; sppovkk[k]='e'; k++; 
sppovkk[k]='n'; k++; sppovkk[k]='i';  k++; sppovkk[k]='y'; k++; sppovkk[k]='a'; k++; sppovkk[k]='_'; k++; sppovkk[k]='k'; k++; 
sppovkk[k]='v'; k++; sppovkk[k]='i';  k++; sppovkk[k]='.'; k++; sppovkk[k]='t'; k++; sppovkk[k]='x'; k++; sppovkk[k]='t'; k++; 
sppovkk[k]='\0'; k++; 
k=0;
skpovkk[k]='K'; k++; skpovkk[k]='o'; k++; skpovkk[k]='e'; k++; skpovkk[k]='f'; k++; skpovkk[k]='f'; k++; skpovkk[k]='i';  k++;
skpovkk[k]='c'; k++; skpovkk[k]='i'; k++; skpovkk[k]='e'; k++; skpovkk[k]='n'; k++; skpovkk[k]='t'; k++; skpovkk[k]='_';  k++;
skpovkk[k]='p'; k++; skpovkk[k]='o'; k++; skpovkk[k]='g'; k++; skpovkk[k]='l'; k++; skpovkk[k]='o'; k++; skpovkk[k]='s';  k++;
skpovkk[k]='c'; k++; skpovkk[k]='h'; k++; skpovkk[k]='e'; k++; skpovkk[k]='n'; k++; skpovkk[k]='i'; k++; skpovkk[k]='y';  k++;
skpovkk[k]='a'; k++; skpovkk[k]='_'; k++; skpovkk[k]='k'; k++; skpovkk[k]='v'; k++; skpovkk[k]='i'; k++; //skpovkk[k]='\0'; k++; 
if ((vm<4) || (vm>10)) { cout << "Net takoy marki KVI!" << endl; k=getchar(); exit(1); }
else if (vm==4) {  skpovkk[k]='4';  k++; }
else if (vm==5) {  skpovkk[k]='5';  k++; }
else if (vm==6) {  skpovkk[k]='6';  k++; }
else if (vm==7) {  skpovkk[k]='7';  k++; }
else if (vm==8) {  skpovkk[k]='8';  k++; }
else if (vm==9) {  skpovkk[k]='9';  k++; }
else if (vm==10) { skpovkk[k]='1';  k++; skpovkk[k]='0';  k++; }
skpovkk[k]='0'; k++; skpovkk[k]='0'; k++; skpovkk[k]='.'; k++; 
skpovkk[k]='t'; k++; skpovkk[k]='x'; k++; skpovkk[k]='t'; k++; skpovkk[k]='\0'; 
strcpy(sdvovk,ssdov); strcat(sdvovk,sdvovkk); k=strlen(sdvovk)+1; szfov[k]='\0'; 
strcpy(sppovk,ssdov); strcat(sppovk,sppovkk); k=strlen(sppovk)+1; szfov[k]='\0'; 
strcpy(skpovk,ssdov); strcat(skpovk,skpovkk); k=strlen(skpovk)+1; szfov[k]='\0'; 
}
}
void osvpamov(int vyve)
{ delete []sndov; delete []ssdov; delete []szfov; 
delete []szfovk; delete []snmov; delete []sfnoov; 
if (!vyve) {
delete []sppovs; delete []sppovsk; delete []skpovs; 
delete []sdvovs; delete []skpovsk; delete []sdvovsk; }
else if (vyve==1) { 
delete []sppovv; delete []sppovvk; delete []sdvovv; 
delete []sdvovvk; delete []skpovv; delete []skpovvk; }
else if (vyve==2) {
delete []sdvovi; delete []sppovi; delete []skpovi;  
delete []sdvovik; delete []sppovik; delete []skpovik; }
else if (vyve==3) {
delete []sdvovk; delete []sppovk; delete []skpovk;  
delete []sdvovkk; delete []sppovkk; delete []skpovkk; }
}
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
double epsisredver(double T, double *tdkusctm, double *dkusctm, int dkusctl, double *dkoscet, double *dkoscem, int dkoscel)
{ int k=1; napstrdir(k,k); 
double *dv=NULL; dv=dliny_voln_ver(dv, 0);
dv=new double[dlar]; 
if (!dv) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlar; k++) dv[k]=0.0;
dv=dliny_voln_ver(dv, 1);
double *npp=new double[dlar];
if (!npp) { cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<dlar; k++) npp[k]=0.0; //cout << "ce = " << dlar << endl;
npp=Kramers_Kronig_ver(dv, dlar, npp); //cout << "k = " << k << endl;
double dkusct=opredKTPTKTochSha(dkusctm, tdkusctm, T, dkusctl); 
if (dkusct>1e0) dkusct=1e0; if (dkusct<0.0) dkusct=0.0; 
double dkosce=opredKTPTKTochSha(dkoscem, dkoscet, T, dkoscel); 
if (dkosce>1e0) dkosce=1e0; if (dkosce<0.0) dkosce=0.0; 
double ume=dkusct*dkosce; 
double *epsil=new double[dlar];
if (!epsil) { cout << snmov << endl; k=getchar(); exit(1); } 
epsil=epsilnu(npp, dlar, epsil); //cout << "ume = " << ume << "\t";
double epssr=usredVelichPlank(dv, epsil, npp, T, dlar, ume); //cout << "T = " << T << "\teps_s = " << epssr << "\tume = " << ume << "\n";
if (dv) delete []dv; if (npp) delete []npp; if (epsil) delete []epsil;
k=1; osvpamov(k); return epssr; }
double epsisredsha(double T, double *tdkusctm, double *dkusctm, int dkusctl, double *dkoscet, double *dkoscem, int dkoscel)
{ int k=0; napstrdir(k,k); 
double *dv=NULL; dv=dliny_voln_sham(dv, 0);
dv=new double[dlarsh]; 
if (!dv) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlarsh; k++) dv[k]=0.0;
dv=dliny_voln_sham(dv, 1);
double *npp=new double[dlarsh];
if (!npp) { cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<dlarsh; k++) npp[k]=0.0;
npp=Kramers_Kronig_sha(dv, dlarsh, npp);
double dkusct=opredKTPTKTochSha(dkusctm, tdkusctm, T, dkusctl); 
if (dkusct>1e0) dkusct=1e0; if (dkusct<0.0) dkusct=0.0; 
double dkosce=opredKTPTKTochSha(dkoscem, dkoscet, T, dkoscel); 
if (dkosce>1e0) dkosce=1e0; if (dkosce<0.0) dkosce=0.0; 
double ume=dkusct*dkosce;
double *epsil=new double[dlarsh];
if (!epsil) { cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<dlarsh; k++) epsil[k]=0.0;
epsil=epsilnu(npp, dlarsh, epsil);
double epssr=usredVelichPlank(dv, epsil, npp, T, dlarsh, ume); //cout << "ns = " << ns << "\teps_sr = " << epssr << "\t";
delete []dv; delete []npp; delete []epsil; 
k=0; osvpamov(k); return epssr; }
double epsisreditom(double T, double *tdkusctm, double *dkusctm, int dkusctl, double *dkoscet, double *dkoscem, int dkoscel, int vi)
{ int k=2; napstrdir(k,vi);
double *dv=NULL; 
dv=dliny_voln_itom(dv, vi, 0);
dv=new double[dlari];
if (!dv) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlari; k++) dv[k]=0.0;
dv=dliny_voln_itom(dv, vi, 1);
double *npp=new double[dlari];
if (!npp) { cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<dlari; k++) npp[k]=0.0;
npp=Kramers_Kronig_itom(dv, dlari, vi, npp);
double dkusct=opredKTPTKTochSha(dkusctm, tdkusctm, T, dkusctl); //for (k=0; k<dkusctl; k++) cout << "dkusctm = " << dkusctm[k] << "\ttdkusctm = " << tdkusctm[k] << "\t"; cout << endl;
if (dkusct>1e0) dkusct=1e0; if (dkusct<0.0) dkusct=0.0; 
double dkosce=opredKTPTKTochSha(dkoscem, dkoscet, T, dkoscel); //for (k=0; k<dkoscel; k++) cout << "dkoscem = " << dkoscem[k] << "\tdkoscet = " << dkoscet[k] << "\t"; cout << endl;
if (dkosce>1e0) dkosce=1e0; if (dkosce<0.0) dkosce=0.0;
double ume=dkusct*dkosce; //cout << "dkosce = " << dkosce << "\tdkusct = " << dkusct << "\tume = " << ume << "\t";
double *epsil=new double[dlari]; 
if (!epsil) { cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<dlari; k++) epsil[k]=0.0;
epsil=epsilnu(npp, dlari, epsil);
double epssr=usredVelichPlank(dv, epsil, npp, T, dlari, ume); //cout << "T = " << T << "\teps_sr = " << epssr << "\t";
delete []dv; delete []npp; delete []epsil; 
k=2; osvpamov(k); return epssr; }
double epsisredkvi(double T, double *tdkusctm, double *dkusctm, int dkusctl, double *dkoscet, double *dkoscem, int dkoscel, int vk)
{ int k=3; napstrdir(k,vk);
double *dv=NULL; 
k=0; dv=dliny_voln_kvi(dv, k);
dv=new double[dlark];
if (!dv) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlark; k++) dv[k]=0.0;
k=1; dv=dliny_voln_kvi(dv, k);
double *npp=new double[dlark];
if (!npp) { cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<dlark; k++) npp[k]=0.0;
npp=Kramers_Kronig_kvi(dv, dlark, npp, vk); 
double dkusct=opredKTPTKTochSha(dkusctm, tdkusctm, T, dkusctl); //for (k=0; k<dkusctl; k++) cout << "dkusctm = " << dkusctm[k] << "\ttdkusctm = " << tdkusctm[k] << "\t"; cout << endl;
if (dkusct>1e0) dkusct=1e0; if (dkusct<0.0) dkusct=0.0; 
double dkosce=opredKTPTKTochSha(dkoscem, dkoscet, T, dkoscel); //for (k=0; k<dkoscel; k++) cout << "dkoscem = " << dkoscem[k] << "\tdkoscet = " << dkoscet[k] << "\t"; cout << endl;
if (dkosce>1e0) dkosce=1e0; if (dkosce<0.0) dkosce=0.0; //cout << "dkosce = " << dkosce << "\tdkusct = " << dkusct << "\t";
double ume=dkusct*dkosce;
double *epsil=new double[dlark];
if (!epsil) { cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<dlark; k++) epsil[k]=0.0;
epsil=epsilnu(npp, dlark, epsil);
double epssr=usredVelichPlank(dv, epsil, npp, T, dlark, ume); 
delete []dv; delete []npp; delete []epsil; //cout << "T = " << T << "\teps_sr = " << epssr << "\tume = " << ume << "\t";
k=3; osvpamov(k); return epssr; }
double reflsha(double T)
{ int k=0; napstrdir(k,k); 
double *dv=NULL;
dv=dliny_voln_sham(dv, 0);
dv=new double[dlarsh];
if (!dv) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlarsh; k++) dv[k]=0.0;
dv=dliny_voln_sham(dv, 1);
double *npp=new double[dlarsh];
if (!npp) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlarsh; k++) npp[k]=0.0;
npp=Kramers_Kronig_sha(dv, dlarsh, npp);
double *ronu=ronusha(npp, dlarsh), t, ume=1e0;
t=usredVelichPlank(dv, ronu, npp, T, dlarsh, ume);
delete []dv; delete []npp; delete []ronu; k=0; osvpamov(k);
return t; }
double reflver(double T)
{ int ide=0, k=1; napstrdir(k,k);
double *dv=NULL; dv=dliny_voln_ver(dv, ide);
ide=1; dv=new double[dlar]; 
if (!dv) { cout << snmov << endl; k=getchar(); exit(1); } 
dv=dliny_voln_ver(dv, ide);
double *npp=new double[dlar]; 
if (!npp) { cout << snmov << endl; k=getchar(); exit(1); } 
npp=Kramers_Kronig_ver(dv, dlar, npp);
double *ronu=ronusha(npp,dlar), t, ume=1e0;
t=usredVelichPlank(dv, ronu, npp, T, dlar, ume);
delete []dv; delete []npp; delete []ronu;
k=1; osvpamov(k); return t; }
double reflitom(double T, int vi)
{ int k=2; napstrdir(k,vi); 
double *dv=NULL;
dv=dliny_voln_itom(dv, vi, 0);
dv=new double[dlari];
if (!dv) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlari; k++) dv[k]=0.0;
dv=dliny_voln_itom(dv, vi, 1);
double *npp=new double[dlari];
if (!npp) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlari; k++) npp[k]=0.0;
npp=Kramers_Kronig_itom(dv, dlari, vi, npp);
double *ronu=ronusha(npp,dlari), t, ume=1e0;
t=usredVelichPlank(dv, ronu, npp, T, dlari, ume);
delete []dv; delete []npp; delete []ronu; k=2; osvpamov(k);
return t; }
double reflkvi(double T, int vm)
{ int vv=3, k=vv; napstrdir(k,vm);
double *dv=NULL;
k=0; dv=dliny_voln_kvi(dv, k);
dv=new double[dlark];
if (!dv) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlark; k++) dv[k]=0.0;
k=1; dv=dliny_voln_kvi(dv, k);
double *npp=new double[dlark];
if (!npp) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlark; k++) npp[k]=0.0;
npp=Kramers_Kronig_kvi(dv, dlark, npp, vm);
double *ronu=ronukvi(npp,dlark), t, ume=1e0;
t=usredVelichPlank(dv, ronu, npp, T, dlark, ume);
delete []dv; delete []npp; delete []ronu; k=vv; osvpamov(k);
return t; }
double *Kramers_Kronig_ver(double *dv, int lear, double *nnu)
{	int k=0; ofstream fo;
fo.open(sppovv,ios_base::out | ios_base::trunc); 
if (!fo.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
	double c0=299792458.0, *nu=new double[lear], *nus=new double[lear], *nnus=new double[lear];
	double *alsr=new double[lear], e=1e-20; //for (k=0; k<lear; k++) alsr[k]=0.0;
	if ((!nu) || (!nnus) || (!nus) || (!alsr)) { cout << snmov << endl; k=getchar(); exit(1); } 
	alsr=Koef_Pogl_ver(alsr); 
	for (k=0; k<lear; k++) if (dv[k]>e) nu[k]=2e0*pi*c0/dv[k]; else nu[k]=0.0; //for (k=0; k<lear; k++) cout << "k = " << k << "\tnu = " << nu[k] << endl; cout << "lear = " << lear << endl; k=getchar();
nus=izmMasChast(nu, enu, lear, nus); //for (k=0; k<lear; k++) cout << "k = " << k << "\tnu = " << nus[k] << endl; cout << "lear = " << lear << endl; k=getchar();
nnus=PokazPrelomAl(nu, nus, alsr, c0, lear, nnus); //for (k=0; k<lear; k++) cout << "k = " << k << "\tnu = " << nnus[k] << endl; cout << "lear = " << lear << endl; k=getchar();
nnu=PokazPrelomPribl(nu, nus, nnus, lear, nnu, enu); //for (k=0; k<lear; k++) cout << "k = " << k << "\tnu = " << nnu[k] << endl; cout << "lear = " << lear << endl; k=getchar();
if (nu) delete[]nu; if (nnus) delete[]nnus; if (nus) delete[]nus; if (alsr) delete[]alsr;
for (k=0; k<lear; k++) { nnu[k]=nnu[k]-1e0+1.543; fo << nnu[k] << endl; } fo.close(); 
return nnu; }
double *Kramers_Kronig_sha(double *dv, int lear, double *nnu)
{	int k; 
	ofstream fo; fo.open(sppovs,ios_base::out | ios_base::trunc); if (!fo.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
	double c0=299792458.0, *nu=new double[lear], *nus=new double[lear], *nnus=new double[lear];
	double *alsr=new double[lear];
	if ((!nu) || (!nnus) || (!nus) || (!alsr)) { cout << snmov << endl; k=getchar(); exit(1); } 
	alsr=Koef_Pogl_sha(alsr);
	for (k=0; k<lear; k++) if (dv[k]>0.0) nu[k]=2e0*pi*c0/dv[k]; else nu[k]=0.0;
nus=izmMasChast(nu, enu, lear, nus);
nnus=PokazPrelomAl(nu, nus, alsr, c0, lear, nnus);
nnu=PokazPrelomPribl(nu, nus, nnus, lear, nnu, enu); 
delete []nu; delete []nnus; delete []alsr; delete []nus;
for (k=0; k<lear; k++) { nnu[k]=nnu[k]-1e0+1.56; fo << nnu[k] << endl; } fo.close();
return nnu; }
double *Kramers_Kronig_itom(double *dv, int lear, int vi, double *nnu)
{	int k; 
	ofstream fo; fo.open(sppovi,ios_base::out | ios_base::trunc); if (!fo.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
	double c0=299792458.0, *nu=new double[lear], *nus=new double[lear], *nnus=new double[lear], ht=1e0, ko=1e-2;
	double *alsr=new double[lear], dob=0.0, w1, w2, w3, rosha=219.0*1e1, rover=25.0*1e1, phi1, phi2, v1, v2;
	alsr=Koef_Pogl_itom(alsr, vi);
	if ((!nu) || (!nus) || (!nnus) || (!alsr)) { cout << snmov << endl; k=getchar(); exit(1); } 
	for (k=0; k<lear; k++) nu[k]=2e0*pi*c0/dv[k]; 
nus=izmMasChast(nu, enu, lear, nus); 
nnus=PokazPrelomAl(nu, nus, alsr, c0, lear, nnus); 
nnu=PokazPrelomPribl(nu, nus, nnus, lear, nnu, enu);
if (!vi) { w1=6e1*ko; w2=ht-w1; }
else if (vi==1) { w1=5e1*ko; w2=ht-w1; }
else if (vi==2) { w1=25.0*ko; w2=ht-w1; }
else if (vi==3) { w1=2e1*ko; w2=ht-w1; }
else { w1=0.0; w2=0.0; }
v1=w1/rover; v2=w2/rosha; 
phi1=v1/(v1+v2); phi2=v2/(v1+v2); 
dob=(1.543-ht)*phi1+(1.56-ht)*phi2;
for (k=0; k<lear; k++) { nnu[k]=nnu[k]+dob; 
fo << nnu[k] << endl; } fo.close(); 
delete []nu; delete []nnus; delete []alsr; delete []nus;
return nnu; }
double *Kramers_Kronig_kvi(double *dv, int lear, double *nnu, int vm)
{	int k=0; ofstream fo; //cout << "vm_KK = " << vm << endl;
fo.open(sppovk,ios_base::out | ios_base::trunc); 
if (!fo.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
	double c0=299792458.0, *nu=new double[lear], *nus=new double[lear], *nnus=new double[lear];
	double *alsr=new double[lear], wv; 
	for (k=0; k<lear; k++) alsr[k]=0.0;
	if ((!nu) || (!nnus) || (!alsr) || (!nus)) { cout << snmov << endl; k=getchar(); exit(1); } 
	alsr=Koef_Pogl_kvi(alsr); 
	for (k=0; k<lear; k++) if (dv[k]>0.0) nu[k]=2e0*pi*c0/dv[k]; else nu[k]=0.0; 
nus=izmMasChast(nu, enu, lear, nus);
nnus=PokazPrelomAl(nu, nus, alsr, c0, lear, nnus); 
nnu=PokazPrelomPribl(nu, nus, nnus, lear, nnu, enu);
if (vm==4) wv=49e-2; else if (vm==5) wv=37e-2; else if (vm==6) wv=31e-2; else if (vm==7) wv=27e-2;
else if (vm==8) wv=23e-2; else if (vm==9) wv=22e-2; else if (vm==10) wv=2e-1; 
else { cout << "Oshibka v vybore marki KVI!"; k=getchar(); exit(1); }
wv=(1.543-1e0)*wv+(1e0-wv)*(1.56-1e0); 
for (k=0; k<lear; k++) { nnu[k]=nnu[k]+wv; fo << nnu[k] << endl; } fo.close();
if (nu) delete []nu; if (nnus) delete []nnus; if (nus) delete []nus; if (alsr) delete []alsr; 
return nnu; }
double *ronusha(double *npp, int lear)
{ int k; double *ronu=new double[lear], n; if (!ronu) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<lear; k++) {
	n=fabs(npp[k]);
    ronu[k]=0.5+((n-1.0)*(3.0*n+1.0))/(6.0*pow(n+1.0,2.0))-(2.0*pow(n,3.0)*(pow(n,2.0)+2.0*n-1.0))/((pow(n,2.0)+1.0)*(pow(n,4.0)-1.0)); 
    ronu[k]=ronu[k]+(8.0*pow(n,4.0)*(pow(n,4.0)+1.0)*log(n))/((pow(n,2.0)+1.0)*pow((pow(n,4.0)-1.0),2.0));
    ronu[k]=ronu[k]+pow(n,2.0)*pow((pow(n,2.0)-1.0),2.0)*log(fabs((n-1.0)/(n+1.0)))/pow((pow(n,2.0)+1.0),3.0); }
return ronu; }
double *ronukvi(double *npp, int lear)
{ int k; double *ronu=new double[lear], n; if (!ronu) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<lear; k++) {
	n=fabs(npp[k]);
    ronu[k]=0.5+((n-1.0)*(3.0*n+1.0))/(6.0*pow(n+1.0,2.0))-(2.0*pow(n,3.0)*(pow(n,2.0)+2.0*n-1.0))/((pow(n,2.0)+1.0)*(pow(n,4.0)-1.0)); 
    ronu[k]=ronu[k]+(8.0*pow(n,4.0)*(pow(n,4.0)+1.0)*log(n))/((pow(n,2.0)+1.0)*pow((pow(n,4.0)-1.0),2.0));
    ronu[k]=ronu[k]+pow(n,2.0)*pow((pow(n,2.0)-1.0),2.0)*log(fabs((n-1.0)/(n+1.0)))/pow((pow(n,2.0)+1.0),3.0); }
return ronu; }
double *epsilnu(double *npp, int lear, double *epsarr)
{ int k; double eps, n;
for (k=0; k<lear; k++) {
	n=fabs(npp[k]);
	eps=(4.0*n+2.0)/3.0/pow((n+1.0),2.0);
	eps=eps+2.0*pow(n,3.0)*(pow(n,2.0)+2.0*n-1.0)/(pow(n,2.0)+1.0)/(pow(n,4.0)-1.0);
	eps=eps-8.0*pow(n,4.0)*(pow(n,4.0)+1.0)*log(n)/(pow(n,2.0)+1.0)/pow((pow(n,4.0)-1.0),2.0);
	eps=eps-pow(n,2.0)*log((n-1.0)/(n+1.0))*pow((pow(n,2.0)-1.0),2.0)/pow((pow(n,2.0)+1.0),3.0); 
	epsarr[k]=eps; }
return epsarr; }
double *PokazPrelomAl(double *nu, double *nus, double *ka, double vl, int arle, double *np)
{ int k, j, p=arle-1, q; double dkpdo, podln, *fn=new double[arle]; 
if (!fn) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<arle; k++) {
    q=0; for (j=0; j<p; j++) {
        dkpdo=(ka[j+1]-ka[j])/(nu[j+1]-nu[j]);
        podln=((nu[j]+nus[k])/(nu[j]-nus[k]));
        podln=fabs(podln);
        fn[q]=dkpdo*log(podln);
        q++; }
    fn[p]=fn[p-1];
    np[k]=1.0+(vl/pi)*trapz(nu,fn,arle)/2.0/nu[k]; }
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
double *Koef_Pogl_ver(double *kp)
{	ifstream fin; double p; int k=0, q=100; char *s=new char[q]; 
	if (!s) { cout << snmov << endl; k=getchar(); exit(1); }
fin.open(skpovv,ios_base::in); for (k=0; k<q; k++) s[k]='\0'; //cout << "k = " << k << endl;
if (!fin.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
k=0; while ((!fin.eof()) && (k<dlar)) { fin.getline(s,q,'\n'); p=atof(s); kp[k]=p; k++; }
fin.close(); delete []s; return kp; }
double *Koef_Pogl_sha(double *kp)
{	ifstream fin; int lear=dlarsh; double p; 
	int k, q=100; char *s=new char[q]; 
	if (!s) { cout << snmov << endl; k=getchar(); exit(1); } 
fin.open(skpovs,ios_base::in); for (k=0; k<q; k++) s[k]='\0';
if (!fin.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
k=0; while ((!fin.eof()) && (k<dlarsh)) { fin.getline(s,q,'\n'); p=atof(s); kp[k]=p; k++; }
fin.close(); delete []s; return kp; }
double *Koef_Pogl_itom(double *kp, int vi)
{	ifstream fin; int lear=dlari; double p; 
	int k, q=100; char *s=new char[q]; 
	if (!s) { cout << snmov << endl; k=getchar(); exit(1); } 
	fin.open(skpovi, ios_base::in); 
for (k=0; k<q; k++) s[k]='\0';
if (!fin.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
k=0; while ((!fin.eof()) && (k<dlari)) { fin.getline(s,q,'\n'); p=atof(s); kp[k]=p; k++; }
fin.close(); delete []s; return kp; }
double *Koef_Pogl_kvi(double *kp)
{	ifstream fin; double p; int k=0, q=100; char *s=new char[q]; if (!s) { cout << snmov << endl; k=getchar(); exit(1); }
fin.open(skpovk,ios_base::in); for (k=0; k<q; k++) s[k]='\0'; //cout << "k_KP = " << k << endl;
if (!fin.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
k=0; while ((!fin.eof()) && (k<dlark)) { fin.getline(s,q,'\n'); p=atof(s); kp[k]=p; k++; }
fin.close(); delete []s; //cout << skpovk << "\t"; for (k=0; k<dlark/50; k++) cout << kp[k] << "\t"; 
return kp; }
double *dliny_voln_ver(double *dv, int ide)
{ ifstream fin; double p; int k, q=100, lear; char *s=new char[q]; 
if (!s) { cout << snmov << endl; k=getchar(); exit(1); } for (k=0; k<q; k++) s[k]='\0';
fin.open(sdvovv,ios_base::in); 
if (!fin.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
k=0; while (!fin.eof()) { fin.getline(s,q,'\n'); p=atof(s); k++; } 
fin.clear(); fin.seekg(0);
lear=k; dlar=lear; 
if (!ide) { delete []s; fin.close(); return NULL; } else {
k=0; while ((!fin.eof()) && (k<dlar)) { fin.getline(s,q,'\n'); p=atof(s); dv[k]=p;  k++; }
fin.close();
for (k=0; k<lear; k++) dv[k]=(1e-2)/dv[k];
delete []s; return dv; } }
double *dliny_voln_sham(double *dv, int ide)
{ ifstream fin; double p; int k, q=100, lear; char *s=new char[q]; 
if (!s) { cout << snmov << endl; k=getchar(); exit(1); } for (k=0; k<q; k++) s[k]='\0';
fin.open(sdvovs,ios_base::in); 
if (!fin.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); } 
k=0; while (!fin.eof()) { fin.getline(s,q,'\n'); p=atof(s); k++; } 
fin.clear(); fin.seekg(0);
lear=k; dlarsh=lear; if (!ide) { fin.close(); delete []s; return NULL; } else {
k=0; while (!fin.eof()) { fin.getline(s,q,'\n'); p=atof(s); dv[k]=p;  k++; }
fin.close();
for (k=0; k<lear; k++) dv[k]=1e-2/dv[k]; //for (k=0; k<lear; k++) cout << dv[k] << "\t"; 
delete []s; return dv; } }
double *dliny_voln_itom(double *dv, int vi, int ide)
{ ifstream fin; double p; int k, q=100, lear; char *s=new char[q]; 
if (!s) { cout << snmov << endl; k=getchar(); exit(1); } for (k=0; k<q; k++) s[k]='\0';
fin.open(sdvovi,ios_base::in);
if (!fin.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); } 
k=0; while (!fin.eof()) { fin.getline(s,q,'\n'); p=atof(s); k++; } 
fin.clear(); fin.seekg(0);
lear=k; dlari=lear; if (!ide) { fin.close(); delete []s; return NULL; } else {
k=0; while (!fin.eof()) { fin.getline(s,q,'\n'); p=atof(s); dv[k]=p;  k++; }
fin.close();
for (k=0; k<lear; k++) dv[k]=1e-2/dv[k]; //for (k=0; k<lear; k++) cout << dv[k] << "\t"; 
delete []s; return dv; } }
double *dliny_voln_kvi(double *dv, int ide)
{ ifstream fin; double p; int k, q=100, lear; char *s=new char[q]; 
if (!s) { cout << snmov << endl; k=getchar(); exit(1); } for (k=0; k<q; k++) s[k]='\0';
fin.open(sdvovk,ios_base::in); 
if (!fin.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); } 
k=0; while (!fin.eof()) { fin.getline(s,q,'\n'); p=atof(s); k++; } 
lear=k; dlark=lear; fin.clear(); k=0; fin.seekg(k);
if (!ide) { fin.close(); delete []s; return NULL; } else if (ide==1) {
k=0; while (!fin.eof()) { fin.getline(s,q,'\n'); p=atof(s); dv[k]=p;  k++; }
fin.close();
for (k=0; k<lear; k++) dv[k]=1e-2/dv[k]; 
delete []s; return dv; } }
void zapisvfile(double *vyvod, int dlina, char *nazf)
{ 	if (dlina>0) { int k; FILE *fo = fopen(nazf, "a+"); 
	if (!fo) { cout << sfnoov << endl; k=getchar(); exit(1); } cout << "dl_zap" << dlina << endl;
	for (k=0; k<dlina; k++) fprintf(fo,"%0.15le\n",vyvod[k]); //fprintf(fo,"%\n\n\n"); 
	fclose(fo); } }
double trapz(double *arx, double *ary, int lenarr)
{ int k; double s=0.0;
for (k=1; k<lenarr; k++) s=s+(ary[k]+ary[k-1])*(arx[k]-arx[k-1])/2e0;
return s; }
double *oslaintestepcher(int dm, int vy, double *sc) //ðàññ÷èòûâàåò îñëàáëåíèå èíòåãðàëüíîé ñòåïåíè ÷åðíîòû äëÿ SiO2, Al2O3, MgO
{	int n=2, k, j, w=0, p=6, q=10, m, l; double **scv=new double*[p], **hsv=new double*[p];
	double *hs=new double[n], **kor=new double*[q], t0, t1, *vscs=new double[dm], *vscm=new double[dm], *vsca=new double[dm], *po=NULL;
if ((!scv) || (!sc) || (!hsv) || (!hs) || (!kor)) { cout << "No memory" << endl; k=getchar(); exit(1);}
	k=0;
	sc[k]=1e0;   k++; sc[k]=0.94;  k++; sc[k]=0.87;  k++; sc[k]=0.801; k++; sc[k]=0.736; k++; 
	sc[k]=0.676; k++; sc[k]=0.635; k++; sc[k]=0.590; k++; sc[k]=0.567; k++; sc[k]=0.543; k++; 
	sc[k]=0.53;  k++; sc[k]=0.525; k++; sc[k]=0.515; k++; sc[k]=0.507; //ñòåïåíü ÷åðíîòû ìàãíåçèòà
	w=0; scv[w]=sc; sc=new double[dm];
	k=0;
	hs[k]=98e-2; k++; hs[k]=2e-2;
	k=0;
	t0=hs[k]; k++; t1=hs[k];
	k=0;
	hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1);
	hsv[w]=hs; w++; hs=new double[n]; //Ìàãíåçèò: 98 % - MgO, 2 % - SiO2
	if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1); }
k=0;
sc[k]=1e0;   k++; sc[k]=0.976; k++; sc[k]=0.949; k++; sc[k]=0.905; k++; sc[k]=0.859; k++; 
sc[k]=0.812; k++; sc[k]=0.774; k++; sc[k]=0.737; k++; sc[k]=0.709; k++; sc[k]=0.681; k++; 
sc[k]=0.661; k++; sc[k]=0.639; k++; sc[k]=0.626; k++; sc[k]=0.620; //ñòåïåíü ÷åðíîòû øàìîòà
scv[w]=sc; sc=new double[dm]; 
k=0; 
hs[k]=56e-2; k++; hs[k]=396e-3; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Øàìîò: 56 % - SiO2, 39,6 % - Al2O3 - øàìîò, ñ÷èòàåì, ÷òî ìàãíèÿ íåò â øàìîòå
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1);}
k=0;
sc[k]=1e0;    k++; sc[k]=98e-2;  k++; sc[k]=951e-3; k++; sc[k]=92e-2;  k++; sc[k]=883e-3; k++; 
sc[k]=853e-3; k++; sc[k]=821e-3; k++; sc[k]=79e-2;  k++; sc[k]=767e-3; k++; sc[k]=746e-3; k++;
sc[k]=73e-2;  k++; sc[k]=715e-3; k++; sc[k]=705e-3; k++; sc[k]=692e-3; //ñòåïåíü ÷åðíîòû êîðóíäîøàìîòà
scv[w]=sc; sc=new double[dm]; 
k=0; 
hs[k]=28e-2; k++; hs[k]=7e-1; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Êîðóíäîøàìîò: 28 % - SiO2, 70 % - Al2O3 - êîðóíäîøàìîò
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1);}
k=0;
sc[k]=1.0000; k++; sc[k]=983e-3; k++; sc[k]=936e-3; k++; sc[k]=867e-3; k++; sc[k]=819e-3; k++; 
sc[k]=721e-3; k++; sc[k]=659e-3; k++; sc[k]=593e-3; k++; sc[k]=541e-3; k++; sc[k]=49e-2;  k++; 
sc[k]=453e-3; k++; sc[k]=429e-3; k++; sc[k]=403e-3; k++; sc[k]=384e-3; //ñòåïåíü ÷åðíîòû êàîëèíîâîãî òåïëîèçîëÿöèîííîãî êèðïè÷à (ÊÒÊ)
scv[w]=sc; sc=new double[dm];
k=0; 
hs[k]=57e-2; k++; hs[k]=4e-1; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //ÊÒÊ: 57 % - SiO2, 40 % - Al2O3
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1); }
k=0;
sc[k]=1.000;  k++; sc[k]=984e-3; k++; sc[k]=941e-3; k++; sc[k]=882e-3; k++; sc[k]=813e-3; k++; 
sc[k]=751e-3; k++; sc[k]=695e-3; k++; sc[k]=641e-3; k++; sc[k]=594e-3; k++; sc[k]=558e-3; k++; 
sc[k]=53e-2;  k++; sc[k]=499e-3; k++; sc[k]=479e-3; k++; sc[k]=462e-3; //ñòåïåíü ÷åðíîòû ìóëëèòà
scv[w]=sc; sc=new double[dm];
k=0; 
hs[k]=28e-2; k++; hs[k]=72e-2; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Ìóëëèò: 28 % - SiO2, 72 % - Al2O3
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1);}
k=0;
sc[k]=1e0;    k++; sc[k]=984e-3; k++; sc[k]=953e-3; k++; sc[k]=917e-3; k++; sc[k]=854e-3; k++; 
sc[k]=808e-3; k++; sc[k]=756e-3; k++; sc[k]=711e-3; k++; sc[k]=578e-3; k++; sc[k]=523e-3; k++; 
sc[k]=495e-3; k++; sc[k]=468e-3; k++; sc[k]=448e-3; k++; sc[k]=429e-3; //ñòåïåíü ÷åðíîòû êðåìíåçåìà
scv[w]=sc; 
k=0; 
hs[k]=985e-3; k++; hs[k]=1e-2; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Êðåìíåçåì: 98,5 % - SiO2, 1 % - Al2O3
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
{	int k=0, j=0, dm=0, vm=vi; napstrdir(vyve, vm); //cout << "vyve = " << vyve << "vm = " << vm << endl; 
	double *dl=NULL, *alfs=NULL, *npp=NULL, *kpr=new double[arle], krpk, temk;
	double *ktr=new double[arle], dkusct=1e0, dkosce=1e0; //ÊÒÏ ïî Ðîññåëàíäó
	double dko2=dkosups, dko4=dkokpiko, dko0=1e0; //ó÷åò ïîðèñòîé ñòðóêòóðû è ñîîòíîøåíèÿ ÊÏ/ÊÎ
	if (!ident) { //äëÿ øàìîòà
			k=0; dl=dliny_voln_sham(dl, k); dl=new double[dlarsh]; if (!dl) { cout << snmov << endl; k=getchar(); exit(1); } 
				for (k=0; k<dlarsh; k++) dl[k]=0.0; k=1; dl=dliny_voln_sham(dl, k); 
				alfs=new double[dlarsh]; if (!alfs) { cout << snmov << endl; k=getchar(); exit(1); } for (k=0; k<dlarsh; k++) alfs[k]=0.0;
				alfs=Koef_Pogl_sha(alfs); npp=new double[dlarsh]; if (!npp) { cout << snmov << endl; k=getchar(); exit(1); } 
				for (k=0; k<dlarsh; k++) npp[k]=0.0; npp=Kramers_Kronig_sha(dl, dlarsh, npp); dm=dlarsh; }
	else if (ident==1) { //äëÿ âåðìèêóëèòà
			k=0; dl=dliny_voln_ver(dl, k); dl=new double[dlar]; if (!dl) { cout << snmov << endl; k=getchar(); exit(1); } 
				for (k=0; k<dlar; k++) dl[k]=0.0; k=1; dl=dliny_voln_ver(dl, k); 
				alfs=new double[dlar]; if (!alfs) { cout << snmov << endl; k=getchar(); exit(1); } for (k=0; k<dlar; k++) alfs[k]=0.0;
				alfs=Koef_Pogl_ver(alfs); npp=new double[dlar]; if (!npp) { cout << snmov << endl; k=getchar(); exit(1); } 
				for (k=0; k<dlar; k++) npp[k]=0.0; npp=Kramers_Kronig_ver(dl, dlar, npp); dm=dlar; } 
	else if (ident==2) { //äëÿ ÈÒÎÌ
				k=0; dl=dliny_voln_itom(dl, vm, k); dl=new double[dlari]; if (!dl) { cout << snmov << endl; k=getchar(); exit(1); } 
				for (k=0; k<dlari; k++) dl[k]=0.0; k=1; dl=dliny_voln_itom(dl, vm, k); 
				alfs=new double[dlari]; if (!alfs) { cout << snmov << endl; k=getchar(); exit(1); } for (k=0; k<dlari; k++) alfs[k]=0.0; 
				alfs=Koef_Pogl_itom(alfs, vm); npp=new double[dlari]; if (!npp) { cout << snmov << endl; k=getchar(); exit(1); } 
				for (k=0; k<dlari; k++) npp[k]=0.0; npp=Kramers_Kronig_itom(dl, dlari, vm, npp); dm=dlari; }
	else if (ident==3) { //äëÿ ÊÂÈ
				k=0; dl=dliny_voln_kvi(dl, k); dl=new double[dlark]; if (!dl) { cout << snmov << endl; k=getchar(); exit(1); } 
				for (k=0; k<dlark; k++) dl[k]=0.0; k=1; dl=dliny_voln_kvi(dl, k); 
				alfs=new double[dlark]; if (!alfs) { cout << snmov << endl; k=getchar(); exit(1); } 
				for (k=0; k<dlark; k++) alfs[k]=0.0; alfs=Koef_Pogl_kvi(alfs); 
				npp=new double[dlark]; if (!npp) { cout << snmov << endl; k=getchar(); exit(1); } 
				for (k=0; k<dlark; k++) npp[k]=0.0; npp=Kramers_Kronig_kvi(dl, dlark, npp, vm); dm=dlark; }
	for (k=0; k<dm; k++) alfs[k]=alfs[k]*dko2*dko4; //cout << "dlark = " << dlark << endl; //for (k=0; k<(dm/50); k++) cout << alfs[k] << "\t"; for (k=0; k<(dm/50); k++) cout << npp[k] << "\t"; 
	double *GnpT=new double[arle], sigma=5.67e-8, *npT=new double[arle]; //ïðîèçâîäíàÿ dn/dT, ôóíêöèÿ <n>(T)
	if ((!kpr) || (!GnpT) || (!ktr) || (!npT)) { cout << snmov << endl; k=getchar(); exit(1); }
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
{ double rpa=(n2*cos(phi)-n1*cos(phis))/(n2*cos(phi)+n1*cos(phis)); //êîýôôèöèåíò îòðàæåíèÿ ïàðàëëåëüíûé ïëîñêîñòè ïàäåíèÿ
double rpe=(n1*cos(phi)-n2*cos(phis))/(n1*cos(phi)+n2*cos(phis)); //êîýôôèöèåíò îòðàæåíèÿ ïåðïåíäèêóëÿðíûé ê ïëîñêîñòè ïàäåíèÿ
return (fabs(pow(rpa,2e0))+fabs(pow(rpe,2e0)))/2e0; }
double ReflSredSha(double Te)
{ int k=0; napstrdir(k,k); 
double *dv=NULL;
dv=dliny_voln_sham(dv, 0);
dv=new double[dlarsh];
dv=dliny_voln_sham(dv, 1);
if (!dv) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlarsh; k++) dv[k]=0.0;
double *npp;
npp=new double[dlarsh];
if (!npp) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlarsh; k++) npp[k]=0.0;
npp=Kramers_Kronig_sha(dv, dlarsh, npp);
double *ronu=PoiSpeKoeOtr(npp,dlarsh), ume=1e0; 
double t=usredVelichPlank(dv, ronu, npp, Te, dlarsh, ume); 
delete []dv; delete []npp; delete []ronu; k=0; osvpamov(k); 
return t; }
double ReflSredVer(double T)
{ int ide=0, k=1; napstrdir(k,k);
double *dv=NULL; dv=dliny_voln_ver(dv, ide);
ide=1; dv=new double[dlar]; 
if (!dv) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlar; k++) dv[k]=0.0;
dv=dliny_voln_ver(dv, ide); 
double *npp=new double[dlar]; 
if (!npp) { cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<dlar; k++) npp[k]=0.0;
npp=Kramers_Kronig_ver(dv, dlar, npp);
double *ronu=PoiSpeKoeOtr(npp,dlar), ume=1e0;
double t=usredVelichPlank(dv, ronu, npp, T, dlar, ume);
delete []dv; delete []npp; delete []ronu; k=1; osvpamov(k);
return t; }
double ReflSreditom(double Te, int vi)
{ int k=2; napstrdir(k,vi);
double *dv=NULL;
dv=dliny_voln_itom(dv, vi, 0);
dv=new double[dlari];
if (!dv) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlari; k++) dv[k]=0.0;
dv=dliny_voln_itom(dv, vi, 1); 
double *npp=new double[dlari]; 
if (!npp) { cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<dlari; k++) npp[k]=0.0;
npp=Kramers_Kronig_itom(dv, dlari, vi, npp);
double *ronu=PoiSpeKoeOtr(npp,dlari), ume=1e0;
double t=usredVelichPlank(dv, ronu, npp, Te, dlari, ume); cout << "Refl = " << t << "\t";
delete []dv; delete []npp; delete []ronu; k=2; osvpamov(k); 
return t; }
double ReflSredkvi(double Te, int vm)
{ int k=3; napstrdir(k,vm);
double *dv=NULL; 
k=0; dv=dliny_voln_kvi(dv, k);
dv=new double[dlark]; 
if (!dv) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlark; k++) dv[k]=0.0;
k=1; dv=dliny_voln_kvi(dv, k);
double *npp=new double[dlark]; 
if (!npp) { cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<dlark; k++) npp[k]=0.0;
npp=Kramers_Kronig_kvi(dv, dlark, npp, vm);
double *ronu=PoiSpeKoeOtr(npp,dlark), ume=1e0;
double t=usredVelichPlank(dv, ronu, npp, Te, dlark, ume);
delete []dv; delete []npp; delete []ronu; k=3; osvpamov(k); 
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
double SeryeStenkiRasIzl(double ww, double hh, double ll, double *T, double *EPS, double *HOs, double *RHOs, double *A, int *ns, int chst)
{ N=chst; double b, t, epsi=1e-10; t=T[ns[0]]; b=T[ns[1]]; if (fabs(t-b)<epsi) return 0.0;
double **Fs=new double*[N], z, r, y; int i, j, f=3; 
double *PIN=new double[N], *POUT=new double[N], *arg=new double[f], *p;
for (i=0; i<N; i++) { p=new double[N]; for (j=0; j<N; j++) p[j]=0.0; Fs[i]=p; }
	double sumq=0.0, sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	int *id=new int[N], iclsd; for (i=0; i<N; i++) id[i]=1; // Ïîâåðõíîñòè 1 è 3 (bottom and top)
    	for (i=0; i<N; i++) if (!(id[i])) PIN[i]=q[i]; else PIN[i]=sigma*pow(T[i],4e0); //Çàïîëíåíèå ìàññèâà PIN ÏÒÏ q è òåìïåðàòóðîé T - ïåðåâîä òåìïåðàòóð â ÏÏÈ, Fs(1,2)=F1-2, 
	iclsd=1; //óãëîâûå êîýôôèöèåíòû; êîíôèãóðàöèÿ çàìêíóòà (iclsd=1), äèàãîíàëüíûå ýëåìåíòû íå íóæíû //for (j=0; j<N; j++) cout << "RHOs ( " << j << " ) = " << RHOs[j] << endl;  for (j=0; j<N; j++) cout << "EPS ( " << j << " ) = " << EPS[j] << endl;
	arg[0]=fabs(hh); arg[1]=fabs(ll); arg[2]=fabs(ww); //h,l,w
	Fs[0][1]=VIEW(39,f,arg); //Fs(1,2)=F1-2, F12
    	arg[0]=fabs(ww); arg[1]=fabs(ll); arg[2]=fabs(hh);     // a,b,c
		r=VIEW(38,f,arg); t=parlpltf(ww,ww,2.0*ww,ll,0.0,ll,hh); b=RHOs[1];
    	Fs[0][2]=r+b*t; //Fs1-3=F1-3+rhos2*F1(2)-3, Fs13
    	Fs[0][3]=Fs[0][1]+RHOs[1]*perppltf(ww,2.0*ww,0.0,hh,ll,0.0,ll); //Fs1-4=F1-4+rhos2*F1(2)-4
    	Fs[0][4]=0.5*(1.0-(1.0-RHOs[1])*Fs[0][1]-Fs[0][2]-Fs[0][3]); //ñìåùåíû äâà ïðÿìîóãîëüíèêà íà ww ïî îñè X
    	Fs[0][5]=Fs[0][4]; //Fs23
    	arg[0]=fabs(ww); arg[1]=fabs(ll); arg[2]=fabs(hh); // h,l,w
    	Fs[1][2]=VIEW(39,3,arg)+RHOs[0]*perppltf(hh,2.0*hh,0.0,ww,ll,0.0,ll);  // Fs2-3=F2-3+rhos1*F1(1)-3 //ñìåùåíû äâà ïðÿìîóãîëüíèêà íà hh ïî îñè X, ò.ê. ó íèõ ðàçíûå ðàçìåðû
	   	arg[0]=fabs(ll); arg[1]=fabs(hh); arg[2]=fabs(ww); //a,b,c
    	Fs[1][3]=VIEW(38,3,arg)+RHOs[0]*parlpltf(hh,hh,2.0*hh,ll,0.0,ll,ww);  // Fs2-4=F2-4+rhos1*F2(1)-4, Fs24
    	Fs[1][4]=0.5*(1.0-(1.0-RHOs[0])*A[0]*Fs[0][1]/A[1]-Fs[1][2]-Fs[1][3]); // Fs25
    	Fs[1][5]=Fs[1][4];
    	arg[0]=fabs(ww); arg[1]=fabs(ll); arg[2]=2.0*fabs(hh);  // a,b,c
    	Fs[2][2]=RHOs[0]*(VIEW(38,3,arg)+RHOs[1]*parlpltf(ww,ww,2.0*ww,ll,0.0,ll,2.0*hh)); //Fs33=rhos1*F3(1)-3 + rhos1*rhos2*F3(12+21)-3 //èçëó÷åíèå ñòåíêè ñàìó íà ñåáÿ
    	Fs[2][3]=Fs[0][3]+RHOs[0]*(perppltf(0.0,ww,hh,2.0*hh,ll,0.0,ll)+RHOs[1]*perppltf(ww,2.0*ww,hh,2.0*hh,ll,0.0,ll)); // Fs3-4=Fs1-4+rhos1*F3(1)-3+rhos1*rhos2*F3(12+21)-4
    	Fs[2][4]=0.5*(1.0-(1.0-RHOs[0])*A[0]*Fs[0][2]/A[2]-(1.0-RHOs[1])*A[1]*Fs[1][2]/A[2]-Fs[2][2]-Fs[2][3]);
    	Fs[2][5]=Fs[2][4];
    	arg[0]=fabs(hh); arg[1]=fabs(ll); arg[2]=2.0*fabs(ww); // a,b,c
    	Fs[3][3]=RHOs[1]*(VIEW(38,3,arg)+RHOs[0]*parlpltf(hh,hh,2.0*hh,ll,0.0,ll,2.0*ww)); //Fs44=rhos2*F4(2)-4+rhos1*rhos2*F4(12+21)-4
    	Fs[3][4]=0.5*(1.0-((1.0-RHOs[0])*A[0]*Fs[0][3]+(1.0-RHOs[1])*A[1]*Fs[1][3]+A[2]*Fs[2][3])/A[3]-Fs[3][3]);
    	Fs[3][5]=Fs[3][4];
    	Fs[4][5]=0.5*(1.0-((1.0-RHOs[0])*A[0]*Fs[0][4]+(1.0-RHOs[1])*A[1]*Fs[1][4]+A[2]*Fs[2][4]+A[3]*Fs[3][4])/A[4]); //Fs56 //for (j=0; j<N; j++) cout << "A ( " << j << " ) = " << A[j] << endl;  //for (i=0; i<N; i++) { y=RHOs[i]+EPS[i]; cout << "Edinitsa = " << y << endl; } //for (i=0; i<N; i++) cout << "Temp ( " << i << " ) = " << T[i] << endl; 
		POUT=GRAYDIFFSPEC(iclsd, A, EPS, RHOs, HOs, Fs, id, PIN, POUT, N); //Ðåøåíèå ñèñòåìû óðàâíåíèé ìåòîäîì Ãàóññà //for (i=0; i<N; i++) cout << POUT[i] << endl; //POUT=GRAYDIFF(iclsd, A, EPS, HOs, Fs, id, PIN, POUT, N);
	sumq=0.0; //ïðîâåðêà ïîëíîé ÏÏÈ
    	for (i=0; i<N; i++) //Âûâîä - ïåðåâîä â òåìïåðàòóðû
		{ if (!(id[i])) T[i]=pow((POUT[i]/sigma), 1e0/4e0); else q[i]=POUT[i];
        	QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; //cout << "Surface " << i << "\tT [K] = " << T[i] << "\tq [W/m2] = " << q[i] << "\tQA [W] = " << QA[i]/ll << endl;
		} r=0.0; z=1.0; for (i=0; i<N; i++) if (ns[i]==i) { r=r+z*q[i]; z=-z; } 
		for (i=0; i<N; i++) { p=Fs[i]; delete []p; } delete []Fs; delete []PIN; 
		delete []POUT; delete []arg; delete []q; delete []QA; delete []id; return r; }
double *GRAYDIFF(int iclsd, double *A, double *EPS, double *HOs, double **Fs, int *ID, double *PIN, double *POUT, int chst)
{
N=chst; int i, j;
double **qm=new double*[N], **em=new double*[N], **pm=new double*[N], *B=new double[N], *p, *idf=new double[N], hf=1e0, s, ikr, epsi=1e-10;
for (i=0; i<N; i++) { p=new double[N]; qm[i]=p; p=new double[N]; em[i]=p; p=new double[N]; pm[i]=p; }
for (i=0; i<N; i++) { s=0.0; for (j=0; j<ID[i]; j++) s=s+hf; idf[i]=s; }
for (i=1; i<N; i++)
	 for (j=0; j<i; j++)
		  if (fabs(A[i])>epsi) 
			  Fs[i][j]=(A[j]/A[i])*Fs[j][i]; else Fs[i][j]=0.0; //Äëÿ çàêðûòîé êîíôèãóðàöèè, íåîáõîäèìî ðàññ÷èòàòü äèàãîíàëüíûå ýëåìåíòû èç ïðàâèëà ñóììèðîâàíèÿ
	if (iclsd==1) { for (i=0; i<N; i++) {
		Fs[i][i]=1e0;
		for (j=0; j<N; j++) {
			if (j!=i) {
				Fs[i][i]=Fs[i][i]-Fs[i][j]; } } } }
	 for (i=0; i<N; i++) { // Çàïîëíåíèå ìàòðèö êîýôôèöèåíòîâ q è e
		for (j=0; j<N; j++) {
			if (i==j) ikr=1e0; else ikr=0.0; // Ñèìâîë Êðîíåêåðà delta_ij //cout << "i = " << i << "\tj = " << j << "\tk = " << ikr << endl;
			if (fabs(EPS[j])>epsi) 
				qm[i][j]=ikr/EPS[j]-(1e0/EPS[j]-1e0)*Fs[i][j];
			em[i][j]=ikr-Fs[i][j]; } }
	 for (i=0; i<N; i++) { // Çàïîëíåíèå ìàòðèöû âûõîäíûõ êîýôôèöèåíòîâ POUT è RHS
			B[i]=-HOs[i];
			for (j=0; j<N; j++) {
					pm[i][j]=qm[i][j]*idf[j]-em[i][j]*(1e0-idf[j]);
					B[i]=B[i]+(em[i][j]*idf[j]-qm[i][j]*(1e0-idf[j]))*PIN[j]; } }
	 POUT=MetodGaussa(pm, B, N, POUT); //for (i=0; i<N; i++) cout << "pout = " << POUT[i] << endl; cout << endl; POUT=GAUSS(pm, B, N, POUT); for (i=0; i<N; i++) cout << "pout = " << POUT[i] << endl; cout << endl;
	 double sumq = 0.0, *T=new double[N], sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	 for (i=0; i<N; i++) { if (!(ID[i])) T[i]=pow((POUT[i]/sigma),1e0/4e0); else q[i]=POUT[i];
     QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; } //printf("Summa vsekh potokov = %0.10lf\n", sumq); 
	 delete []T; delete []q; delete []QA;
	 for (i=0; i<N; i++) { p=qm[i]; delete []p; p=em[i]; delete []p; p=pm[i]; delete []p; } 
	 delete []qm; delete []em; delete []pm; delete []B; delete []idf; return POUT;
}
double *GRAYDIFFSPEC(int iclsd, double *A, double *EPS, double *RHOs, double *HOs, double **Fs, int *ID, double *PIN, double *POUT, int chst)
{ N=chst; int i, j; 
double **qm=new double*[N], **em=new double*[N], **pm=new double*[N], *B=new double[N], *p, *idf=new double[N], hf=1e0, s, ikr, y, epsi=1e-10; // Âû÷èñëåíèå íåäîñòàþùèõ óãëîâûõ êîýôôèöèåíòîâ - ëåâîé íèæíåé ÷àñòè ìàòðèöû èç ñâîéñòâà âçàèìíîñòè
for (i=0; i<N; i++) { p=new double[N]; qm[i]=p; p=new double[N]; em[i]=p; p=new double[N]; pm[i]=p; }
for (i=0; i<N; i++) { s=0.0; for (j=0; j<ID[i]; j++) s=s+hf; idf[i]=s; }
for (i=1; i<N; i++)
	 for (j=0; j<i; j++)
		  if (fabs(A[i])>epsi) 
			  Fs[i][j]=(A[j]/A[i])*Fs[j][i]; else Fs[i][j]=0.0; //Äëÿ çàêðûòîé êîíôèãóðàöèè, íåîáõîäèìî ðàññ÷èòàòü äèàãîíàëüíûå ýëåìåíòû èç ïðàâèëà ñóììèðîâàíèÿ
	if (iclsd==1) { for (i=0; i<N; i++) {
		Fs[i][i]=1e0;
		for (j=0; j<N; j++) {
			if (j!=i) {
				Fs[i][i]=Fs[i][i]-(1e0-RHOs[j])*Fs[i][j]; } }
		if (RHOs[i]<1e0) Fs[i][i]=Fs[i][i]/(1e0-RHOs[i]); else RHOs[i]=0.0; } }
	 for (i=0; i<N; i++) { // Çàïîëíåíèå ìàòðèö êîýôôèöèåíòîâ q è e
		for (j=0; j<N; j++) {
			if (i==j) ikr=1e0; else ikr=0.0; // Ñèìâîë Êðîíåêåðà delta_ij //cout << "i = " << i << "\tj = " << j << "\tk = " << ikr << endl;
			if (fabs(EPS[j])>epsi) 
				qm[i][j]=ikr/EPS[j]-((1.0-RHOs[j])/EPS[j]-1.0)*Fs[i][j];
			em[i][j]=ikr-(1.0-RHOs[j])*Fs[i][j]; } }
	 for (i=0; i<N; i++) { // Çàïîëíåíèå ìàòðèöû âûõîäíûõ êîýôôèöèåíòîâ POUT è RHS
			B[i]=-HOs[i];
			for (j=0; j<N; j++) {
					pm[i][j]=qm[i][j]*idf[j]-em[i][j]*(1e0-idf[j]);
					B[i]=B[i]+(em[i][j]*idf[j]-qm[i][j]*(1e0-idf[j]))*PIN[j]; } } //cout << "Fs" << endl; for (i=0; i<N; i++) { for (j=0; j<N; j++) printf("%0.4lf\t",Fs[i][j]); cout << endl; } //for (i=0; i<N; i++) { y=0.0; for (j=0; j<N; j++) y=y+(1e0-RHOs[j])*Fs[i][j]; cout << "stroka ( " << i << " ) = " << y << endl; } //j=getchar();
	 POUT=MetodGaussa(pm, B, N, POUT); //for (i=0; i<N; i++) cout << "pout = " << POUT[i] << endl; cout << endl; POUT=GAUSS(pm, B, N, POUT); for (i=0; i<N; i++) cout << "pout = " << POUT[i] << endl; cout << endl; POUT=reshMetObrMatr(pm, B, N, POUT); for (i=0; i<N; i++) cout << "pout = " << POUT[i] << endl; cout << endl;
	 double sumq = 0.0, *T=new double[N], sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	for (i=0; i<N; i++) { if (!(ID[i])) T[i]=pow((POUT[i]/sigma),1e0/4e0); else q[i]=POUT[i];
    QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; } //printf("Summa = %0.10lf\n", sumq); 
	delete []T; delete []q; delete []QA;
	for (i=0; i<N; i++) { p=qm[i]; delete []p; p=em[i]; delete []p; p=pm[i]; delete []p; } 
	 delete []qm; delete []em; delete []pm; delete []B; delete []idf; return POUT; } //Èíâåðòèðîâàíèå ìàòðèöû êîýôôèöèåíòîâ POUT è óìíîæåíèå íà ìàòðèöó RHS, ÷òîáû ïîëó÷èòü POUT
double *SEMIGRAY(int iclsd, double *A, double **EPS, double **RHOs, double **HOs, double ***Fs, int *ID, double *q, double *T, double *PIN, double *POUT, int chst, int chel)
{ N=chst; int i, j, k, *id1=new int[N];
double ikr, sigma=5.670E-8, *epsl=new double[N], *rhosl=new double[N], **FSl=new double*[N], *q1=new double[N], *HOs2=new double[N], *p; // Ðàñ÷åò ïëîòíîñòè òåïëîâîãî ïîòîêà  äëÿ âíåøíåãî èçëó÷åíèÿ (äëÿ êàæîãî äèàïàçîíà)
for (j=0; j<N; j++) { p=new double[N]; FSl[j]=p; }
    for (i=0; i<N; i++) { id1[i]=1; PIN[i]=0.0; } // Óñòàíîâèòü çíà÷åíèÿ eps, rhos è Fs äëÿ k-îãî äèàïàçîíà
	k=0; for (i=0; i<N; i++) { epsl[i]=EPS[k][i]; rhosl[i]=RHOs[k][i];
        for (j=i; j<N; j++) FSl[i][j]=Fs[k][i][j]; }
	GRAYDIFFSPEC(iclsd, A, epsl, rhosl, HOs[k], FSl, id1, PIN, q1, N); //Âû÷èñëåíèå íåèçâåñòíûõ ïëîòíîñòåé ïîòîêîâ òåïëîâîé ýíåðãèè è òåìïåðàòóð ïðè q2=q-q1
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
double VIEW(int NO, int NARG, double *ARG)
{ double VIEW=0.0, A, B, C, X, Y, RTX, RTY, RT, H, L, W, HH, WW, W2, H2, HW2, HW, H12, W12, C1, C2, C3, epsi=1e-10, e=epsi; int k;
if ((NO<38) || (NO>39)) { cout << "Illegal value for number (NO =" << NO << ") " << endl; k=getchar(); exit(1); }
if (NO==38) { if (NARG!=3) {cout << "Wrong number of input parameters (NARG =" << NARG << ") for NO = " << NO << endl; k=getchar(); exit(1);}
      A=ARG[0]; B=ARG[1]; C=ARG[2]; 
      X=A/C; Y=B/C;
      if ((X<e) || (Y<e)) VIEW=0.0; else {
      RTX=pow(1.0+X*X,1.0/2.0); RTY=pow(1.0+Y*Y,1.0/2.0); RT=pow(1.0+X*X+Y*Y,1.0/2.0);
	  VIEW=(log(RTX*RTY/RT)+X*RTY*atan(X/RTY)+Y*RTX*atan(Y/RTX)-X*atan(X)-Y*atan(Y))*2.0/(pi*X*Y); } }
else if (NO==39) { if(NARG!=3) { cout << "Wrong number of input parameters (NARG =" << NARG << ") for NO = " << NO << endl; k=getchar(); exit(1); }
      H=ARG[0]; L=ARG[1]; W=ARG[2];
      HH=H/L; WW=W/L; W2=WW*WW; H2=HH*HH; HW2=H2+W2; HW=sqrt(H2+W2); H12=H2+1.0; W12=W2+1.0;
      C1=W12*H12/(H12+W2); C2=W2*(H12+W2)/W12/HW2; C3=H2*(H12+W2)/H12/HW2;
      VIEW=(WW*atan(1.0/WW)+HH*atan(1.0/HH)-HW*atan(1.0/HW)+0.25*(log(C1)+W2*log(C2)+H2*log(C3)))/(pi*WW); } 
return VIEW; }
double perppltf(double X1, double X2, double Y1, double Y2, double Z1, double Z2, double Z3)
{ double perppltf, A=0.0, F=0.0, *ARG=new double[3], epsi=1e-10, e=epsi;
int NARG=3;
      ARG[0]=Y2; ARG[1]=Z3; ARG[2]=X2;
      A=X2*Z3;
      if (fabs(A*Y2)>e) F=A*VIEW(39,NARG,ARG); else F=0.0;
      ARG[0]=Y1; A=X2*Z3;
      if (fabs(A*Y1)>e) F=F-A*VIEW(39,NARG,ARG);
      ARG[0]=Y2; ARG[2]=X1; A=X1*Z3;
      if (fabs(A*Y2)>e) F=F-A*VIEW(39,NARG,ARG);
      ARG[0]=Y1; A=X1*Z3;
      if (fabs(A*Y1)>e) F=F+A*VIEW(39,NARG,ARG);
      ARG[0]=Y2; ARG[1]=Z2; A=X1*Z2;
      if (fabs(A*Y2)>e) F=F+A*VIEW(39,NARG,ARG);
      ARG[0]=Y1; A=X1*Z2;
      if (fabs(A*Y1)>e) F=F-A*VIEW(39,NARG,ARG);
      ARG[0]=Y2; ARG[2]=X2; A=X2*Z2;
      if (fabs(A*Y2)>e) F=F-A*VIEW(39,NARG,ARG);
      ARG[0]=Y1; A=X2*Z2;
      if (fabs(A*Y1)>e) F=F+A*VIEW(39,NARG,ARG);
      ARG[0]=Y2; ARG[1]=(Z3-Z1); A=X2*(Z3-Z1);
      if (fabs(A*Y2)>e) F=F-A*VIEW(39,NARG,ARG);
      ARG[0]=Y1; A=X2*(Z3-Z1);
      if (fabs(A*Y1)>e) F=F+A*VIEW(39,NARG,ARG);
      ARG[0]=Y2; ARG[2]=X1; A=X1*(Z3-Z1);
      if (fabs(A*Y2)>e) F=F+A*VIEW(39,NARG,ARG);
      ARG[0]=Y1; A=X1*(Z3-Z1);
      if (fabs(A*Y1)>e) F=F-A*VIEW(39,NARG,ARG);
      ARG[0]=Y2; ARG[1]=(Z2-Z1); ARG[2]=X2; A=X2*(Z2-Z1);
      if (fabs(A*Y2)>e) F=F+A*VIEW(39,NARG,ARG);
      ARG[0]=Y1; A=X2*(Z2-Z1); 
	  if (fabs(A*Y1)>e) F=F-A*VIEW(39,NARG,ARG);
      ARG[0]=Y2; ARG[2]=X1; A=X1*(Z2-Z1);
      if(fabs(A*Y2)>e) F=F-A*VIEW(39,NARG,ARG);
      ARG[0]=Y1; A=X1*(Z2-Z1);
      if(fabs(A*Y1)>e) F=F+A*VIEW(39,NARG,ARG);
      if ((fabs(Z1)>e) && (X2!=X1)) perppltf=F/(2.0*(X2-X1)*Z1); else perppltf=0.0;
	  delete []ARG; return perppltf; }
double parlpltf(double X1, double X2, double X3, double Y1, double Y2, double Y3, double C)
{ int NARG=3; double parlpltf, A=0.0, F=0.0, *ARG=new double[NARG], epsi=1e-10, e=epsi, t;
      ARG[0]=X3; ARG[1]=Y3; ARG[2]=C; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=A*VIEW(38,NARG,ARG); else F=0.0;
      ARG[1]=Y2; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG);
      ARG[1]=Y3-Y1; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG);
	  ARG[1]=Y2-Y1; A=ARG[0]*ARG[1]; t=VIEW(38,NARG,ARG);
      if(fabs(A)>e) F=F+A*VIEW(38,NARG,ARG);
      ARG[0]=X2; ARG[1]=Y3; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG);
      ARG[1]=Y2; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F+A*VIEW(38,NARG,ARG);
      ARG[1]=Y3-Y1; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F+A*VIEW(38,NARG,ARG);
      ARG[1]=Y2-Y1; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG);
      ARG[0]=X3-X1; ARG[1]=Y3; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG);
      ARG[1]=Y2; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F+A*VIEW(38,NARG,ARG);
      ARG[1]=Y3-Y1; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F+A*VIEW(38,NARG,ARG);
      ARG[1]=Y2-Y1; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG);
      ARG[0]=X2-X1; ARG[1]=Y3; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F+A*VIEW(38,NARG,ARG);
      ARG[1]=Y2; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG);
      ARG[1]=Y3-Y1; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG);
      ARG[1]=Y2-Y1; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F+A*VIEW(38,NARG,ARG);
      if ((fabs(X1)>e) && (fabs(Y1)>e)) parlpltf=F/(4.0*(X1*Y1));
      delete []ARG; return parlpltf; }
double *GAUSS(double **A, double *B, int chst, double *X)
{ N=chst; int I, *L=new int[N], K, J, LK;
double *S=new double[N], SMAX, RMAX, R, XMULT, SUM, epsi=1e-10;
if ((!S) || (!L)) { cout << "No memory!" << endl; K=getchar(); exit(1); }
      for (I=0; I<N; I++)
			{ L[I]=I; SMAX=0.0; //L[I] - íîìåð ñòðîêè
			for (J=0; J<N; J++)
				if (fabs(A[I][J])>SMAX) 
					SMAX=fabs(A[I][J]);
			S[I]=SMAX; } //S[I] - ìàññèâ ìàêñèìàëüíûõ ýëåìåíòîâ â ñòðîêå
      for (K=0; K<N-1; K++) 
        { RMAX=0.0; J=K;
        for (I=K; I<N; I++) 
		{ if (fabs(S[L[I]])>epsi) 
				R=fabs(A[L[I]][K])/S[L[I]]; else R=0.0; //äåëèì ïåðâûå ýëåìåíòû íà S[I]
		if (R>RMAX) { J=I; RMAX=R; } } //íàõîäèì 'ýëåìåíò ñ ìàêñèìàëüíûì îòíîøåíèåì K-ãî ýëåìåíòà ê ìàêñèìàëüíîìó â K-îé ñòðîêå, çàïèñûâàåì åãî íîìåð
		LK=L[J]; L[J]=L[K]; L[K]=LK; //ìåíÿåì ìåñòàìè L[K] - òåêóùèé è L[J] - ñ ìàêñèìàëüíûì
		for (I=K+1; I<N; I++)
			{ XMULT = A[L[I]][K]/A[LK][K]; //äåëèì êàæäûé ýëåìåíò K-îãî ñòîëáöà íà ìàêñèìàëüíûé, íà÷èíàÿ ñ (K+1)-ãî ýëåìåíòà ñòðîêè 
					if (A[LK][K]==0.0)
						XMULT=0.0;
				for (J=K+1; J<N; J++)
					A[L[I]][J]=A[L[I]][J]-XMULT*A[LK][J]; //âû÷èòàåì, ÷òîáû K-ûé ýëåìåíò â K-îì ñòîëáöå îáíóëèëñÿ
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
double F0_lamT(double lT0)
{ int k=20, j=0; double *F0=new double[k], *lamT=new double[k], dolya;
if ((!F0) || (!lamT)) { cout << "No memory!" << endl; j=getchar(); exit(1); }
lamT[j]=0.18889e-2; F0[j]=0.05059; j++; lamT[j]=0.22222e-2; F0[j]=0.10503; j++; lamT[j]=0.24444e-2; F0[j]=0.14953; j++;
lamT[j]=0.27222e-2; F0[j]=0.21033; j++; lamT[j]=0.28889e-2; F0[j]=0.24803; j++; lamT[j]=0.31667e-2; F0[j]=0.31067; j++;
lamT[j]=0.33889e-2; F0[j]=0.35933; j++; lamT[j]=0.36111e-2; F0[j]=0.40585; j++; lamT[j]=0.38889e-2; F0[j]=0.46031; j++;
lamT[j]=0.41111e-2; F0[j]=0.50066; j++; lamT[j]=0.44444e-2; F0[j]=0.55573; j++; lamT[j]=0.47778e-2; F0[j]=0.60449; j++;
lamT[j]=0.51667e-2; F0[j]=0.65402; j++; lamT[j]=0.56111e-2; F0[j]=0.70211; j++; lamT[j]=0.61667e-2; F0[j]=0.75146; j++;
lamT[j]=0.68889e-2; F0[j]=0.80152; j++; lamT[j]=0.78889e-2; F0[j]=0.85171; j++; lamT[j]=0.93889e-2; F0[j]=0.90031; j++;
lamT[j]=1.25556e-2; F0[j]=0.95094; j++; lamT[j]=2.33333e-2; F0[j]=0.99051; j++;
dolya=opredKTPTKTochSha(F0, lamT, lT0, k); if (dolya<0.0) dolya=0.0; if (dolya>=1.0) dolya=1.0; delete []F0; delete []lamT; return dolya; }
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
