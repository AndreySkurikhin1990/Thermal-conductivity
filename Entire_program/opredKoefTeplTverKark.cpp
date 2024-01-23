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
//-----
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
double *DulnevZern(double, double *, double *, double *, double *, double, double, double, double, int, int, 
	double *, double *, double *, double *, double *, int, int);
double DulnevSloForSha(double, double, double, double, double, double, double);
double MoleSostTeplVozd(double, double, double, double);
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
double urovPodder(double, double, double);
double sign(double);
double DulKoefTep1CylPer(double, double, double, double);
double opredKoefAkkomodLandau(double);
double opredPlotnVozd(double);
double *arrTempAir();
double *arrPlotnAir();
double opredDulnLam1(double, double, double, double, double, double, int, double);
double *MetodDulnevaSigalovoyPolyDispVer(double *, int, double, double *, int, double, double *, double *, double *, double *, int, double);
double *MetodDulnevaSigalovoyDopVer(double *, int, double, double *, int, double, double *, double *, double *, double *, int, double);
double *MetodDulnevaSigalovoyBezIspVer(double *, int, double, double *, int, double, double *, double *, double *, double *, int, double);
double DulnevSloForVer(double, double, double, double, double, double, double, int);
double opredUrovPodderM03(double);
double *VasFraySha(double *, int);
double SysPlastShakhPor(double, double, double);
double MetodStarostina(double, double, double);
double MetodRusselNeprVozd(double, double, double);
double MetodRusselNeprTverTelo(double, double, double);
double *VasFrayObsh(double *, double *, double, int);
double Poiskkx(double);
double *DopFunctOpredKTPTverChas(double *, double *, int, double, double *, int, int, double *, double, double *, 
	double *, double *, double *, int, double);
double **vybFunRasPorPoRazSha(double, int, int, char *, int);
double opredKTPTKToch(double *, double *, double, int);
double urovPod(double);
double *proverkakvi(double *, double *, double *, double, double, int);
double **RasPorPoRazkvi(int, char *);
double **RasPorPoRazitom(int, char *);
double **RasPorPoRazVer(double, int, int, int, int, char *, int);
double **opredSredKTPTK();
double **NapMasVozd(char *);
double *opredKTPTverKarkNach(double *, double *, double, double, double, double, int, int, int, int, double *, double *, int, 
	double *, int, int, double *, double *, double *, int, char *, int, int, int, int);
double *oprKTPTverKark(double *, double *, double, double, double, double, int, int, double *, double *tkuscv, double *stchv, double *srra, double *raspr, int dmv, double srp, double *ktpvotk, 
	double *vtetk, double *Prvotk, int dmkv, char *snm, int vfv, int vysove, int vpkf, int vybves, int vybmar, int fl);
double oprSrRazPor(int, int, double, int, int, int, int, int, int, char *);
double oprMaxRazPor(int, int, double, int, int, int, int, int, int, char *);
double *proverkaKTPTKNaAdek(double *, double *, int);
double **napMasEKTPVer(int, int, int, double *, int, double, int, int, int, double **, int, int, char *);
char **napolStrok(int, int);
double **pereraschetTverKark(double **);
//-----//Управляющие функции
double *opredKTPTverKarkNach(double *tem, double *laefm, double por, double wsio, double walo, double wmgo, int vfv, int vysove, 
	int vpkf, int isrp, double *kusc, double *tkusc, int dmkusc, double *stch, int lentem, int fl, double *ktpvotk, double *vtetk, 
	double *Prvoz, int dmkvoz, char *snm, int vybves, int vrsh, int vystsha, int vybmar)
{
	int k=0, dm=0, l=k+1, m=l+1, n=m+1, q=n+1, v=q+1; 
	double **mu=NULL, *ktptk=NULL, r=0.0, ko=1e6;
	if (!vybves) mu=vybFunRasPorPoRazSha(por, vrsh, vystsha, snm, vybves);
	if (vybves==1) mu=RasPorPoRazVer(por, vfv, vysove, isrp, vpkf, snm, vybves); //0 - старые, 1 - новые значения
	if (vybves==2) { 
		if ((vybmar>=0) && (vybmar<=3)) mu=RasPorPoRazitom(vybmar, snm); 
		else { cout << "Marka ITOM ne naydena!" << endl; k=getchar(); exit(1); } }
	if (vybves==3) { 
		if ((vybmar>=4) && (vybmar<=10)) mu=RasPorPoRazkvi(vybmar, snm); 
		else { cout << "Net takoy marki keramovermikulita!" << endl; k=getchar(); exit(1); } }
	mu=pereraschetTverKark(mu);
	double *raspr=mu[k], *srra=mu[l], hf=1e0, e=1e-6, t=0.0;
	double *prgr=mu[m], *legr=mu[n], *mez=mu[q], srp=0.0, marp=0.0;
k=0; srp=mez[k]; k++; marp=mez[k]; r=marp; if (r<hf) r=r*ko; 
t=e; k=0; while (t<r) { k++; t=t+hf; } dm=k; 
cout << "dm = " << dm << "\tsrp = " << srp << "\tmrp = " << marp << endl;
ktptk=oprKTPTverKark(tem, laefm, por, wsio, walo, wmgo, lentem, dmkusc, kusc, tkusc, stch, srra, 
raspr, dm, srp, ktpvotk, vtetk, Prvoz, dmkvoz, snm, vfv, vysove, vpkf, vybves, vybmar, fl);
if (raspr) delete[]raspr; if (srra) delete[]srra; if (legr) delete[]legr; 
if (prgr) delete[]prgr; if (mez) delete[]mez; if (mu) delete[]mu; 
return ktptk; }
double *oprKTPTverKark(double *tem, double *laefm, double por, double wsio, double walo, double wmgo, int n, int dmkusc, 
	double *kuscv, double *tkuscv, double *stchv, double *srra, double *raspr, int dmv, double srp, double *ktpvotk, 
	double *vtetk, double *Prvotk, int dmkv, char *snm, int vfv, int vysove, int vpkf, int vybves, int vybmar, int fl)
{ 
	int k=0, cvp=22, cve=cvp+7*2, cved=cve+1; if (vybves==1) cved=cve+4;
	int na=0, j=0, q=0, *no=new int[cved], f=(cve-cvp)/2, w=0;
	double *srk=NULL, **sr=new double*[cved], *lavo=new double[n], rmi=0.0, rma=0.0, ep=1e-8, x0=0.0;
	double *srzn=new double[n], s=0.0, p=0.0, srch=0.0, r=0.0, ht=1e0, lvmi=ht, ko=1e-3, pormax=9e-1, ot=ht/3e0;
if ((!no) || (!sr) || (!lavo) || (!srzn)) { cout << snm << endl; j=getchar(); exit(1); } 
for (k=0; k<cved; k++) { no[k]=0; sr[k]=NULL; }
for (k=0; k<n; k++) { lavo[k]=opredKTPTKToch(ktpvotk, vtetk, tem[k], dmkv); //cout << "te = " << tem[k] << "\tlavo = " << lavo[k] << "\t";
if (lvmi>lavo[k]) lvmi=lavo[k]; }
if (vybves==1) {
if (!vfv) { rma=2e0; rmi=7e-1; }
else if (vfv==1) { rma=8e0; rmi=4e0; }
else if (vfv==2) { rma=1.6; rmi=0.35; } 
srch=(rmi+rma)*ko/2e0; }
else srch=fabs(srp)/pow(por,ot)/2e0;
w=0; p=0.0; for (k=na; k<cved; k++) {
srk=new double[n]; if (!srk) { cout << snm << endl; j=getchar(); exit(1); } for (j=0; j<n; j++) srk[j]=0.0;
	if (k==na) srk=DulnevZern(por, tem, srk, laefm, lavo, wsio, walo, wmgo, srch, n, dmkusc, kuscv, tkuscv, stchv, Prvotk, 
		vtetk, dmkv, vybves);
	if ((k<cvp) && (k>na)) srk=opredTvChaKoeTepSrInSpSha(k, por, tem, laefm, n, srk, lavo, srch);
    if ((k<cve) && (k>=cvp)) srk=DulnevKoefTepVerm(k, por, tem, n, srk, laefm, lavo, srp, f, srra, raspr, vybves, vfv, dmv, cvp); 
	if (!vybves) { if ((k<cved) && (k>=cve)) srk=DopFunctOpredKTPTverChas(srk, tem, n, por, laefm, k-cve, vfv, lavo, srch, stchv, Prvotk, 
		ktpvotk, vtetk, dmkv, ep); }
	if (vybves>=1) if ((k<cved) && (k>=cve)) srk=VasFrayObsh(laefm, lavo, por, n);
	srk=proverkaKTPTKNaAdek(srk, laefm, n); //for (j=0; j<n; j++) cout << "j = " << j << "\tk = " << k << "\tsrk = " << srk[j] << endl; cout << endl;
	q=1; for (j=0; j<n; j++) if (srk[j]<lvmi) { q=-1; break; } 
	if (q>0) { sr[w]=srk; no[w]=k; w++; p=p+ht; } else { if (srk) delete[]srk; } }
	cved=w; for (j=0; j<n; j++) { 
		s=0.0; for (k=0; k<cved; k++) { w=no[k]; r=sr[k][j]; s=s+r; } 
		if (s>lvmi) srzn[j]=s/p; } 
	if (fl<0) { cout << "cved = " << cved << "\tl = " << lvmi << "\tsrp = " << srp << endl;
	for (k=0; k<cved; k++) { w=no[k]; srk=sr[k]; if (srk) { for (j=0; j<n; j++) cout << "k = " << w+1 << "\ttem = " << tem[j] << "\tlam_tk = " << srk[j] << endl; cout << endl; } } 
	cout  << "Nomera nenulevyh elementov" << endl; for (k=0; k<cved; k++) cout << no[k]+1 << " "; cout << endl; 
	for (k=0; k<n; k++) cout << "T ( " << k << " ) = " << tem[k] << "\tlam_tk = " << srzn[k] << "\tlam_e = " << laefm[k] << endl; }
for (k=0; k<cved; k++) { srk=sr[k]; if (srk) delete []srk; } if (sr) delete []sr; 
if (lavo) delete []lavo; if (no) delete []no;
return srzn; }
//Функции расчета КТП твердого каркаса и вспомогательные для них
double *proverkaKTPTKNaAdek(double *srk, double *laefm, int n)
{
int f=1, k=0;
double e=1e-8, r=0.0, t=0.0, p=0.0;
for (k=0; k<n; k++) {
    r=srk[k]; t=laefm[k];
    if (r<t) { f=-1; break; }
    if (k>0) { p=srk[k-1];
    if (fabs(p-r)<e) { f=-1; break; }
    if (p>r) { f=-1; break; } } }
if (f<0) for (k=0; k<n; k++) srk[k]=0.0;
return srk;
}
double opredFiBolN(double x, double y)
{ 
	int n=11, k=0; double *y1=new double[n], *y01=new double[n], *fib=new double[n], fb01, fb1, fibo; 
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
if (y1) delete[]y1; if (y01) delete[]y01; if (fib) delete []fib;
return fibo; }
double DulnevSloForSha(double m2, double T, double eps, double lavo, double lam1, double srch, double Prvoz)
{ 
	double r=srch, d=2.0*r, epf=1e-2, hf=1e0, ot=hf/3e0; 
	double Nk=pow(pow(m2,2.0)-10.0*m2+9.0,0.5); Nk=(m2+3.0+Nk)/2.0/m2; //y2=3.3e-3*pow(1-m2,-2./9.); pud=0; hsl=30e-3; rona=0.2; ro1=rona/(1-m2); y2=y2*pow(pud+9.8*ro1*(1-m2)*hsl,1./3.); y1=1e-20; y2=y1;
	double eta=1e-4, y1=1e-3, y2=y1/pow(eta,0.5); //y1=(10.0+50.0)*1e-4/2.0;
	double y3=2.0*pow(Nk-1.0,0.5)/Nk, y4=y3/pow(hf-m2,ot); //hsh=2e-3/2.0;
	double hsh=1e-10, fb=0.0, fbn=fb, fbnn=fb;
fbn=opredFiBolN(y2,y2/y3); fbnn=opredFiBolNN(lavo/lam1,m2);
if (fbn>hf) if (fbnn<hf) fb=fbnn; else fb=0.0; else if (fbnn>hf) fb=fbn; else 
if (fabs(fbn-fbnn)<epf) fb=(fbn+fbnn)/2.0; else if (fbn>fbnn) fb=fbnn; else fb=fbn;
double A=pow(y2,2.0)-pow(y1,2.0), F=pow(hf-pow(y2,2.0),0.5);
double D=pow(hf-pow(y3,2.0),0.5), E=pow(y4,2.0)-pow(y3,2.0);
double delsrsz=d*(hsh+hf/Nk), lamszm=MoleSostTeplVozd(T, delsrsz, lavo, Prvoz);
double epspr=eps/(2.0-eps), sig=5.668e-8;
double lamszl=4.0*epspr*sig*delsrsz*pow(T,3.0), lamsz=lamszm+lamszl;
double nusz=lamsz/lam1, w=(lavo/lamszm-nusz*D)/(hf-nusz);
double y12=pow(y1,2.0)/(hsh/2.0+(hf-hsh/2.0)*fb), DF=fabs(w-D)/fabs(w-F);
double nug=nusz, numz=MoleSostTeplVozd(T, hsh*r, lavo, Prvoz); numz=numz/lam1;
double nu2sp=nug; DF=(D-F+w*log(DF))*2.0*nug/(hf-nug); 
double AF=A/(hf-hsh/2.0-F+hsh/2.0/numz), ADF=hf/(AF+DF);
ADF=hf/(D/pow(y3,2.0)+ADF); ADF=ADF+nu2sp*E+y12;
ADF=ADF*lam1/pow(y4,2.0);
return ADF; }
double opredFiBolNN(double nu, double m2)
{ 
	double tocras=1e-8, hf=1e0, cb=1e2, ca=-1e1, ra=hf, ep=tocras, cc=hf, fa=hf, fb=hf, fc=hf; 
	int h=0, kit=100;
while ((ra>ep) && (h<kit)) {
cc=(ca+cb)/2.0;
fa=2.0*pow(ca,3.0)-3.0*pow(ca,2.0)+hf-m2;
fb=2.0*pow(cb,3.0)-3.0*pow(cb,2.0)+hf-m2;
fc=2.0*pow(cc,3.0)-3.0*pow(cc,2.0)+hf-m2;
if ((fc*fb>0.0) && (fa*fc<0.0)) cb=cc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) ca=cc; 
ra=fabs(fa-fb); h++; }
double fi1=fabs(cc/(2.0-cc)); fi1=pow(fi1,0.5);
fi1=pow(nu,2.0)*(7.5-11.0*nu+4.5*pow(nu,2.0))*(hf-fi1)+fi1;
return fi1; }
double MoleSostTeplVozd(double T, double de, double lamg, double Pr)
{ 
	double pi=acos(-1e0), gam=7e0/5e0, H=101325.0, H0=1e5, kB=1.380658e-23, hf=1e0, angs=1e-10;
	double n=H/kB/T, ko=1e-2, wo2=21.0*ko, do2=6e-1, wn2=79.0*ko, dn2=65e-2;
	double d=(do2*wo2+dn2*wn2)*angs, sig=pi*pow(d,2.0)/4.0;
	double dli=hf/pow(2e0,0.5)/n/sig;
	double a=opredKoefAkkomodLandau(T), cz=8.42e-3, Ty=113.0, la0=cz/H0/(hf+Ty/T);
	double Kn=la0*H0/H/de, B=4e0*(gam/(gam+hf))*((2e0-a)/a)*(Kn/Pr);
	lamg=lamg/(hf+B);
return lamg; }
double urovPodder(double pro, double ref, double urpo)
{ 
	double fl=0.0; int f=1; 
	if (pro<0) f=0; else if (pro<(urpo*ref)) f=0; else f=1;
if (f) fl=pro; else fl=0.0; 
return fl; }
double urovOtsechen(double pro, double ref, double urot)
{ 
	double fl=0.0; int f=1; 
	if (pro<0.0) f=-1; else if (pro>(urot*ref)) f=-1; else f=1;
if (f>0) fl=pro; else fl=0.0; 
return fl; }
double ReshNeyaUravMaxKoTepl(double lae, double lan, double po) //lan - КТП непрерывной фазы (воздух), ищем КТП твердого скелета
{ 
	double lgl=1e-9, pgl=1e5, hf=1e0, ladb=pgl, lada=lgl, ladc=0.0, ra=pgl;
	double tocras=1e-8, ep=tocras, fa=hf, fb=hf, fc=hf; 
	int h=0, hko=1000; //1
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
{ 
	double lgl=1e-9, pgl=1e5, hf=1e0, ladb=pgl, lada=lgl, ladc=hf, ra=pgl;
	double tocras=1e-8, ep=tocras, fa=hf, fb=hf, fc=hf; 
	int h=0, hko=1000; //2
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
{ 
	double tocras=1e-8, lgl=1e-9, pgl=1e5, hf=1e0, ladb=pgl, lada=lgl, ladc=1e5, ra=pgl, ep=tocras;
	double fa=hf, fb=hf, fc=hf; 
	int h=0, hko=1000; //3
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2.0;
fa=pow(lada,(hf-po))*pow(lan,po)-lae;
fb=pow(ladb,(hf-po))*pow(lan,po)-lae;
fc=pow(ladc,(hf-po))*pow(lan,po)-lae;
if ((fc*fb>0) && (fa*fc<0)) ladb=ladc; 
if ((fc*fa>0) && (fb*fc<0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Lie = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplMaxwell(double lae, double lan, double po)
{ 
	double lgl=1e-9, pgl=1e5, hf=1e0, ladb=pgl, lada=lgl, ladc=hf;
	double tocras=1e-8, ra=pgl, ep=tocras, fa=hf, fb=hf, fc=hf; 
	int h=0, hko=1000; //4
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0;
fa=lan*(lada+2e0*lan-2e0*(hf-po)*(lan-lada))/(lada+2e0*lan+(hf-po)*(lan-lada))-lae;
fb=lan*(ladb+2e0*lan-2e0*(hf-po)*(lan-ladb))/(ladb+2e0*lan+(hf-po)*(lan-ladb))-lae;
fc=lan*(ladc+2e0*lan-2e0*(hf-po)*(lan-ladc))/(ladc+2e0*lan+(hf-po)*(lan-ladc))-lae;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Maxw = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplBurgFrik(double lae, double lan, double po)
{ 
	int h=0, k=3, j=1, w=2, hko=1000;
	double lgl=1e-9, pgl=1e5, hf=1e0, ladb=pgl, lada=lgl, ladc=hf, ra=pgl;
	double tocras=1e-8, ep=tocras, fa=hf, fb=hf, fc=hf, *f=new double[k];
	double fba=hf, fbb=hf, fbc=hf, zn=-hf; 
if (!f) { cout << "No memory" << endl; k=getchar(); exit(1); } 
k=0; f[k]=hf/8e0; f[j]=f[k]; f[w]=hf-f[k]-f[j]; //5
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0;
fba=0.0; fbb=0.0; fbc=0.0;
for (k=0; k<3; k++) {
     fba=fba+pow(hf+((lada/lan-hf)*f[k]),zn);
     fbb=fbb+pow(hf+((ladb/lan-hf)*f[k]),zn);
     fbc=fbc+pow(hf+((ladc/lan-hf)*f[k]),zn); }
fba=fba/3e0; fbb=fbb/3e0; fbc=fbc/3e0;
fa=lan*(hf+(hf-po)*(fba*lada/lan-hf))/(hf+(hf-po)*(fba-hf))-lae;
fb=lan*(hf+(hf-po)*(fbb*ladb/lan-hf))/(hf+(hf-po)*(fbb-hf))-lae;
fc=lan*(hf+(hf-po)*(fbc*ladc/lan-hf))/(hf+(hf-po)*(fbc-hf))-lae;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } 
if (f) delete[]f; //cout << "lam_Bur_Fri = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplShumanVoss(double lae, double lan, double po)
{ 
	double lgl=1e-9, pgl=1e5, hf=1e0, ladb=pgl, lada=lgl, ladc=1e5, ra=pgl;
	double tocras=1e-8, ep=tocras, fa=hf, fb=hf, fc=hf; //6
	double lamaa=hf, lamab=hf, lamac=hf, p=hf; 
	int h=0, hko=1000;
p=PoiskPShumanVoss(po); //cout << "Shum Voss p = " << p << endl;
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2.0;
lamaa=lan*lada/(lan+p*(lan-lada));
lamaa=lamaa*(hf+p*(hf+p)*(lan-lada)*log(lan*(hf+p)/p/lada)/(lan+p*(lan-lada)));
fa=lan*pow(po,3.0)+(hf-pow(po,3.0))*lamaa-lae;
lamab=lan*ladb/(lan+p*(lan-ladb));
lamab=lamab*(hf+p*(hf+p)*(lan-ladb)*log(lan*(hf+p)/p/ladb)/(lan+p*(lan-ladb)));
fb=lan*pow(po,3.0)+(hf-pow(po,3.0))*lamab-lae;
lamac=lan*ladc/(lan+p*(lan-ladc));
lamac=lamac*(hf+p*(hf+p)*(lan-ladc)*log(lan*(hf+p)/p/ladc)/(lan+p*(lan-ladc)));
fc=lan*pow(po,3.0)+(hf-pow(po,3.0))*lamac-lae;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Shu_Vos = " << ladc << endl;
return ladc; }
double PoiskPShumanVoss(double po)
{ 
	double lgl=1e-9, pgl=1e5, hf=1e0, pb=pgl, pa=lgl, pc=0.0;
	double tocras=1e-8, fa=hf, fb=hf, fc=hf, ra=pgl, ep=tocras; 
	int h=0, hko=1000;
while ((ra>ep) && (h<hko)) {
pc=(pa+pb)/2.0;
fa=(pow(pa,2.0)+pa)*log(hf+hf/pa)-pa-po;
fb=(pow(pb,2.0)+pb)*log(hf+hf/pb)-pb-po;
fc=(pow(pc,2.0)+pc)*log(hf+hf/pc)-pc-po;
if ((fc*fb>0.0) && (fa*fc<0.0)) pb=pc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) pa=pc;
ra=fabs(pa-pb); h++; } 
return pc; }
double ReshNeyaUravKoTeplGoringChirchillSpher(double lae, double lan, double po)
{ 
	double tocras=1e-8, lgl=1e-9, pgl=1e5, hf=1e0, ladb=pgl, lada=lgl, ra=pgl, ep=tocras;
	double vd=hf-po, ladc=hf, fa=hf, fb=hf, fc=hf; 
	int h=0, hko=1000; //7
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0; 
fa=lae/lan-(2e0+lada/lan-2e0*vd*(hf-lada/lan))/(2e0+lada/lan+vd*(hf-lada/lan));
fb=lae/lan-(2e0+ladb/lan-2e0*vd*(hf-ladb/lan))/(2e0+ladb/lan+vd*(hf-ladb/lan));
fc=lae/lan-(2e0+ladc/lan-2e0*vd*(hf-ladc/lan))/(2e0+ladc/lan+vd*(hf-ladc/lan));
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Gor_Chi_Sph = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplGoringChirchillCyl(double lae, double lan, double po)
{ 
	double lgl=1e-9, pgl=1e5, hf=1e0, ladb=pgl, ladc=0.0, lada=lgl, ra=pgl;
	double tocras=1e-8, ep=tocras, vd=hf-po, fa=hf, fb=hf, fc=hf; 
	int h=0, hko=1000; //8
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0; 
fa=((hf+lada/lan)/(hf-lada/lan)-vd)/((hf+lada/lan)/(hf-lada/lan)+vd)-lae/lan;
fb=((hf+ladb/lan)/(hf-ladb/lan)-vd)/((hf+ladb/lan)/(hf-ladb/lan)+vd)-lae/lan;
fc=((hf+ladc/lan)/(hf-ladc/lan)-vd)/((hf+ladc/lan)/(hf-ladc/lan)+vd)-lae/lan;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Gor_Chi_Cyl = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplBruggman(double lae, double lan, double po)//Модель Бруггмана-Ханаи - только для круглых частиц
{ 
	double lgl=1e-9, pgl=1e5, hf=1e0, ladb=pgl, lada=lgl, ladc=0.0, ra=pgl;
	double tocras=1e-8, ep=tocras, fa=hf, fb=hf, fc=hf, ot=hf/3e0; 
	int h=0, hko=1000; //9
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2.0;
fa=po-(lada-lae)*pow((lan/lae),ot)/(lada-lan);
fb=po-(ladb-lae)*pow((lan/lae),ot)/(ladb-lan);
fc=po-(ladc-lae)*pow((lan/lae),ot)/(ladc-lan);
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Brug = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplMaxwAcken(double lae, double lan, double po)
{ 
	double tocras=1e-8, lgl=1e-9, pgl=1e5, ladb=pgl, lada=lgl, ra=pgl, ep=tocras;
	double hf=1e0, vd=hf-po, fa=hf, fb=hf, fc=hf, ladc=hf; 
	int h=0, hko=1000; //10
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2.0;
fa=(hf-lan/lada)/(2.0*lan/lada+hf);
fa=(hf+2.0*vd*fa)/(hf-vd*fa)-lae/lan;
fb=(hf-lan/ladb)/(2.0*lan/ladb+hf);
fb=(hf+2.0*vd*fb)/(hf-vd*fb)-lae/lan;
fc=(hf-lan/ladc)/(2.0*lan/ladc+hf);
fc=(hf+2.0*vd*fc)/(hf-vd*fc)-lae/lan;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "lam_Maxw_Aken = " << ladc << endl;
return ladc; }
double KoTeplFrickEllipsoidChasti(double lae, double lan, double po)
{ 
	int h=0, hko=1000, nko=hko, k=3;
	double tocras=1e-8, lgl=1e-9, pgl=1e5, hf=1e0, ladb=pgl, lada=lgl, ra=pgl, ep=tocras, s=0.0; //11
	double AsRa=99e-2, b=2.0, c=b, a=AsRa*b, sc=s, ov=hf/2e0;
	double sa=s, sb=s, M=s, ladc=s, fa=s, fb=s, fc=s;
	double p0=s, e=s, lb=s, *xa=new double[k], *labc=new double[k];
	double *pea=new double[k], *peb=new double[k], *pec=new double[k];
	if ((!xa) || (!labc) || (!pea) || (!peb) || (!pec)) { cout << "No memory" << endl; k=getchar(); exit(1); } 
if (a>b) { e=pow(a,2.0)-pow(b,2.0); p0=pow(e,ov)/a; 
M=-2.0*p0+log(fabs((hf+p0)/(hf-p0))); M=M/pow(e,(3.0*ov)); }
else { e=pow(b,2.0)-pow(a,2.0); p0=pow(e,ov)/a; 
M=-p0+atan(p0); M=M*2.0/pow(e,(3.0*ov)); }
lb=2.0/a/pow(b,2.0)-M; lb=lb*ov;
k=0; labc[k]=M; k++; labc[k]=lb; k++; labc[k]=lb;
for (k=0; k<3; k++)
    xa[k]=2.0/a/b/c/labc[k]-hf;
while ((ra>ep) && (h<nko)) {
ladc=(lada+ladb)/2.0;
sa=0.0; sb=0.0; sc=0.0;
for (k=0; k<3; k++) {
    pea[k]=(hf+xa[k])/(xa[k]+lan/lada); sa=sa+pea[k];
    peb[k]=(hf+xa[k])/(xa[k]+lan/ladb); sb=sb+peb[k];
    pec[k]=(hf+xa[k])/(xa[k]+lan/ladc); sc=sc+pec[k]; }
fa=lada+po*(lan-lada)*sa/3.0-lae; //объемная доля растворенного вещества (эллипсоидов)
fb=ladb+po*(lan-ladb)*sb/3.0-lae;
fc=ladc+po*(lan-ladc)*sc/3.0-lae;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } 
if (xa) delete[]xa; if (labc) delete[]labc; if (pea) delete[]pea; 
if (peb) delete[]peb; if (pec) delete[]pec; //cout << "lam_Fri_Ellips = " << ladc << endl;
return ladc; }
double opredRdeltaSrVzaPro(double por, double laef, double lavo, int vyb)
{ 
	double L=1e0, Pa=1e1, Pdelt=1e2, delb=PoiskDelta(por,L), a=delb/Pa; //12, 13 - грубая и точная оценки соот-но
	double hf=1e0, x=0, la1=x, delm=x, lm=x, fb=hf, la1u=x, epf=1e-2;
	double pi=acos(-1e0), fbn=hf, fbnn=hf, ep=1e-3, y=x;
la1=raschTeplVzaimoPronKomp(delb,L,laef,lavo); //поиск КТП  твердого каркаса - грубая оценка
delm=L/Pdelt; //дельта малая - толщина трещин
lm=2.0*L-delm; //эль малая
x=2.0*delb/pow(pi,hf/2.0)/lm; y=a/delb;
fbn=opredFiBolN(x,y);
if (la1>ep) { fbnn=lavo/la1; fbnn=opredFiBolNN(fbnn, por); } else fbnn=0.0;
if (fbn>hf) 
	if (fbnn<hf) fb=fbnn; 
	else fb=0.0; 
else 
	if (fbnn>hf) fb=fbn; 
	else 
		if (fabs(fbn-fbnn)<epf) fb=(fbn+fbnn)/2.0; 
		else 
			if (fbn>fbnn) fb=fbnn; 
			else fb=fbn; //lam2 - пора, lam1 - твердый каркас //cout << "fb = " << fb << "\tx = " << x << "\ty = " << y << endl;
la1u=utochLamTverKark(delm,a,fb,lm,la1,delb,lavo,laef,L); //уточнение КТП твердого каркаса - тонкая оценка
if (vyb<13) { //cout << "lam_grub_otsen = " << la1 << endl; 
return la1; } else { //cout << "lam_tonk_ots = " << la1u << endl; 
return la1u; } }
double PoiskDelta(double m2, double L) //поиск дельты большой
{ 
	double hf=1e0, tocras=1e-8, ca=0.0, cb=hf, cc=ca, ra=fabs(ca-cb), ep=tocras, fa=cc, fb=cc, fc=cc; 
	int h=0, hko=1000, nit=hko; //cc - относительный размер бруса
while ((ra>ep) && (h<nit)) {
cc=(ca+cb)/2.0;
fa=2.0*pow(ca,3.0)-3.0*pow(ca,2.0)+hf-m2; 
fb=2.0*pow(cb,3.0)-3.0*pow(cb,2.0)+hf-m2; 
fc=2.0*pow(cc,3.0)-3.0*pow(cc,2.0)+hf-m2; 
if ((fc*fb>0.0) && (fa*fc<0.0)) cb=cc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) ca=cc; 
ra=fabs(ca-cb); h++; } 
return (cc*L); }
double raschTeplVzaimoPronKomp(double del, double L, double laef, double la_v) //взаимпопроникающие компоненты
{ 
	double lgl=1e-9, pgl=1e5, hf=1e0, R1a=hf, R1b=hf, R1c=hf, R2a=hf, R2b=hf, R2c=hf;
	double fa=hf, fb=hf, fc=hf, R=hf/laef/L, R3=hf/la_v/del;
	double R4=L/la_v/pow((L-del),2.0), la1b=pgl, la1a=lgl, la1c=hf;
	double tocras=1e-8, ra=pgl, ep=tocras; 
	int h=0, hko=1000;
while ((ra>ep) && (h<hko)) {
la1c=(la1a+la1b)/2.0;    
R1a=L/la1a/pow(del,2.0);
R1b=L/la1b/pow(del,2.0);
R1c=L/la1c/pow(del,2.0);
R2a=hf/la1a/(L-del);
R2b=hf/la1b/(L-del);
R2c=hf/la1c/(L-del);
fa=hf/R-hf/R1a-2.0/(R2a+R3)-hf/R4;
fb=hf/R-hf/R1b-2.0/(R2b+R3)-hf/R4;
fc=hf/R-hf/R1c-2.0/(R2c+R3)-hf/R4;
if ((fc*fb>0.0) && (fa*fc<0.0)) la1b=la1c; 
if ((fc*fa>0.0) && (fb*fc<0.0)) la1a=la1c; 
ra=fabs(la1a-la1b); h++; }
return la1c; }
double utochLamTverKark(double delm, double a, double phib, double lm, double la1, double delb, double lavo, double laef, double L) //структура с переменным сечением компонент
{ 
	double lgl=1e-9, pgl=1e5, pi=acos(-1e0), hf=1e0, R=hf/laef/L, b=2.0*a/pow(pi,0.5), delsh=delm/L;
	double r=2.0*delb/pow(pi,0.5), la2=lavo, y=b/r, laz=lavo;
	double R6=delm/2.0/laz/(pow(delb,2.0)-pow(a,2.0)), R4=hf, R3=hf;
	double phim=hf, la1a=hf, la1b=hf, ra=hf, ep=hf, la1c=hf, fa=hf, fb=hf, fc=hf;
	double lamsha=hf, lamshb=hf, lamshc=hf, R1a=hf, R1b=hf, R1c=hf, Rp0a=hf;
	double Rp0b=hf, Rp0c=hf, R7a=hf, R7b=hf, R7c=hf, R2a=hf, R2b=hf, R2c=hf; 
	double R5a=hf, R5b=hf, R5c=hf, Rdsa=hf, Rdsb=hf, Rdsc=hf, Rdelmina=hf;
	double tocras=1e-8, Rdelmaxa=hf, Rdelminb=hf, Rdelmaxb=hf, Rdelminc=hf, Rdelmaxc=hf; //laz - КТП компоненты, заполняющей трещины (воздух)
R4=L/la2/(pow((L-delb),2.0)); 
R3=hf/delb/la2;
phim=phib*pow(y,2e0)/(phib-pow(y,2e0));
la1a=lgl; la1b=la1+pgl; ra=pgl; ep=tocras; 
int h=0, hko=1000, nit=hko;
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
R5a=hf/(L-delb)/la1a;
R5b=hf/(L-delb)/la1b;
R5c=hf/(L-delb)/la1c;
Rdsa=hf/(R1a+R2a)+hf/(R6+R7a); Rdsa=hf/Rdsa;
Rdsb=hf/(R1b+R2b)+hf/(R6+R7b); Rdsb=hf/Rdsb;
Rdsc=hf/(R1c+R2c)+hf/(R6+R7c); Rdsc=hf/Rdsc;
Rdelmina=RdeltaMin(R1a,R6,la1a,delb,L,delsh);
Rdelmaxa=RdeltaMax(la1a,a,delb,L,delsh,R1a,R6);
Rdelminb=RdeltaMin(R1b,R6,la1b,delb,L,delsh);
Rdelmaxb=RdeltaMax(la1b,a,delb,L,delsh,R1b,R6);
Rdelminc=RdeltaMin(R1c,R6,la1c,delb,L,delsh);
Rdelmaxc=RdeltaMax(la1c,a,delb,L,delsh,R1c,R6);
fa=hf/Rdsa+hf/R4+2e0/(R3+R5a)-hf/R;
if ((Rdsa>=Rdelmina) && (Rdsa<=Rdelmaxa))
fa=hf/Rdsa+hf/R4+2e0/(R3+R5a)-hf/R;
fb=hf/Rdsb+hf/R4+2e0/(R3+R5b)-hf/R;
if ((Rdsb>=Rdelminb) && (Rdsb<=Rdelmaxb))
fb=hf/Rdsb+hf/R4+2e0/(R3+R5b)-hf/R;
fc=hf/Rdsc+hf/R4+2e0/(R3+R5c)-hf/R; 
if ((Rdsc>=Rdelminc) && (Rdsc<=Rdelmaxc))
fc=hf/Rdsc+hf/R4+2e0/(R3+R5c)-hf/R;  
if ((fc*fb>0.0) && (fa*fc<0.0)) la1b=la1c; 
if ((fc*fa>0.0) && (fb*fc<0.0)) la1a=la1c; 
ra=fabs(la1a-la1b); h++; } 
return la1c; }
double RdeltaMin(double R1, double R6, double la1, double delb, double L, double delsh)
{ 
	double hf=1e0, rdma=(hf/R1+hf/R6); 
	rdma=hf/rdma+(L-delsh/2.0)/la1/pow(delb,2.0);
return rdma; }
double RdeltaMax(double la1, double a, double delb, double L, double dsh, double R1, double R6)
{ 
	double hf=1e0, R2=L*(hf-dsh/2.0)/la1/pow(a,2.0);
	double R7=L*(hf-dsh/2.0)/la1/(pow(delb,2.0)-pow(a,2.0));
	double rdma=hf/(R1+R2)+1.0/(R6+R7);
	if (fabs(rdma)>0.0) rdma=hf/rdma;
return rdma; }
double ReshNeyaUravKoTeplGoringChirchillCylPlotSl(double laef, double lavo, double po, double srch)
{ 
	double tocras=1e-8, lgl=1e-9, pgl=1e5, x0=srch/2e0, lan=lavo, ra=0.0, ep=tocras, hf=1e0, ot=hf/3e0, ov=hf/2e0;
	int h=0, k=0, hko=1000; //cout << "x0 = " << x0 << "\tlaef = " << laef << "\tlavo = " << lan << endl; //x0=(2.+0.7)*1e-3/2./2.;
	double pi=acos(-1e0), kBa=ra, kBb=ra, kBc=ra, CBa=ra, CBb=ra, CBc=ra;
	double fa=ra, fb=ra, fc=ra, kC=ra, ladb=ra, ladc=ra, lada=ra; //14
kC=1e0; ladb=pgl; lada=lgl; ra=1e5; h=0; 
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0;
kBa=lada/kC/(lan-lada);
kBa=sign(kBa)*pow(fabs(kBa),ot); 
kBa=fabs(kBa);
CBa=pow(kBa,2e0)-kBa*x0+pow(x0,2e0);
CBa=fabs(CBa);
CBa=log(pow(CBa,ov)/fabs(kBa+x0));
CBa=CBa+pow(3e0,ov)*atan((2e0*x0-kBa)/pow(3e0*kBa,ov));
CBa=CBa-pow(3e0,ov)*atan(1e0/pow(3e0,ov));
CBa=CBa*pi*lan*lada/6e0/(lan-lada)/kC/kBa;
CBa=CBa+lada*(hf-pi*pow(x0,2e0)/4e0);
kBb=ladb/kC/(lan-ladb);
kBb=sign(kBb)*pow(fabs(kBb),ot); 
kBb=fabs(kBb);
CBb=pow(kBb,2e0)-kBb*x0+pow(x0,2e0);
CBb=fabs(CBb);
CBb=log(pow(CBb,ov)/fabs(kBb+x0));
CBb=CBb+pow(3e0,ov)*atan((2e0*x0-kBb)/pow(3e0*kBb,ov));
CBb=CBb-pow(3e0,ov)*atan(1e0/pow(3e0,ov));
CBb=CBb*pi*lan*ladb/6e0/(lan-ladb)/kC/kBb;
CBb=CBb+ladb*(hf-pi*pow(x0,2e0)/4e0);
kBc=ladc/kC/(lan-ladc);
kBc=sign(kBc)*pow(fabs(kBc),ot);
kBc=fabs(kBc);
CBc=pow(kBc,2e0)-kBc*x0+pow(x0,2e0);
CBc=fabs(CBc);
CBc=log(pow(CBc,ov)/fabs(kBc+x0));
CBc=CBc+pow(3e0,ov)*atan((2e0*x0-kBc)/pow(3e0*kBc,ov));
CBc=CBc-pow(3e0,ov)*atan(1e0/pow(3e0,ov));
CBc=CBc*pi*lan*ladc/6e0/(lan-ladc)/kC/kBc;
CBc=CBc+ladc*(hf-pi*pow(x0,2e0)/4e0);
fa=CBa-laef;
fb=CBb-laef;
fc=CBc-laef;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "la Gor Cher Cyl Plot Sloi = " << ladc << endl;
return ladc; }
double MetodHashinaShtrikmanaMin(double laef, double lan, double por) //оценка вилки - для изотропных в среднем образцов
{ 
	double tocras=1e-8, lgl=1e-9, pgl=1e5, lada=lgl, ladb=pgl, ra=pgl, ep=tocras, hf=1e0;
	double ta=hf, tb=hf, tc=hf, fa=hf, fb=hf, fc=hf, ladc=hf; 
	int h=0, hko=1000; //15
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2.0;
ta=hf/(lada-lan)+(hf-por)/3.0/lada; ta=lan+por/ta;
tb=hf/(ladb-lan)+(hf-por)/3.0/ladb; tb=lan+por/tb;
tc=hf/(ladc-lan)+(hf-por)/3.0/ladc; tc=lan+por/tc;
fa=ta-laef;
fb=tb-laef;
fc=tc-laef;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "la Hash Shtrik min = " << ladc << endl;
return ladc; }
double MetodHashinaShtrikmanaMax(double laef, double lan, double por)
{ 
	double tocras=1e-8, lgl=1e-9, pgl=1e5, lada=lgl, ladb=pgl, ra=pgl, ep=tocras, hf=1e0;
	double ta=hf, tb=hf, tc=hf, fa=hf, fb=hf, fc=hf, ladc=hf; 
	int h=0, hko=1000; //16
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0;    
ta=hf/(lada-lan)+por/3e0/lada; ta=lada+(hf-por)/ta;
tb=hf/(ladb-lan)+por/3e0/ladb; tb=ladb+(hf-por)/tb;
tc=hf/(ladc-lan)+por/3e0/ladc; tc=ladc+(hf-por)/tc;
fa=ta-laef;
fb=tb-laef;
fc=tc-laef;
if ((fc*fb>0.0) && (fa*fc<0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "la Hash Shtrik max = " << ladc << endl;
return ladc; }
double ReshNeyaUravKoTeplSwift(double lae, double lan)
{ 
	int h=0, hko=1000; //17
	double tocras=1e-8, lgl=1e-9, pgl=1e5, hf=1e0, lada=hf, ladb=hf, ra=pgl, ep=hf, ladc=hf;
	double fa=hf, fb=hf, fc=hf, d=(-2e0), mo=(-1e0), pi=acos(mo); 
lada=lgl; ladb=pgl; ra=pgl; ep=tocras; fa=0.0; fb=0.0; fc=0.0;
while ((ra>ep) && (h<hko)) {
ladc=(lada+ladb)/2e0;
fa=(log(fabs(lan/lada)))*pow((lan/lada-hf),d);
fa=pow((lan/lada-hf),mo)-fa;
fa=0.577*pi*fa+0.093-lae/lan;
fb=(log(fabs(lan/ladb)))*pow((lan/ladb-hf),d);
fb=pow((lan/ladb-hf),mo)-fb;
fb=0.577*pi*fb+0.093-lae/lan;
fc=(log(fabs(lan/ladc)))*pow((lan/ladc-hf),d);
fc=pow((lan/ladc-hf),mo)-fc;
fc=0.577*pi*fc+0.093-lae/lan;
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(lada-ladb); h++; } //cout << "la Swift = " << ladc << endl;
return ladc; }
double *opredTvChaKoeTepSrInSpSha(int vyb, double po, double *tem, double *laefm, int lear, double *srla, double *lavoz, double srch) //определение КТП твердой части иными способами
{ 
	double hf=1e0, lam=hf, lavo=hf, laef=hf, uo=urovPod(po);
	double tocras=1e-8, up=hf, ep=tocras, m2=po; 
	int k=0, f=0; //if (vyb<2) { for (k=0; k<lear; k++) cout << "tem = " << tem[k] << " lae = " << laefm[k] << " lavo = " << lavoz[k] << endl; }
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
{ 
	double hf=1e0; 
	if (x>0.0) return hf; 
	else if (x<0.0) return (-hf); 
	else return 0.0; }
double *DulnevKoefTepShamot(int vyb, double por, double *tem, int n, double *srk, double *laefm, double *lavo, 
	double srp, int koel, double *srra, double *rpr) //программа определяет КТП твердого каркаса шамота методом Дульнева
{ 
	int k=0; double *uk=NULL; //cout << "\tvyb = " << vyb << "\tkoel = " << koel; 
if (vyb<koel) uk=BezUchetaRaspedPorSha(por, tem, laefm, lavo, vyb, n, srk, srp); //1 - Кубы, адиаб., 2 - Кубы, изот., 3 - По Одолевскому 1, 4 - По Одолевскому 2, 5 - Цилиндры, адиаб., 6 - Цилиндры, изот.
else uk=UchetRapsredPorPoRazmSha(por, tem, laefm, lavo, vyb, n, srk, k, srra, rpr); //7 - Кубы, адиаб., 8 - Кубы, изотерм., 9 - Цилиндры, адиаб., 10 - Цилиндры, изот., 11 - По Одолевскому 1, 12 - По Одолевскому 2
return uk;
} 
double *DulnevKoefTepVerm(int vyb, double por, double *tem, int n, double *srk, double *laefm, double *lavo, double srp, 
	int koel, double *srra, double *raspr, int vybves, int vpv, int dmi, int nac) //определяет КТП твердого каркаса методом Дульнева
{ 
	int k=0; double *ukaz=NULL; //cout << "\tvyb = " << vyb << "\tkoel = " << koel; 
if ((vyb-nac)<koel) ukaz=BezUchetaRaspedPorSha(por, tem, laefm, lavo, vyb-nac, n, srk, srp); //1 - Кубы, адиаб., 2 - Кубы, изот., 3, 4 - По Одолевскому 1, 2, 5 - Цилиндры, адиаб., 6 - Цилиндры, изот.
else ukaz=UchetRapsredPorPoRazmVer(por, tem, laefm, lavo, vyb-nac, n, srk, k, srra, raspr, vybves, vpv, dmi); //7 - Кубы, адиаб., 8 - Кубы, изотерм., 9 - Цилиндры, адиаб., 10 - Цилиндры, изот., 11 - По Одолевскому 1, 12 - По Одолевскому 2
return ukaz; } 
double *BezUchetaRaspedPorSha(double por, double *tem, double *lameff, double *lamvoz, int vy, int n, double *lamvy, double srp) //без учета распределения пор
{ 
	double hf=1e0, ot=hf/3e0, m2=por, l=srp, lb=l/pow(por, ot), lavo=hf;
	double laef=hf, uo=urovPod(por), up=hf; int k=0; //cout << " l = " << l << "\tpor = " << por << "\tlb = " << lb << endl;;
for (k=0; k<n; k++) { 
laef=lameff[k]; lavo=lamvoz[k]; lamvy[k]=0.0;
if (!vy) lamvy[k]=DulKoefTep1Adi(lavo, por, laef); //куб, адиаб. //22 //адиабатное разбиение, полость - прямоугольный параллелепипед
if (vy==1) lamvy[k]=DulKoefTep1Izoterm(lavo, por, laef); //куб, изотерм. //23 //изотермическое разбиение, полость - прямоугольный параллелепипед
if (vy==2) lamvy[k]=DulKoefTep1OdolevMatr(lavo, por, laef); //по Одолевскому (матрицы) //24 //по Одолевскому (матрицы)
if (vy==3) lamvy[k]=DulKoefTep1OdolevStatSm(lavo, por, laef); //по Одолевскому (стат. смесь) //25 //по Одолевскому (стат. смесь)
if (vy==4) lamvy[k]=DulKoefTep1AdiCyl(lavo, por, laef, lb); //цилиндр., адиаб. //26 //адиабатное разбиение, полость - цилиндр
if (vy==5) lamvy[k]=DulKoefTep1IzotermCyl(lavo, por, laef, lb); //цилиндр., изотерм. //27 //изотермическое разбиение, полость - цилиндр 
if (vy==6) lamvy[k]=DulKoefTep1CylPer(lavo, por, laef, lb); //цилиндр, перпендикулярно тепловому потоку //28 //полость - цилиндр, тепловой поток направлен перпендикулярно оси
} 
lamvy=proverkakvi(lamvy, lameff, lamvoz, up, uo, n); //if (!vy) for (k=0; k<n; k++) cout << "Dul KTP Kub Adiab 22 = " << lamvy[k] << endl; //if (vy==1) for (k=0; k<n; k++) cout << "Dul KTP Kub Izoterm 23 = " << lamvy[k] << endl; //if (vy==2) for (k=0; k<n; k++) cout << "Dul KTP Odolev 24 = " << lamvy[k] << endl; if (vy==3) for (k=0; k<n; k++) cout << "Dul KTP Odolev 25 = " << lamvy[k] << endl; if (vy==4) for (k=0; k<n; k++) cout << "Dul Adiad Cyl 26 = " << lamvy[k] << endl; if (vy==5) for (k=0; k<n; k++) cout << "Dul Izoterm Cyl 27 = " << lamvy[k] << endl; if (vy==6) for (k=0; k<n; k++) cout << "Lam Cyl Per 28 = " << lamvy[k] << endl;
return lamvy; }
double *UchetRapsredPorPoRazmSha(double po, double *T, double *lameff, double *lamvoz, int vy, int n, 
	double *lamvy, int no, double *srra, double *rpr) //учет распределения пор
{ 
	int dmsrps=23, k=0, m=dmsrps, j=0;
	double hf=1e0, *vo=new double[m], *pw=new double[m], mx=hf, lcub=hf;
	double tocras=1e-8, lb=hf, ep=tocras, w=1e2, ot=hf/3e0, e=1e-6, t=e; 
if ((!vo) || (!pw)) { cout << "No memory" << endl; k=getchar(); exit(1); } 
if (po<hf) w=hf;
for (k=0; k<m; k++) { pw[k]=0.0; vo[k]=0.0; }
j=0; mx=0.0; t=mx;
for (k=0; k<m; k++) { 
    pw[k]=rpr[k]*po/w; //объемная доля поры заданного размера в полном (во всем) объеме пор
    if (fabs(pw[k])>e) vo[k]=pow(srra[k], 3e0)/pw[k]; else vo[k]=0.0; //все поры - кубы
	mx=mx+vo[k]; t=t+hf; } 
lcub = pow(mx, ot)/t; //оценка максимального значения l - размера образца
lb = lcub / pow(po, ot); //оценка размера поры
double lavo=hf, lame=hf, lamadit=hf, lamizot=hf, lamodot1=hf, up=hf, m2=hf;
double lamodot2=hf, lamcyladit=hf, lamcylizot=hf, lamcylpert=hf, uo=urovPod(po);
for (k=0; k<n; k++) {
    lavo=lamvoz[k]; lame=lameff[k]; lamadit=lame; lamizot=lame; lamodot1=lame; 
	lamodot2=lame; lamcyladit=lame; lamcylizot=lame; lamcylpert=lame;
    for (j=0; j<m; j++) { if (pw[j]>0.0) {
if (vy == 7) lamadit = DulKoefTep1Adi(lavo, pw[j], lamadit); //куб, адиаб. //29
if (vy == 8) lamizot = DulKoefTep1Izoterm(lavo, pw[j], lamizot); //куб., изотерм. //30
if (vy == 9) lamcyladit = DulKoefTep1AdiCyl(lavo, pw[j], lamcyladit, lb); //цил., адиаб. //31
if (vy == 10) lamcylizot = DulKoefTep1IzotermCyl(lavo, pw[j], lamcylizot, lb); //цил., изотерм. //32
if (vy == 11) lamodot1 = DulKoefTep1OdolevMatr(lavo, pw[j], lamodot1); //по Одолевскому (матрицы) //33
if (vy == 12) lamodot2 = DulKoefTep1OdolevStatSm(lavo, pw[j], lamodot2); //по Одолевскому (стат. смесь) //34
if (vy == 13) lamcylpert = DulKoefTep1CylPer(lavo, pw[j], lamcylpert, lb); } } //цил., перп. тепл. потоку //35
if (vy == 7) lamvy[k]=lamadit; 
if (vy == 8) lamvy[k]=lamizot;
if (vy == 9) lamvy[k]=lamcyladit; 
if (vy == 10) lamvy[k]=lamcylizot;
if (vy == 11) lamvy[k]=lamodot1; 
if (vy == 12) lamvy[k]=lamodot2;
if (vy == 13) lamvy[k]=lamcylpert; }
lamvy=proverkakvi(lamvy, lameff, lamvoz, up, uo, n); //if (vy == 13) for (k=0; k<n; k++) cout << " lam voz = " << lamvoz[k] << endl; if (vy == 7) for (k=0; k<n; k++) cout << "Dul KTP Adi Uchet (" << k << ") = " << lamvy[k] << endl;
if (pw) delete[]pw; if (vo) delete[]vo;
return lamvy; }
double *UchetRapsredPorPoRazmVer(double po, double *T, double *lameff, double *lamvoz, int vy, int n, double *lamvy, int no, 
	double *srra, double *raspr, int vybves, int pkv, int dmi) //учет распределения пор
{
	int k=0, j=0, m=dmi, vn=8;
	double tocras=1e-8, t=0.0, ht=1e0, w=1e2, ep=tocras, ot=ht/3e0, e=1e-6; //cout << "m = " << m << "\tpo = " << po << "\tt = " << t << endl;
t=0.0; for (k=0; k<m; k++) t=t+raspr[k]; //cout << "t = " << t << endl;
for (k=0; k<m; k++) raspr[k]=raspr[k]/t; //t=0.0; for (k=0; k<m; k++) t=t+raspr[k]; 
double *vo=new double[m], *pw=new double[m], mx=0.0, lcub=0.0, lb=0.0;
if ((!vo) || (!pw)) { cout << "No memory" << endl; k=getchar(); exit(1); } 
for (k=0; k<m; k++) { pw[k]=0.0; vo[k]=0.0; }
if (po<ht) w=ht;
t=0.0; mx=t;
for (k=0; k<m; k++) { 
    pw[k]=raspr[k]*po/w; //объемная доля поры заданного размера в полном (во всем) объеме пор
    if (fabs(pw[k])>e) vo[k]=pow(srra[k], 3e0)/pw[k]; else vo[k]=0.0; //все поры - кубы
	mx=mx+vo[k]; t=t+ht; } //if (vy==vn) { for (k=0; k<m; k++) cout << "raspr = " << raspr[k] << "\tsrra = " << srra[k] << "\t"; cout << endl; for (k=0; k<m; k++) cout << "pw = " << pw[k] << "\tvo = " << vo[k] << "\t"; cout << endl; }
if (fabs(t)>e) lcub = pow(mx, ot)/t; //оценка максимального значения l - размера образца
lb = lcub / pow(po, ot); //оценка размера поры //cout << "lcub = " << lcub << "\tlb = " << lb << "\tpo = " << po << endl;
double uo=urovPod(po), up=1e0, lavo=0.0, lame=0.0;
for (k=0; k<n; k++) {
	lame=lameff[k]; lavo=lamvoz[k]; lamvy[k]=lame; //cout << "ktp vozd = " << lavo << "\tte = " << T[k] << endl;
    for (j=0; j<m; j++) {
if (pw[j]>0.0) {
if (vy == 7) lamvy[k]=DulKoefTep1Adi(lavo, pw[j], lamvy[k]); //куб, адиаб. //29
if (vy == 8) lamvy[k]=DulKoefTep1Izoterm(lavo, pw[j], lamvy[k]); //куб., изотерм. //30
if (vy == 9) lamvy[k]=DulKoefTep1OdolevMatr(lavo, pw[j], lamvy[k]); //по Одолевскому (матрицы) //31
if (vy == 10) lamvy[k]=DulKoefTep1OdolevStatSm(lavo, pw[j], lamvy[k]); //по Одолевскому (стат. смесь) //32
if (vy == 11) lamvy[k]=DulKoefTep1AdiCyl(lavo, pw[j], lamvy[k], lb); //цил., адиаб. //33
if (vy == 12) lamvy[k]=DulKoefTep1IzotermCyl(lavo, pw[j], lamvy[k], lb); //цил., изотерм. //34
if (vy == 13) lamvy[k]=DulKoefTep1CylPer(lavo, pw[j], lamvy[k], lb); } } } //цил., перп. тепл. потоку //35 //cout << "uo = " << uo << "\tup = " << up << "\tpor = " << po << "\tm = " << m << endl; //if (vy == 13) for (k=0; k<n; k++) cout << "Lam po Cyl Perpen 35 = " << lamvy[k] << endl; //cout << endl;
lamvy=proverkakvi(lamvy, lameff, lamvoz, up, uo, n); //if (vy == 8) { for (k=0; k<n; k++) cout << "Lam Kub Izot 30 = " << lamcubizo[k] << "\tlam = " << lamvy[k] << endl; cout << endl; } //if (vy == 13) { for (k=0; k<n; k++) cout << "Lam Cyl Per 35 = " << lamcylper[k] << "\tlam = " << lamvy[k] << endl; cout << endl; } //if (vy == 9) for (k=0; k<n; k++) cout << "Lam po Odolevskomu Matr 31 = " << lamvy[k] << endl; if (vy == 10) for (k=0; k<n; k++) cout << "Lam po Odolevskomu Stat Sm 32 = " << lamvy[k] << endl;//if (vy == 11) for (k=0; k<n; k++) cout << "Lam po Cyl Adiab 33 = " << lamvy[k] << endl; if (vy == 12) for (k=0; k<n; k++) cout << "Lam po Cyl Izoterm 34 = " << lamvy[k] << endl; //if (vy == 13) for (k=0; k<n; k++) cout << "Lam po Cyl Perpen 35 = " << lamvy[k] << endl;
if (pw) delete[]pw; if (vo) delete[]vo;
return lamvy; }
double DulKoefTep1Adi(double lam2, double m2, double lame) //lam1 - КТП твердого каркаса, lam2 - КТП воздуха
{ 
	double tocras=1e-8, lgl=1e-9, pgl=1e5, hf=1e0, lamb=pgl, lama=lgl, ep=tocras, ra=pgl, ot=hf/3e0, dt=2e0/3e0;
	double fa=hf, fb=hf, fc=hf, nua=hf, nub=hf, nuc=hf, lamc=hf; //22
int k=0, hko=1000, nit=hko;
while ((ra>ep) && (k<nit)) {
lamc=(lama+lamb)/2e0; 
nua=lam2/lama; nub=lam2/lamb; nuc=lam2/lamc;
fa=(nua-(nua-hf)*pow(m2,ot)*(hf-pow(m2,dt)))/(nua-pow(m2,ot)*(nua-hf))-lame/lama;
fb=(nub-(nub-hf)*pow(m2,ot)*(hf-pow(m2,dt)))/(nub-pow(m2,ot)*(nub-hf))-lame/lamb;
fc=(nuc-(nuc-hf)*pow(m2,ot)*(hf-pow(m2,dt)))/(nuc-pow(m2,ot)*(nuc-hf))-lame/lamc;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } //cout << "lamc = " << lamc << "\tfc = " << fc << "\tra = " << ra << "\tk = " << k << "\n"; 
return lamc; }
double DulKoefTep1Izoterm(double lam2, double m2, double lame)
{ 
	double tocras=1e-8, lgl=1e-9, pgl=1e5, hf=1e0, lamb=pgl, lama=lgl, ep=tocras, ra=pgl, ot=hf/3e0, dt=2e0/3e0;
	double fa=hf, fb=hf, fc=hf, lamc=hf, nua=hf, nub=hf, nuc=hf; //23 //30
int k=0, hko=1000, nit=hko;
while ((ra>ep) && (k<nit)) {
lamc=(lama+lamb)/2e0;
nua=lam2/lama; nub=lam2/lamb; nuc=lam2/lamc;
fa=(hf+(nua-hf)*(pow(m2,dt)))/(hf+pow(m2,dt)*(nua-hf))*(hf-pow(m2,ot))-lame/lama;
fb=(hf+(nub-hf)*(pow(m2,dt)))/(hf+pow(m2,dt)*(nub-hf))*(hf-pow(m2,ot))-lame/lamb;
fc=(hf+(nuc-hf)*(pow(m2,dt)))/(hf+pow(m2,dt)*(nuc-hf))*(hf-pow(m2,ot))-lame/lamc;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } //cout << "lamc = " << lamc << "\tfc = " << fc << "\tra = " << ra << "\tk = " << k << "\n"; 
return lamc; }
double DulKoefTep1OdolevMatr(double lam2, double m2, double lame) //матричная гетерогенная система - одна фаза образует связную матрицу при любой объемной концентрации этой фазы, система имеет включения в виде кубов, центры которых образуют простую кубическую решетку, ребра параллельны
{ 
	double tocras=1e-8, lgl=1e-9, pgl=1e5, hf=1e0, lamb=pgl, lama=lgl, ep=tocras, ra=pgl, lamc=hf; 
	int k=0, hko=1000, nit=hko; //lam2 - КТП воздуха
	double nua=hf, nub=hf, nuc=hf, fa=hf, fb=hf, fc=hf;
while ((ra>ep) && (k<nit)) {
lamc=(lama+lamb)/2e0;
nua=hf-lama/lam2;
nub=hf-lamb/lam2;
nuc=hf-lamc/lam2;
fa=hf-(hf-m2)/(hf/nua-m2/3.0)-lame/lam2;
fb=hf-(hf-m2)/(hf/nub-m2/3.0)-lame/lam2;
fc=hf-(hf-m2)/(hf/nuc-m2/3.0)-lame/lam2;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } 
return lamc; }
double DulKoefTep1OdolevStatSm(double lam2, double m2, double lame)
{ 
	double tocras=1e-8, lgl=1e-9, pgl=1e5, hf=1e0, lamb=pgl, lama=lgl, ep=tocras, ra=pgl, lamc=hf, ov=hf/2e0; 
	int k=0, hko=1000, nit=hko; //нет регулярной структуры, статистическая смесь, частицы распределены хаотически
	double v1=hf, v2=hf, nua=hf, nub=hf, nuc=hf, fa=hf, fb=hf, fc=hf; //lam2 - КТП воздуха, ищем КТП твердого скелета
while ((ra>ep) && (k<nit)) {
lamc=(lama + lamb)/2.0;
v1=m2; v2=hf-m2;
nua=((3.0*v1-hf)*lama+(3.0*v2-hf)*lam2)/4.0;
nub=((3.0*v1-hf)*lamb+(3.0*v2-hf)*lam2)/4.0;
nuc=((3.0*v1-hf)*lamc+(3.0*v2-hf)*lam2)/4.0;
fa=nua+pow(pow(nua,2.0)+lama*lam2/2.0,ov)-lame; //lam2 - КТП в порах
fb=nub+pow(pow(nub,2.0)+lamb*lam2/2.0,ov)-lame; //lame - ЭКТП
fc=nuc+pow(pow(nuc,2.0)+lamc*lam2/2.0,ov)-lame;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } 
return lamc; }
double DulKoefTep1AdiCyl(double lam2, double m2, double lame, double d)
{ 
	double tocras=1e-8, lgl=1e-9, pgl=1e5, hf=1e0, fra=0.9, lamb=pgl, lama=lgl, ep=tocras, h=d;
	double ra=pgl, lamc=hf, h1=hf, d1=hf, F1=hf, R1a=hf, R1b=hf, R1c=hf, R2=hf, F12=hf, R3a=hf;
	double R3b=hf, R3c=hf, pi=acos(-1e0), F=hf, R=hf, fa=hf, fb=hf, fc=hf;
	int k=0, hko=1000, nit=hko;
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
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } 
return lamc; }
double DulKoefTep1IzotermCyl(double lam2, double m2, double lame, double d)
{ 
	double pi=acos(-1e0), tocras=1e-8, lgl=1e-9, pgl=1e5, hf=1e0, fra=0.9, lamb=pgl;
	double lama=lgl, ep=tocras, h=d, ra=pgl, lamc=hf, h1=hf, d1=hf, F1=hf, R1a=hf, R1b=hf;
	double R1c=hf, R2=hf, F12=hf, R3a=hf, R3b=hf, R3c=hf, F=hf, R=hf, fa=hf, fb=hf, fc=hf, ov=hf/2e0;
int k=0, hko=1000, nit=hko;
while ((ra>ep) && (k<nit)) {
    lamc=(lama+lamb)/2e0;
    h1=fra*d;
    d1=d*pow(4.0*m2/pi/fra,ov);
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
{ 
	double pi=acos(-1e0), tocras=1e-8, lgl=1e-9, pgl=1e5, hf=1e0, ov=hf/2e0, ra=pgl, lamb=pgl, lam1=lavo;
	double lama=lgl, dt=hf, h=hf, Qo=laef*h*dt, rb=pow(por/pi,ov)*lb; //lam1 - КТП воздуха, lam2 - КТП ТТ
	double lamc=hf, bb=hf, bm=hf, int1=hf, ab=hf, ep=tocras, am=hf;
	double Q1=hf, Q2=hf, fa=hf, fb=hf, fc=hf; //тепловой поток перпендикулярен оси цилиндра
	int k=0, hko=1000, nit=hko; 
while ((ra>ep) && (k<nit)) {
    lamc=(lama+lamb)/2e0;
	bb=dt*lam1*lama*h*rb; 
	bm=2e0*(lama-lam1)*rb; am=lam1*lb;
    if (fabs(am)>=fabs(bm)) { 
		ab=pow(pow(am,2e0)-pow(bm,2e0),ov); 
		int1=2e0*(atan((am+bm)/ab)-atan(bm/ab))/ab;}
	else { 
		ab=pow(pow(bm,2e0)-pow(am,2e0),ov); 
		int1=(log(fabs(am+bm-ab)/fabs(am+bm+ab))-log(fabs(bm-ab)/fabs(bm+ab)))/ab; }
    Q1=bb*(pi/2e0)/bm-am*bb/bm*int1; Q1=2e0*Q1; 
	Q2=lama*h*(lb-2e0*rb)/lb; 
	fa=Qo-Q1-Q2;
    bb=dt*lam1*lamb*h*rb; 
	bm=2e0*(lamb-lam1)*rb;
	if (fabs(am)>=fabs(bm)) { 
		ab=pow(pow(am,2e0)-pow(bm,2e0),ov); 
		int1=2e0*(atan((am+bm)/ab)-atan(bm/ab))/ab; }
    else { ab=pow(pow(bm,2e0)-pow(am,2e0),ov); 
	int1=(log(fabs(am+bm-ab)/fabs(am+bm+ab))-log(fabs(bm-ab)/fabs(bm+ab)))/ab; }
    Q1=bb*(pi/2e0)/bm-am*bb/bm*int1; Q1=2e0*Q1; 
	Q2=lamb*h*(lb-2e0*rb)/lb; fb=Qo-Q1-Q2;
    bb=dt*lam1*lamc*h*rb; 
	bm=2e0*(lamc-lam1)*rb;
	if (fabs(am)>=fabs(bm)) { 
		ab=pow(pow(am,2e0)-pow(bm,2e0),ov); 
		int1=2e0*(atan((am+bm)/ab)-atan(bm/ab))/ab; }
	else { ab=pow(pow(bm,2e0)-pow(am,2e0),ov); 
	int1=(log(fabs(am+bm-ab)/fabs(am+bm+ab))-log(fabs(bm-ab)/fabs(bm+ab)))/ab; }
    Q1=bb*(pi/2e0)/bm-am*bb/bm*int1; Q1=2e0*Q1; 
	Q2=lamc*h*(lb-2e0*rb)/lb; fc=Qo-Q1-Q2;
    if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; 
	if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
	k++; ra=fabs(lama-lamb); } 
return lamc; }
double *DulnevZern(double por, double *te, double *srk, double *laefm, double *lavo, double wsio, double walo, 
	double wmgo, double srch, int n, int dm, double *kuscv, double *tkuscv, double *stchv, double *Prvotk, double *tvtk, 
	int dmpv, int vybves)
{
	int k=0, f=0; 
	double hf=1e0, uo=0.0, up=hf, epsil=uo, ts=uo, la_voz=uo;
	double la_e=uo, lam=uo, m2=uo, koscve=hf, stchve=uo, Prvoz=hf;
uo=urovPod(por); //cout << "por = " << por << "\tuo = " << uo << "\tsrch = " << srch << endl; 
for (k=0; k<n; k++) {
ts=te[k];
stchve=opredKTPTKToch(stchv, te, ts, n);
epsil=stchve*koscve;
la_voz=opredKTPTKToch(lavo, te, ts, n);
la_e=opredKTPTKToch(laefm, te, ts, n);
Prvoz=opredKTPTKToch(Prvotk, tvtk, ts, dmpv);
lam=opredDulnLam1(por, ts, epsil, la_voz, la_e, srch, vybves, Prvoz); //cout << "eps = " << epsil << "\tPr_v = " << Prvoz << "\tla_voz = " << la_voz << "\tla_e = " << la_e << "\t"; 
srk[k]=lam; }
srk=proverkakvi(srk, laefm, lavo, up, uo, n); //for (k=0; k<n; k++) cout << "stchve ( " << k << " ) = " << stchv[k] << "\tlaefm = " << laefm[k] << "\tlavo = " << lavo[k] << "\tlam = " << srk[k] << endl;; //k=getchar(); 
return srk; }
double opredDulnLam1(double po, double T, double eps, double lavo, double lae, double srch, int vv, double Prvoz)
{ 
	double tocras=1e-7, hf=1e0, ladb=1e5, lada=1e-8, ra=hf, ep=tocras, ladc=0.0;
	double fa=0.0, fb=fa, fc=fa;
	int h=0, hko=1000, kit=hko;
	if ((vv<0) || (vv>=4)) { cout << "Oshibka v vybore materiala!" << endl; h=getchar(); exit(1); }
while ((ra>ep) && (h<kit)) {
ladc=(lada+ladb)/2e0;
if (vv>0) { 
	fa=DulnevSloForVer(po, T, eps, lavo, lada, srch, Prvoz, h);
	fb=DulnevSloForVer(po, T, eps, lavo, ladb, srch, Prvoz, h);
	fc=DulnevSloForVer(po, T, eps, lavo, ladc, srch, Prvoz, h); }
if (!vv) {
	fa=DulnevSloForSha(po, T, eps, lavo, lada, srch, Prvoz);
	fb=DulnevSloForSha(po, T, eps, lavo, ladb, srch, Prvoz);
	fc=DulnevSloForSha(po, T, eps, lavo, ladc, srch, Prvoz); }
if (vv<0) { cout << "Oshibka!" << endl; h=getchar(); exit(1); }
fa=fa-lae; fb=fb-lae; fc=fc-lae; //if (!h) { cout << "fa = " << fa << "\tfb = " << fb << "\tfc = " << fc << "\t"; }
if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc; 
ra=fabs(fa-fb); h++; } //cout << "h = " << h << "\t";
return ladc; }
double DulnevSloForVer(double m2, double T, double eps, double lavo, double lam1, double srch, double Prvoz, int nom)
{ 
	double r=srch/2e0, d=srch, epf=1e-2, hf=1e0, ot=hf/3e0, ov=hf/2e0; 
	double Nk=pow(pow(m2,2e0)-1e1*m2+9e0,ov); Nk=(m2+3e0+Nk)/2e0/m2; //y2=3.3e-3*pow(1-m2,-2./9.); pud=0; hsl=30e-3; rona=0.2; ro1=rona/(1-m2); y2=y2*pow(pud+9.8*ro1*(1-m2)*hsl,1./3.); y1=1e-20; y2=y1;
double eta=1e-4, y1=1e-3, y2=y1/pow(eta,ov); //y1=(10+50)*1e-4/2;
double y3=2e0*pow(Nk-hf,ov)/Nk, y4=y3/pow(hf-m2,ot); //hsh=2e-3/2;
double hsh=1e-3, fb=0.0, fbn=0.0, fbnn=0.0; //if (!nom) cout << "T = " << T << "\tr = " << r << "\tNk = " << Nk << "\ty3 = " << y3 << "\ty4 = " << y4 << "\t";
if (fabs(y3)>0.0) fbn=opredFiBolN(y2,y2/y3); if (fabs(lam1)>0.0) fbnn=opredFiBolNN(lavo/lam1,m2);
if (fbn>hf) 
	if (fbnn<hf) fb=fbnn; 
	else fb=0.0; 
else if (fbnn>hf) fb=fbn; 
else if (fabs(fbn-fbnn)<epf) fb=(fbn+fbnn)/2.0; 
else if (fbn>fbnn) fb=fbnn; else fb=fbn; 
double A=pow(y2,2e0)-pow(y1,2e0), F=pow(hf-pow(y2,2e0),ov);
double D=pow(hf-pow(y3,2e0),ov), E=pow(y4,2e0)-pow(y3,2e0);
double delsrsz=d*(hsh+hf/Nk), lamszm=MoleSostTeplVozd(T, delsrsz, lavo, Prvoz);
double epspr=eps/(2e0-eps), sig=5.668e-8;
double lamszl=4e0*epspr*sig*delsrsz*pow(T,3e0), lamsz=lamszm+lamszl;
double nusz=lamsz/lam1, w=(lavo/lamszm-nusz*D)/(hf-nusz);
double y12=pow(y1,2e0)/(hsh/2e0+(hf-hsh/2e0)*fb), DF=fabs(w-D)/fabs(w-F);
double nug=nusz, numz=MoleSostTeplVozd(T, hsh*r, lavo, Prvoz); numz=numz/lam1;
double nu2sp=nug; DF=(D-F+w*log(DF))*2.0*nug/(hf-nug); 
double AF=A/(hf-hsh/2e0-F+hsh/2e0/numz), ADF=hf/(AF+DF);
ADF=hf/(D/pow(y3,2e0)+ADF); 
ADF=ADF+nu2sp*E+y12; 
ADF=ADF*lam1/pow(y4,2e0); //if (!nom) cout << "ADF = " << ADF << "\tfb = " << fb << "\tA = " << A << "\tF = " << F << "\n";
return ADF; }
double *DopFunctOpredKTPTverChas(double *srk, double *tem, int n, double por, double *laefm, int p, int vfv, double *lavo, double r, double *stche, double *Prv, double *lamvoz, double *temvoz, int dmkov, double ep)
{ 
if (!p) srk=MetodDulnevaSigalovoyPolyDispVer(tem, n, por, laefm, vfv, r, stche, Prv, lamvoz, temvoz, dmkov, ep);
if (p==1) srk=MetodDulnevaSigalovoyDopVer(tem, n, por, laefm, vfv, r, stche, Prv, lamvoz, temvoz, dmkov, ep);
if (p==2) srk=MetodDulnevaSigalovoyBezIspVer(tem, n, por, laefm, vfv, r, stche, Prv, lamvoz, temvoz, dmkov, ep); 
if (p==3) srk=VasFrayObsh(laefm, lavo, por, n);
if ((p<0) && (p>3)) { cout << "Irregular number!"; p=getchar(); exit(1); }
double up=1e0, uo=urovPod(por); srk=proverkakvi(srk, laefm, lavo, up, uo, n);
return srk; }
double *MetodDulnevaSigalovoyPolyDispVer(double *tem, int lete, double m22, double *laefm, int vfv, double r, double *stch, double *Prvo, double *ktpvotk, double *tevotk, int dmkv, double tochras)
{ 
	double pi=acos(-1e0), ko=1e-2, porex=95.0*ko, hf=1e0, ot=hf/3e0, ov=hf/2e0;
	double mg=26.0*ko, m20=porex, m2=m20-m22; //m20 - общая пористость, m2 - межзерновая пористость (внешняя), m22 - внутренняя пористость
	if (m2<mg) m2=36.5*ko; //межзерновая пористость по Нижегородову
double A=pow(74.0*ko/(hf-m2),ot); 
double gu=9.81, hv=3e1*1e-3, rov=125.0, sigb=2.156e5 /*предел прочности при сжатии материала частиц, Па*/, ka=7e0/5e0, T, ladc, ladb, lada; 
double H0=101325.0, H=1e5, kB=1.380658e-23, angs=1e-10, wo2=21.0*ko, wn2=79.0*ko, do2=6e-1, dn2=65e-2; 
double d=(do2*wo2+dn2*wn2)*angs, sigse=pi*pow(d,2e0)/4e0, E=8.428e6; //E - модуль Юнга вермикулита
double akk=9e-1, ko1=0.2595/(1e0-0.2595), ko2=0.2595-ko1, laef=0.0; //ak - коэффициент аккомодации
double fp=0.2595/m2-(ko1*m2+ko2), ep=0.0, uo=0.0, up=1e0, ra=0.0;
int k=0, h=0, hko=1000, kit=hko; 
double *ladi=new double[lete], n=0.0, dliSvoPro=0.0, cPr=0.0, N=0.0, tc=3e0/4e0;
double dem1=0.0, dem2=0.0, dem4=0.0, B=0.0, lamg1=0.0, lamg2=0.0, lamg4=0.0, lamg5=0.0, izco=0.0;
double sigm1=0.0, sigm2=0.0, sigm3=0.0, sigm3a=0.0, sigm3b=0.0, sigm3c=0.0, sigm4=0.0, sigm5=0.0;
double dem5=0.0, lam0=0.0, sigmr=0.0, debo=0.0, rr=0.0, fa=0.0, fb=0.0, fc=0.0, ch1=0.0, ch2=0.0;
double ch3=0.0, ch123=0.0, ch4a=0.0, ch4b=0.0, ch4c=0.0, gi=0.0, ch5=0.0, ch6=0.0, sig=(5.67e-2);
if (!ladi) { cout << "No memory" << endl; k=getchar(); exit(1); } 
uo=urovPod(m22);
for (k=0; k<lete; k++) {
ladb=1e5; lada=1e-5; ra=1e2; ep=tochras; h=0;
while ((ra>ep) && (h<kit)) {
ladc=(lada+ladb)/2e0; T=tem[k];
laef=opredKTPTKToch(laefm, tem, T, lete);
izco=2e0*sig*pow(stch[k],2e0)*tc;
	lam0=opredKTPTKToch(ktpvotk, tevotk, T, dmkv); //температура в К
    n=H/kB/T;
    dliSvoPro=hf/pow(2e0,ov)/n/sigse; //длина свободного пробега молекулы воздуха
	cPr=opredKTPTKToch(Prvo, tevotk, T, dmkv); //число Прандтля
    N=pow(pow(m2,2e0)-1e1*m2+9.0,ov); N=(m2+3e0+N)/2e0/m2;
    debo=rov*gu*hv; //удельная нагрузка
    dem1=r*(A-pi/4e0); /*dem1=2.0*r*(A-pi/4.0);*/ if (dem1<0.0) dem1=1e0;
    dem2=7e0*r*pow(debo,7e0/18e0)*pow(sigb,2e0)/pow(N,2e0)/pow(E,8e0/9e0);
    dem4=2e0*r*(pow(2e0,ov)*A-hf); if (dem4<0.0) dem4=1e0; dem5=2e0*pow(2e0,ov)*r*A;
    akk=opredKoefAkkomodLandau(T); 
    B=4e0*ka/(ka+hf)/cPr*((2e0-akk)/akk)*H0*dliSvoPro;
    lamg1=lam0/(hf+B/H/dem1); lamg2=lam0/(hf+B/H/dem2); 
	lamg4=lam0/(hf+B/H/dem4); lamg5=lam0/(hf+B/H/dem5);
    sigm1=lamg1*pi*r*(A*log(A/(A-hf))-hf);
    sigm2=lamg2*pow(N,ov)*pow(E,4e0/9e0)*pow(debo,1e0/18e0)*r/3.2/pow(sigb,ov);
    sigm2=sigm2*fp; sigm3=(2e-2)*pow(N,ov)*pow(E,8e0/9e0)*pow(debo,11e0/18e0)*r/pow(sigb,3e0/2e0);
    sigm3a=lada*sigm3; sigm3b=ladb*sigm3; sigm3c=ladc*sigm3;
    sigm4=lamg4*pi*r*(pow(2e0,ov)*A-hf)/2e0; //большая пора - цилиндр
    if (m2>=74e-2) { sigm5=lamg5*2e0*r*(2e0*pow(A,2e0)-pi)/pow(2e0,ov)/A; rr=1e0+pow(A,2e0)/pi; }
    else { sigm5=0.0; rr=1e0+pow(pow(2e0,ov)*A-hf,2e0)/2e0; }
    gi=izco*pow(T/1e2,3e0); 
    sigmr=pi*pow(r,2e0)*gi*rr;
    fa=(pow(2e0,ov)*(sigm1+sigm2+sigm4/2+sigm5+sigmr)+sigm3a)/A/r-laef;
    fb=(pow(2e0,ov)*(sigm1+sigm2+sigm4/2+sigm5+sigmr)+sigm3b)/A/r-laef;
    fc=(pow(2e0,ov)*(sigm1+sigm2+sigm4/2+sigm5+sigmr)+sigm3c)/A/r-laef;
    if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
	if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc;
	ra=fabs(fa-fb); h++; }
ladi[k]=ladc; } 
return ladi; }
double *MetodDulnevaSigalovoyDopVer(double *tem, int lete, double m22, double *laefm, int vfv, double r, double *stch, double *Prvo, double *ktpvo, double *tevotk, int dmkv, double tochras) //m22 - пористость вермикулита фракции 2-0,7 мм, из статьи "Теплопроводность моно- и полидисперсных материалов"
{ 
	double pi=acos(-1e0), ko=1e-2, porex=95.0*ko, hf=1e0, ot=hf/3e0, ov=hf/2e0;
	double m20=porex, m2=m20-m22, mg=26.0*ko; 
	if (m2<mg) m2=36.5*ko; //общая пористость, межзерновая пористость по Нижегородову
double A=pow(74.0*ko/(hf-m2),ot) /*alpha - размер воздушной оболочки (ореола)*/, gu=9.8, hv=30e-3, rov=125., sigb=2.156e5; //предел прочности при сжатии материала частиц, Па
double ka=7e0/5e0, H0=101325.0, H=1e5, kB=1.380658e-23, angs=1e-10, wo2=21.0*ko, wn2=79.0*ko, do2=6e-1, dn2=65e-2;
double d=(do2*wo2+dn2*wn2)*angs, sigse=pi*pow(d,2.0)/4.0, sig=(5.67e-2), tc=3e0/4e0;
double E=8.428e6 /*модуль Юнга вермикулита*/, akk=0.9 /*коэффициент аккомодации*/, fa=0.0, fb=0.0, fc=0.0, up=hf, uo=hf, izco=hf;
double ko1=0.2595/(hf-0.2595), ko2=0.2595-ko1, fp=0.2595/m2-(ko1*m2+ko2), lada=hf, ladb=hf, ladc=hf, ra=hf, ep=hf; 
int k=0, h=0, hko=1000, kit=hko;
double *ladi=new double[lete], sigm5=0.0, n=0.0, dliSvoPro=0.0, cPr=0.0, N=0.0, dem1=0.0, dem2=0.0, dem3=0.0, dem4=0.0, dem5=0.0;
double B=0.0, lamg1=0.0, sigm1=0.0, sigm2=0.0, sigm3a=0.0, sigm3b=0.0, sigm3c=0.0, sigm4=0.0; 
double gi=0.0, rr=0.0, sigmr=0.0, ch1=0.0, ch2=0.0, ch3=0.0, ch123=0.0, ch4a=0.0, ch4b=0.0;
double ch4c=0.0, ch5=0.0, ch6=0.0, Zs=0.0, Zss=0.0, V=0.0, G=0.0, lamrs=0.0, lamrss=0.0, X=0.0, T=0.0, laef=0.0, lam0=0.0;
if (!ladi) { cout << "No memory" << endl; k=getchar(); exit(1); } 
uo=urovPod(m22);
for (k=0; k<lete; k++) {
ladb=1e5; lada=1e-5; ra=1e2; ep=tochras; h=0; 
while ((ra>ep) && (h<kit)) {
ladc=(lada+ladb)/2e0; T=tem[k];
	laef=opredKTPTKToch(laefm, tem, T, lete); 
	izco=2e0*sig*pow(stch[k],2e0)*tc;
	lam0=opredKTPTKToch(ktpvo, tevotk, T, dmkv); //температура в К
    n=H/kB/T; dliSvoPro=1e0/pow(2.0,ov)/n/sigse; //длина свободного пробега молекулы воздуха
    cPr=opredKTPTKToch(Prvo, tevotk, T, dmkv); //число Прандтля
	N=pow(pow(m2,2e0)-10.0*m2+9e0,ov); N=(m2+3e0+N)/2e0/m2;
    dem1=2e0*r*(A-pi/4e0); if (dem1<0.0) dem1=2e0*r*(A-2e0/3e0); if (dem1<0.0) dem1=hf;
	dem3=2e0*r*(pow(2e0,ov)*A-hf); if (dem3<0.0) dem3=hf; dem5=2.0*pow(2.0,ov)*r*A;
    akk=opredKoefAkkomodLandau(T); 
    B=4e0*ka/(ka+hf)*((2e0-akk)/akk)*H0*dliSvoPro/cPr;
    X=4.45*(A*log(A/fabs(A-hf))-hf)/(hf+B/H/dem1); //X - молекулярный перенос тепла через ореолы вокруг частицы
    Zs=2.23*(pow(2e0,ov)*A-hf)/(hf+B/H/dem3); //Z, Z' - перенос тепла через поры
    Zss=2.23/(hf+B/H/dem3)/(pow(2e0,ov)*A-hf);
    V=(pow(A,2e0)-pi/2e0)/(hf+B/H/dem5)/A; //V - через дополнительные пути (столбы) при пористости более 75 %
    G=izco*pow(T/1e2,3e0);
    lamrs=4.45*G*r*(hf+0.5*pow(pow(2e0,ov)*A-hf,2e0))/A;
    lamrss=4.45*G*r*(hf+pow(A,2e0)/pi);
    if (m2>=74e-2) { fa=lam0*(X+Zss+V)/A+lamrss+lada-laef; 
	fb=lam0*(X+Zss+V)/A+lamrss+ladb-laef; fc=lam0*(X+Zss+V)/A+lamrss+ladc-laef; }
    else if (m2>=26e-2) { fa=lam0*(X+Zs)/A+lamrs+lada-laef;    
    fb=lam0*(X+Zs)/A+lamrs+ladb-laef; fc=lam0*(X+Zs)/A+lamrs+ladc-laef; }
    if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
	if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc;
	ra=fabs(fa-fb); h++; }
ladi[k]=ladc; } 
return ladi; }
double *MetodDulnevaSigalovoyBezIspVer(double *tem, int lete, double m22, double *laefm, int vfv, double r, double *stch, double *Prvo, double *ktpvo, double *tevotk, int dmkv, double tochras) //m22 - пористость вермикулита
{
	double pi=acos(-1e0), ko=1e-2, porex=95.0*ko, hf=1e0, ot=hf/3e0, ov=hf/2e0;
	double m20 = porex /*общая пористость*/, m2=m20-m22, mg=26.0*ko; 
	if (m2<mg) m2=36.5*ko; //межзерновая пористость по Нижегородову //из статьи "Теплопроводность зернистых систем"
	double A=pow(74.0*ko/(hf-m2),ot) /*alpha - размер воздушной оболочки (ореола)*/, gu=9.8, hv=30e-3, rov=125.0, sigb=2.156e5; //предел прочности при сжатии материала частиц, Па
	double ka=7e0/5e0, H0=101325.0, H=1e5, kB=1.380658e-23, angs=1e-10, wo2=21.0*ko, wn2=79.0*ko, do2=6e-1, dn2=65e-2;
	double d=(do2*wo2+dn2*wn2)*angs, sigse=pi*pow(d,2e0)/4e0;
	double E=8.428e6 /*модуль Юнга вермикулита*/, Gms=20.58e6 /*модуль сдвига*/, akk=9e-1 /*коэффициент аккомодации*/;
	double fa=0.0, fb=0.0, fc=0.0, up=hf, uo=0.0;
	double ko1=0.2595/(hf-0.2595), ko2=0.2595-ko1, fp=0.2595/m2-(ko1*m2+ko2);
	double dem1=0.0, dem2=0.0, dem4=0.0, dem5=0.0, sig=(5.67e-2); 
int h=0, k=0, hko=1000, kit=hko;
double ra=0.0, ep=0.0, lada=0.0, ladc=0.0, ladb=0.0, lam0=0.0, n=0.0, N=0.0;
double dliSvoPro=0.0, B=0.0, ch1=0.0, ch2=0.0, ch3=0.0, ch123=0.0, ch4a=0.0;
double ch4b=0.0, ch4c=0.0, gi=0.0, ch5=0.0, ch6=0.0, izco=0.0;
double *ladi=new double[lete], cPr=0.0, lamg1=0.0, sigm1=0.0, sigm2=0.0;
double sigm3=0.0, sigm4=0.0, sigm5=0.0, T=0.0, debo=0.0; 
double sigm3a=0.0, sigm3b=0.0, sigm3c=0.0, rr=0.0, sigmr=0.0, laef=0.0;
if (!ladi) { cout << "No memory" << endl; k=getchar(); exit(1); } 
uo=urovPod(m22);
for (k=0; k<lete; k++) {
ladb=1e5; lada=1e-5; ra=1e2; ep=tochras; h=0; kit=1e3;
while ((ra>ep) && (h<kit)) {
ladc=(lada+ladb)/2e0; T=tem[k]; 
laef=opredKTPTKToch(laefm, tem, T, lete);
    lam0=opredKTPTKToch(ktpvo, tevotk, T, dmkv); //температура в К
	izco=2e0*sig*pow(stch[k],2e0)*3e0/4e0;
    n=H/kB/T;
    dliSvoPro=1e0/pow(2e0,ov)/n/sigse; //длина свободного пробега молекулы воздуха
    cPr=opredKTPTKToch(Prvo, tevotk, T, dmkv); //число Прандтля
    N=pow(pow(m2,2e0)-1e1*m2+9e0,ov); N=(m2+3e0+N)/2e0/m2;
    debo=rov*gu*hv; //удельная нагрузка
    dem1=r*(A-pi/4e0); if (dem1<0.0) dem1=hf;
    dem2=7e0*r*pow(debo,7e0/18e0)*pow(sigb,ov)/pow(N,ov)/pow(E,8e0/9e0);
    dem4=2e0*r*(pow(2e0,ov)*A-hf); if (dem4<0.0) dem4=hf;
    dem5=2e0*pow(2e0,ov)*r*A;
    akk=opredKoefAkkomodLandau(T); 
    B=4e0*ka/(ka+hf)/cPr*((2e0-akk)/akk)*H0*dliSvoPro;
    ch1=4.45*(A*log(A/fabs(A-hf))-hf)/(hf+B/H/dem1);
    ch2=fp*0.44*pow(N,ov)*pow(E,4e0/9e0)*pow(debo,1e0/18e0)/pow(sigb,ov)/(hf+B/H/dem2);
    ch3=2.23*(pow(2e0,ov)*A-hf)/(hf+B/H/dem4);
    ch4a=(2e-2)*fp*lada*pow(N,ov)*pow(E,8e0/9e0)*pow(debo,11e0/18e0)/pow(sigb,3e0/2e0)/A;
    ch4b=(2e-2)*fp*ladb*pow(N,ov)*pow(E,8e0/9e0)*pow(debo,11e0/18e0)/pow(sigb,3e0/2e0)/A;
    ch4c=(2e-2)*fp*ladc*pow(N,ov)*pow(E,8e0/9e0)*pow(debo,11e0/18e0)/pow(sigb,3e0/2e0)/A;
    gi=izco*pow(T/1e2,3e0);
    if ((m2>=2595e-4) && (m2<=74e-2)) {
        ch5=4.45*gi*r*(hf+pow(pow(2e0,ov)*A-hf,2e0)/2e0)/A;
		ch123=lam0*(ch1+ch2+ch3)/A;
		fa=ch123+ch4a+ch5-laef; fb=ch123+ch4b+ch5-laef; fc=ch123+ch4c+ch5-laef; }
    if (m2>74e-2) {
        ch3=2.23/(hf+B/H/dem4)/(pow(2e0,ov)*A-hf);
        ch6=(pow(A,2e0)-1.58)/A/(hf+B/H/dem5);
        ch5=4.45*gi*r*(hf+pow(A,2e0)/pi)/A;
        fa=(ch1+ch2+ch3+ch6)*lam0/A+ch4a+ch5-laef;
        fb=(ch1+ch2+ch3+ch6)*lam0/A+ch4b+ch5-laef;
		fc=(ch1+ch2+ch3+ch6)*lam0/A+ch4c+ch5-laef; }
    if ((fc*fb>0.0) && (fa*fc<0.0)) ladb=ladc; 
	if ((fc*fa>0.0) && (fb*fc<0.0)) lada=ladc;
	ra=fabs(fa-fb); h++; }
ladc=urovOtsechen(ladc, laef, uo); 
ladc=urovPodder(ladc, laef, up);
ladi[k]=ladc; } 
return ladi; }
double opredKoefAkkomodLandau(double T)
{ 
	double pi=acos(-1e0), PP=6.6260755e-34/2.0/pi, kB=1.380658e-23, NA=6.0221409e23, R=kB*NA, mu=29e-3;
	double a=3e-9, gamv=7.0/5.0, m=mu/NA, ro=opredPlotnVozd(T);
	double c=pow(gamv*R*T/mu,0.5), alpha=pow(a*c,2.0); alpha=(kB*T)/alpha; alpha=pow(alpha,3.0/2.0);
alpha=alpha/6.0/ro; alpha=alpha/pow(2.0*pi*m,0.5); 
return alpha; }
double opredPlotnVozd(double T)
{ 
	int dmvtk=28; double *te=arrTempAir(), *ro=arrPlotnAir(), r=opredKTPTKToch(ro, te, T, dmvtk); 
if (te) delete[]te; if (ro) delete[]ro;
return r; }
double *arrTempAir()
{ 
	int k=0, n=0, dmvtk=28;
	double *te=new double[dmvtk], tem0=273.15; 
if (!te) { cout << "No memory!"; k=getchar(); exit(1); } 
te[k]=0.0;  k++; te[k]=1e1;  k++; te[k]=2e1;  k++; te[k]=3e1;  k++; te[k]=4e1;  k++; te[k]=5e1;  k++;
te[k]=6e1;  k++; te[k]=7e1;  k++; te[k]=8e1;  k++; te[k]=9e1;  k++; te[k]=1e2;  k++; te[k]=12e1; k++;
te[k]=14e1; k++; te[k]=16e1; k++; te[k]=18e1; k++; te[k]=2e2;  k++; te[k]=25e1; k++; te[k]=3e2;  k++;
te[k]=35e1; k++; te[k]=4e2;  k++; te[k]=5e2;  k++; te[k]=6e2;  k++; te[k]=7e2;  k++; te[k]=8e2;  k++;
te[k]=9e2;  k++; te[k]=1e3;  k++; te[k]=11e2; k++; te[k]=12e2; k++; n=k;
for (k=0; k<n; k++) te[k]=te[k]+tem0;
return te; }
double *arrPlotnAir()
{ 
	int k=0, dmvtk=28;
	double *ar=new double[dmvtk]; 
if (!ar) { cout << "No memory!"; k=getchar(); exit(1); }
ar[k]=1.293; k++; ar[k]=1.247; k++; ar[k]=1.205; k++; ar[k]=1.165; k++; ar[k]=1.128; k++; ar[k]=1.093; k++;
ar[k]=1.06; k++; ar[k]=1.029; k++; ar[k]=1.0; k++; ar[k]=0.972; k++; ar[k]=0.946; k++; ar[k]=0.898; k++;
ar[k]=0.854; k++; ar[k]=0.815; k++; ar[k]=0.779; k++; ar[k]=0.746; k++; ar[k]=0.674; k++; ar[k]=0.615; k++;
ar[k]=0.566; k++; ar[k]=0.524; k++; ar[k]=0.456; k++; ar[k]=0.404; k++; ar[k]=0.362; k++; ar[k]=0.329; k++;
ar[k]=0.301; k++; ar[k]=0.277; k++; ar[k]=0.257; k++; ar[k]=0.239; k++; return ar; }
double opredUrovPodderM03(double por)
{ 
	int n=3, k=0; double p=0.0, s=p, *pam=new double[n], hf=1e0;
if (!pam) { cout << "No memory!" << endl; k=getchar(); exit(1); }
k=2; pam[k]=-8.93227577791347; k--; pam[k]=8.89444783404518; k--; pam[k]=1.06833435021354;
p=0.0; s=0.0; for (k=0; k<n; k++) { s=s+pam[k]*pow(por,p); p=p+hf; }
if (pam) delete []pam; 
s=hf/(hf-s*por);
return s; }
double *VasFraySha(double *tem, int n)
{ 
	int ke=9, k=0; double *te=new double[ke], *lsk=new double[ke], *vfs=new double[n];
if ((!te) || (!lsk) || (!vfs)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
k=0;	te[k]=1e2; k++; te[k]=18e1; k++; te[k]=22e1; k++; te[k]=3e2; k++; 
		te[k]=4e2; k++; te[k]=6e2;  k++; te[k]=8e2;  k++; te[k]=1e3; k++; te[k]=12e2; k++;
k=0;	lsk[k]=1.52; k++; lsk[k]=2.25; k++; lsk[k]=2.50; k++; lsk[k]=2.85; k++;
		lsk[k]=3.24; k++; lsk[k]=3.28; k++; lsk[k]=3.46; k++; lsk[k]=4.06; k++; lsk[k]=4.98; k++;
for (k=0; k<n; k++) { vfs[k]=0.0; vfs[k]=opredKTPTKToch(lsk, te, tem[k], ke); } 
if (lsk) delete []lsk; if (te) delete []te; 
return vfs; }
double SysPlastShakhPor(double laef, double lam2, double po)
{ 
	double tocras=1e-8, lgl=1e-9, pgl=1e5, hf=1e0, lamb=pgl, lama=lgl, ep=tocras;
	double ra=pgl, fa=hf, fb=hf, fc=hf, lamc=hf; //18
	int k=0, hko=1000, nit=hko;
while ((ra>ep) && (k<nit)) {
lamc=(lama+lamb)/2e0;
if (po<=5e-1) {
fa=lam2*(4.0*po/(hf+lam2/lama)+lama*(hf-2.0*po)/lam2)-laef;
fb=lam2*(4.0*po/(hf+lam2/lamb)+lamb*(hf-2.0*po)/lam2)-laef;
fc=lam2*(4.0*po/(hf+lam2/lamc)+lamc*(hf-2.0*po)/lam2)-laef; }
else if (po<hf) {
fa=lam2*(4.0*(hf-po)/(hf+lam2/lama)+2.0*(po-hf))-laef;
fb=lam2*(4.0*(hf-po)/(hf+lam2/lamb)+2.0*(po-hf))-laef;
fc=lam2*(4.0*(hf-po)/(hf+lam2/lamc)+2.0*(po-hf))-laef; }
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } //cout << "la Shakh poryadok = " << lamc << endl;
return lamc; }
double MetodStarostina(double laef, double lam2, double po)
{
	double tocras=1e-8, lgl=1e-9, pgl=1e5, lamb=pgl, lama=lgl, ep=tocras, hf=1e0, ra=pgl;
	double fa=hf, fb=hf, fc=hf, lamc=hf, dt=2e0/3e0; //19
	int k=0, hko=1000, nit=hko;
while ((ra>ep) && (k<nit)) {
lamc=(lama+lamb)/2e0;
fa=pow(lama,2.0)*pow(po,dt)+(lam2-lama)*lama;
fa=fa/(lama+(pow(po,dt)-po)*(lam2-lama))-laef;
fb=pow(lamb,2.0)*pow(po,dt)+(lam2-lamb)*lamb;
fb=fb/(lamb+(pow(po,dt)-po)*(lam2-lamb))-laef;
fc=pow(lamc,2.0)*pow(po,dt)+(lam2-lamc)*lamc;
fc=fc/(lamc+(pow(po,dt)-po)*(lam2-lamc))-laef;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } //cout << "la metod Starostina = " << lamc << endl;
return lamc; }
double MetodRusselNeprTverTelo(double laef, double lam2, double po)
{ 
	double tocras=1e-8, lgl=1e-9, pgl=1e5, hf=1e0, lamb=pgl, lama=lgl, ep=tocras;
	double ra=pgl, fa=hf, fb=hf, fc=hf, lamc=hf, dt=2e0/3e0; //20
	int k=0, hko=1000, nit=hko;
while ((ra>ep) && (k<nit)) {
lamc=(lama+lamb)/2e0;
fa=lama*po+(lama/lam2)*(hf-pow(po,dt));
fa=fa/(po-pow(po,dt)+(lam2/lama)*(hf-pow(po,dt)+po))-laef;
fb=lamb*po+(lamb/lam2)*(hf-pow(po,dt));
fb=fb/(po-pow(po,dt)+(lam2/lamb)*(hf-pow(po,dt)+po))-laef;
fc=lamc*po+(lamc/lam2)*(hf-pow(po,dt));
fc=fc/(po-pow(po,dt)+(lam2/lamc)*(hf-pow(po,dt)+po))-laef;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } //cout << "la metod Russelya nepr tver telo = " << lamc << endl;
return lamc; }
double MetodRusselNeprVozd(double laef, double lam2, double po)
{ 
	double tocras=1e-8, lgl=1e-9, pgl=1e5, hf=1e0, lamb=pgl, lama=lgl, ep=tocras, ra=pgl;
	double fa=hf, fb=hf, fc=hf, lamc=hf, dt=2e0/3e0; //21
	int k=0, hko=1000, nit=hko;
while ((ra>ep) && (k<nit)) {
lamc=(lama+lamb)/2e0;
fa=lama*pow(hf-po,dt)+hf-pow(hf-po,dt);
fa=fa/((lama/lam2)*(pow(hf-po,dt)-hf+po)+(2e0-pow(hf-po,dt)-po))-laef;
fb=lamb*pow(hf-po,dt)+hf-pow(hf-po,dt);
fb=fb/((lamb/lam2)*(pow(hf-po,dt)-hf+po)+(2e0-pow(hf-po,dt)-po))-laef;
fc=lamc*pow(1e0-po,dt)+hf-pow(hf-po,dt);
fc=fc/((lamc/lam2)*(pow(hf-po,dt)-hf+po)+(2e0-pow(hf-po,dt)-po))-laef;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
k++; ra=fabs(lama-lamb); } //cout << "la metod Russelya nepr vozd = " << lamc << endl;
return lamc; }
double *VasFrayObsh(double *laefm, double *lavo, double po, int n)
{
	int k=0, hko=1000, nit=hko, h=0; 
	double hf=1e0, mo=-hf, pi=acos(mo), tocras=1e-8, ot=hf/3e0, dt=2e0/3e0, ov=hf/2e0, L=hf;
	double x=Poiskkx((hf-po)/4e0), db=L*x;
	double hlm=x/(ov-x), hlb=hlm/(hf+hlm), lgl=1e-9, pgl=1e5, lamb=pgl, lama=lgl, ep=tocras;
	double ra=pgl, fa=hf, fb=hf, fc=hf, lamc=hf, mu=1e-2, E=8.428e6;
	double Pn=125.0*9.8*3e1*1e-3, Aa=hf, Ab=hf, Ac=hf, kk=(1.5+2e0)/2.0, kb=(2.2+2.9)/2e0;
	double km=(4e0+5e0)/2e0, nugza=hf, nugzb=hf, nugzc=hf, lamz=hf, nu=2e0*(hf-pow(mu,2e0))/E;
	double kc=(35e-2+45e-2)/2e0, rp=(725e-3)*pow(nu*Pn*L/2e0,1e0/3e0), lamka=hf, lamkb=hf, lamkc=hf;
	double lampna=hf, lampnb=hf, lampnc=hf, Q=pow(74e-2/(hf-po),ot), hsh=L*km*1e-3, sfk=pi*pow(rp,2e0);
	double lamksha=hf, lamkshb=hf, lamkshc=hf, nuga=hf, nugb=hf, nugc=hf, ko=1e-2, porex=95.0*ko, m20=porex, m2=po, uo=urovPod(po);
	double laef=hf, lavoz=hf, *lto=new double[n], up=hf; if (!lto) { cout << "No memory!" << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) {
	laef=laefm[k]; lavoz=lavo[k]; h=0;
while ((ra>ep) && (h<nit)) {
lamc=(lama+lamb)/2e0;
lamz=lavoz; nugza=lamz/lama; nugzb=lamz/lamb; nugzc=lamz/lamc;
if (Pn<3e5) { lampna=lama*pow(Pn,2e0/3e0)*kc/75.0/Q; 
lampnb=lamb*pow(Pn,dt)*kc/75.0/Q; 
lampnc=lamc*pow(Pn,dt)*kc/75.0/Q; }
else { lampna=lama*pow(Pn/E,4e0/9e0)*kb/Q; 
lampnb=lamb*pow(Pn/E,4e0/9e0)*kb/Q;
lampnc=lamc*pow(Pn/E,4e0/9e0)*kb/Q; }
lamksha=2e0*pow(2e0,ov)*sfk*lama/L/Q/hsh/kk; 
lamkshb=2e0*pow(2e0,ov)*sfk*lamb/L/Q/hsh/kk; 
lamkshc=2e0*pow(2e0,ov)*sfk*lamc/L/Q/hsh/kk;
lamka=lamksha+lampna;
lamkb=lamkshb+lampnb;
lamkc=lamkshc+lampnc;
Aa=pow(hlb,2e0)*1e3*nugza/4e0/kk/km+lamka/lama;
Ab=pow(hlb,2e0)*1e3*nugzb/4e0/kk/km+lamkb/lamb;
Ac=pow(hlb,2e0)*1e3*nugzc/4e0/kk/km+lamkc/lamc;
nuga=nugza; nugb=nugzb; nugc=nugzc;
fa=1e0/(pow(hf/hlb,2e0)+Aa)+nuga*pow(hf-hlb,2e0)+2e0/(hf+hlm+hf/(nuga*hlb))-laef/lama;
fb=1e0/(pow(hf/hlb,2e0)+Ab)+nugb*pow(hf-hlb,2e0)+2e0/(hf+hlm+hf/(nugb*hlb))-laef-lamb;
fc=1e0/(pow(hf/hlb,2e0)+Ac)+nugc*pow(hf-hlb,2e0)+2e0/(hf+hlm+hf/(nugc*hlb))-laef/lamc;
if ((fc*fb>0.0) && (fa*fc<0.0)) lamb=lamc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) lama=lamc; 
h++; ra=fabs(lama-lamb); } 
lto[k]=lamc; }
return lto; }
double Poiskkx(double kx)
{
	double tocras=1e-8, hf=1e0, xa=0.0, xb=5e-1, xc=hf, ra=fabs(xa-xb), ep=tocras, fa=hf, fb=hf, fc=hf; 
	int h=0, hko=1000, nit=hko; //xc - относительный размер бруса
while ((ra>ep) && (h<nit)) {
xc=(xa+xb)/2e0;
fa=4.0*pow(xa,3.0)-3.0*pow(xa,2.0)+kx; 
fb=4.0*pow(xb,3.0)-3.0*pow(xb,2.0)+kx; 
fc=4.0*pow(xc,3.0)-3.0*pow(xc,2.0)+kx; 
if ((fc*fb>0.0) && (fa*fc<0.0)) xb=xc; 
if ((fc*fa>0.0) && (fb*fc<0.0)) xa=xc; 
ra=fabs(xa-xb); h++; } 
return (xc); }
double urovPod(double po)
{ 
	double ko=1e-2, porex=95.0*ko, m2=porex, uo=0.0, porist=28.0*ko, pormakvi=53.0*ko;
	double porg=5e-1, porveit=6e-1, pormikvi=36.0*ko, hf=1e0;
if (po>porveit) m2=porex; else if (po>pormikvi) m2=pormakvi; else m2=po;
if (po<porist) uo=opredUrovPodderM03(po); else if (po<porg) uo=hf/(hf-po); else uo=hf/(hf-m2);
return uo; }
double *proverkakvi(double *ta, double *laefm, double *lavo, double up, double uo, int n)
{ 
	int k=0, f=0; double lvmi=1e3, e=1e-3; 
for (k=0; k<n; k++) if (lvmi>lavo[k]) lvmi=lavo[k];
for (k=0; k<n; k++) {
	ta[k]=urovOtsechen(ta[k], laefm[k], uo);
	ta[k]=urovPodder(ta[k], lavo[k], up); }
f=1; for (k=0; k<n; k++) if ((ta[k]<lvmi) && (f>0)) { f=-1; break; }
if (f<0) for (k=0; k<n; k++) ta[k]=0.0; 
return ta; }
double **NapMasVozd(char *snm)
{
	int k=0, dlma=28, c=3; double *ktpvoz=NULL, *tvoz=NULL, te0=273.15, *Prvo=NULL;
	ktpvoz=new double[dlma]; tvoz=new double[dlma]; Prvo=new double[dlma]; 
	if ((!ktpvoz) || (!tvoz) || (!Prvo)) { cout << snm << endl; k=getchar(); exit(1); }
	ktpvoz[k]=2.44; k++; ktpvoz[k]=2.51; k++; ktpvoz[k]=2.59; k++; ktpvoz[k]=2.67; k++; ktpvoz[k]=2.76; k++; ktpvoz[k]=2.83; k++;
	ktpvoz[k]=2.90; k++; ktpvoz[k]=2.96; k++; ktpvoz[k]=3.05; k++; ktpvoz[k]=3.13; k++; ktpvoz[k]=3.21; k++; ktpvoz[k]=3.34; k++;
	ktpvoz[k]=3.49; k++; ktpvoz[k]=3.64; k++; ktpvoz[k]=3.78; k++; ktpvoz[k]=3.93; k++; ktpvoz[k]=4.27; k++; ktpvoz[k]=4.60; k++;
	ktpvoz[k]=4.91; k++; ktpvoz[k]=5.21; k++; ktpvoz[k]=5.74; k++; ktpvoz[k]=6.22; k++; ktpvoz[k]=6.71; k++; ktpvoz[k]=7.18; k++;
	ktpvoz[k]=7.63; k++; ktpvoz[k]=8.07; k++; ktpvoz[k]=8.50; k++; ktpvoz[k]=9.15;
	for (k = 0; k < dlma; k++) ktpvoz[k] = ktpvoz[k] / 1e2; //КТП воздуха
	k = 0;
	tvoz[k]=0.0;  k++; tvoz[k]=1e1;  k++; tvoz[k]=2e1;  k++; tvoz[k]=3e1;  k++; tvoz[k]=4e1;  k++; tvoz[k]=5e1;  k++;
	tvoz[k]=6e1;  k++; tvoz[k]=7e1;  k++; tvoz[k]=8e1;  k++; tvoz[k]=9e1;  k++; tvoz[k]=1e2;  k++; tvoz[k]=12e1; k++;
	tvoz[k]=14e1; k++; tvoz[k]=16e1; k++; tvoz[k]=18e1; k++; tvoz[k]=2e2;  k++; tvoz[k]=25e1; k++; tvoz[k]=3e2;  k++;
	tvoz[k]=35e1; k++; tvoz[k]=4e2;  k++; tvoz[k]=5e2;  k++; tvoz[k]=6e2;  k++; tvoz[k]=7e2;  k++; tvoz[k]=8e2;  k++;
	tvoz[k]=9e2;  k++; tvoz[k]=1e3;  k++; tvoz[k]=11e2; k++; tvoz[k]=12e2;
	for (k=0; k<dlma; k++) tvoz[k]=tvoz[k]+te0; //температуры, при которых определяется КТП воздуха
	k=0; 
	Prvo[k]=0.707; k++; Prvo[k]=0.705; k++; Prvo[k]=0.703; k++; Prvo[k]=0.701; k++; Prvo[k]=.699; k++;
	Prvo[k]=0.698; k++; Prvo[k]=0.696; k++; Prvo[k]=0.694; k++; Prvo[k]=0.692; k++; Prvo[k]=0.69; k++;
	Prvo[k]=0.688; k++; Prvo[k]=0.686; k++; Prvo[k]=0.684; k++; Prvo[k]=0.682; k++; Prvo[k]=0.681; k++;
	Prvo[k]=0.68;  k++; Prvo[k]=0.677; k++; Prvo[k]=0.674; k++; Prvo[k]=0.676; k++; Prvo[k]=0.678; k++;
	Prvo[k]=0.687; k++; Prvo[k]=0.699; k++; Prvo[k]=0.706; k++; Prvo[k]=0.713; k++; Prvo[k]=0.717; k++;
	Prvo[k]=0.719; k++; Prvo[k]=0.722; k++; Prvo[k]=0.724; //Число Прандтля воздуха
	double **unau=new double*[c]; if (!unau) { cout << snm << endl; k=getchar(); exit(1); }
	k=0; unau[k]=ktpvoz; k++; unau[k]=tvoz; k++; unau[k]=Prvo;
	return unau;
}
double **opredSredKTPTK()
{
	int cemtk=11, k=1, f=cemtk, j=0, q=j, w=j, b=11, l=j, fl=1, na=j, *ceti=new int[f], u=j, v=7, n=f;
	int nnvyuv=2, kvf=0, jvsv=0, qvmi=0, qvuv=0, d=0, dmkov=3;
	int vyfrve=0, vysove=0, vymivmf=1, vyukve=1, c=0, nnvyfv=2, nnvysv=3, nnvmivmf=3, vybves=1, vybmar=0;
	double *ttk=NULL, cf=0.0, hf=1e0, e=1e-6, *cet=NULL, tna=2e2+273.15, dt=1e2, *ktptk=NULL;
	double **muv=new double*[b], *tem=new double[f], *cemt=new double[f], **mu=NULL, h=3e1*1e-3;
	double *nvyfv=new double[nnvyfv], *nvysv=new double[nnvysv], *nvmivmf=new double[nnvmivmf], **muc=new double*[b];
	double *nvyuv=new double[nnvyuv], **mut=new double*[b], *tepv=new double[cemtk], tf=0.0, *po=NULL; 
	char **mus=napolStrok(vybves, vybmar), *snm=mus[j], *sf=mus[k];
	if ((!nvyfv) || (!nvysv) || (!nvmivmf) || (!nvyuv)) 
	{ cout << snm << endl; k=getchar(); exit(1); } 
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
					mu=napMasEKTPVer(vyfrve, vysove, vymivmf, tem, k, h, vyukve, n, k, muv, f, dmkov, snm);
					k=0; ktptk=mu[k]; k++; ttk=mu[k]; k++; cet=mu[k]; k++; 
					for (u=k; u<v; u++) { po=mu[u]; if (po) delete[]po; } //cout << "\tce = " << cet[0] << endl;
					muv[b]=ktptk; mut[b]=ttk; muc[b]=cet; b++; 
					if (vysove>0) break; } }
			else if (vyfrve==1) {
				for (qvuv=0; qvuv<nnvyuv; qvuv++) {
					vyukve=nvyuv[qvuv];
					mu=napMasEKTPVer(vyfrve, vysove, vymivmf, tem, k, h, vyukve, n, k, muv, f, dmkov, snm);
					k=0; ktptk=mu[k]; k++; ttk=mu[k]; k++; cet=mu[k]; k++; 
					for (u=k; u<v; u++) { po=mu[u]; if (po) delete[]po; } //cout << "\tce = " << cet[0] << endl;
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
	if (mus) delete[]mus; if (snm) delete[]snm; if (sf) delete[]sf;
	f=3; muv=new double*[f]; k=0; muv[k]=ktpvy; k++; muv[k]=tevy; k++; muv[k]=cemvy;
	return muv;
}
double **pereraschetTverKark(double **mu)
{
	int k=0, n=0;
	double *raspr=NULL, *srra=NULL, *prgr=NULL, *legr=NULL, *mez=NULL, ko=1e-6, srp=0.0, mrp=0.0, hf=1e0, r=0.0;
	k=0; raspr=mu[k]; k++; srra=mu[k]; k++; prgr=mu[k]; k++; legr=mu[k]; k++; mez=mu[k];
	k=0; srp=mez[k]; k++; mrp=mez[k];
	k=0; r=ko; while (r<mrp) { r=r+hf; k++; } n=k;
	if (mrp<hf) ko=hf;
	for (k=0; k<n; k++) { srra[k]=srra[k]/ko; prgr[k]=prgr[k]/ko; legr[k]=legr[k]/ko; }
	k=0; mu[k]=raspr; k++; mu[k]=srra; k++; mu[k]=prgr; k++; mu[k]=legr; k++; mu[k]=mez;
	return mu;
}