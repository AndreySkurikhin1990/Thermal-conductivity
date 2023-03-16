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
#define vmivmf 1 //выбор метода измерений для фракции 2-0,7 мм: 0 - нестационарный, 1 - стационарный
#define vyfv 1 //выбор фракции: 0 - фракция 2-0,7 мм, 1 - фракция 8-4 мм, 2 - фракция 1,6-0,35 мм
#define vyuv 2 //выбор укладки: 1 - плоскопараллельная, 2 - вертикальная
#define vysv 2 //выбор состояния: 0 - исходное, 1 - после повторных измерений, 2 - после прокаливания при 1000 град С
using namespace std;
//-------
const double pi=acos(-1e0), porex=95.0*1e-2, tocras=1e-8;
const double ssiv84=467.0*1e-3, salv84=129.0*1e-3, smgv84=282.0*1e-3;
const double ssiv207=458.0*1e-3, salv207=14.0*1e-2, smgv207=29.0*1e-2;
const double porist=28.0*1e-2, tem0=273.15, ssiv=368.0*1e-3, salv=132.0*1e-3, smgv=217.0*1e-3, enu=1e-3;
const int hko=1000, vpkf=0, vpmf=0, isrp=0, dmkooscv=14, dmvtk=28, dkoscvl=6, dsov=60, dmkusc=14;
const int cemn=11, dso=12, dmko=4, cemdu=6, cemdum=2; //dmko - длина массива коэффициентов
const double ys0=3e1*1e-3, dete=1e2, te0=273.15, tn=22.0+te0, tnac=2e2+te0, templa=134.0*1e1+te0; 
//-------
int dlar=0, cem=cemn;
char *skpovvk=NULL, *sdvovvk=NULL, *sppovvk=NULL, *sppovv=NULL, *sdvovv=NULL;
char *skpovv=NULL, *ssdov=NULL, *szfov=NULL, *sfnoov=NULL;
double *Prvo=NULL, *vtetk=NULL, *ktpvotk=NULL, *stch=NULL, *etev=NULL;
//-------
double opredFiBolN(double, double);
double opredFiBolNN(double, double);
double opredDulnLam1(double, double, double, double, double, double);
double opredKoefAkkomodLandau(double);
double MoleSostTeplVozd(double, double, double);
double **napMasEKTPVerNac();
double novNapMas(int, int, int, int, int, int, int, int, int);
double *oslaintestepcher(int, int, double *);
double *RasKorOV(double, double, double, double, double, double);
double UsredMasOV(double **, int, int);
double *oslaintestepcher(int, int, double *);
double *DulnevZern(double, double *, double *, int, double *, double *, double, double, double, double, int, double, double, int, double *, double *, double *);
double DulnevSloForVer(double, double, double, double, double, double);
double opredPrVozd(double);
double *koefoslab(double, double, double, double *, int, double *);
double *oslaintestepcher(int, int, double *);
double *arrTempAir();
double *arrPlotnAir();
double opredUrovPodderM03(double);
double *proverkakvi(double *, double *, double *, double, double, int);
double urovPodder(double, double, double);
double **osvpam(double **, int);
double **polmat(double **, double *, int, int);
double *reshMetKram(double **, double *, int);
double vychopred(double  **, int);
double Determinant(double **, int);
double usredVelichPlank(double *, double *, double *, double, int, double);
double *Koef_Pogl_ver(double *);
double *PokazPrelomPribl(double *, double *, double *, int, double *, double);
double *izmMasChast(double *, double, int, double *);
double *PokazPrelomAl(double *, double *, double *, double, int, double *);
double *epsilnu(double *, int, double *);
double urovPod(double);
double urovOtsechen(double, double, double);
double trapz(double *, double *, int);
void osvpamov(int);
double *oprKTPTverKarkVermi(double *, double *laefm, double por, double wsio, double walo, double wmgo, int, int, double, double, int, double *, double *, double *, int, int);
double opredTeploprovVozd(double);
double opredPlotnVozd(double);
double opredKTPTKTochSha(double *, double *, double, int);
void napstrdir(int, int);
void NapMasVozdSha(double *, double *, int);
double *Kramers_Kronig_ver(double *, int, double *);
double **GetMatr(double **, double **, int, int, int);
double epsisredver(double, double *, double *, int, double *, double *, int);
double *dliny_voln_ver(double *, int);
double *opredKTPTverKarkVerm();
void initpervznach(int);
void NapMasPr();
double opredDulnLam1(double, double, double, double, double, double);
double *koefPribSha(double *, double *, int, double *, char *);
double **oprkoefKTPiskhchao(int, int, double *, double, int, double **, int, int, int);
double **opredtemphc(double *, double *, double *, double *, double *, int, int, double, double *, double *, double *, double **, char *);
double **vydelPol(int, int, double **, double **, int, int);
double **opredTempHolGor(double *, double *, int, int, double, int, double **, int, int, double *, int);
double **napMasEKTPVer(int, int, int, double *, int, double, int, int, int, double **, int, int);
double **arrTem_Netzsch(double **); //массив температур - экспериментальные данные на Netzsch - хаотичная засыпка фракции 2-0,7 мм (исходный)
double *arrKTP_Netzsch(); //массив КТП вермикулита
double *arrKTP_2020();
double *arrTem1_2020();
double **arrTem2_2020(double **);
double *arrTem3_2020();
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
double *danIskh207(double *, double*, int, int, int);
double **PoiskZavVelTem(int, double **, int, int, double *, double);
//-------
void main()
{
int k=0;
double *kttkv=opredKTPTverKarkVerm();
for (k=0; k<cem; k++) cout << "tem = " << etev[k] << "\tkttkv = " << kttkv[k] << endl;
k=getchar();
}
double *opredKTPTverKarkVerm()
{ int k=0, q=0, cemv=cem, vybves=1, vyma=0; 
double **mu=napMasEKTPVerNac(), *tem=NULL, *laefm=NULL, por=0.0, *po=NULL, *tkuscv=NULL, *kuscv=NULL;
double *dkoscvm=NULL, *dkoscvt=NULL, *stchsrver=NULL;
double nf=0.0, t=1e-1, hf=1e0, wsio=0.0, walo=0.0, wmgo=0.0, tena=0.0, dtos=0.0, g=0.0;
double tnoscv=3e2, dtoscv=1e2; 
k=5; po=mu[k]; nf=po[q]; k--; tem=mu[k]; k--; laefm=mu[k]; cemv=cem;
wsio=ssiv; walo=salv; wmgo=smgv; 
if ((!vyfv) || (vyfv==2)) { wsio=ssiv207; walo=salv207; wmgo=smgv207; } 
else if (vyfv==1) { wsio=ssiv84; walo=salv84; wmgo=smgv84; }
por=novNapMas(vybves, vyma, vmivmf, vyfv, vyuv, vysv, vpmf, vpkf, cemv); 
q=0; k=1; tena=tem[q]; dtos=tem[k]-tena;
dkoscvm = new double[dkoscvl]; dkoscvt = new double[dkoscvl]; stchsrver=new double[cemv];
	if ((!dkoscvm) || (!dkoscvt) || (!stchsrver)) { cout << "No memory!" << endl; k = getchar(); exit(1); }
	double tnd=6e2, dtd=2e2, tm=0.0, ko=1e-2; k=0; dkoscvt[k]=tnd; for (k=1; k<dkoscvl; k++) dkoscvt[k]=dkoscvt[k-1]+dtd;
	k = 0; dkoscvm[k] = 4.68; k++; dkoscvm[k] = 4.69; k++; dkoscvm[k] = 5.65; k++; 
	dkoscvm[k] = 13.17; k++; dkoscvm[k] = 20.2; k++; dkoscvm[k] = 27.81;
	for (k = 0; k < dkoscvl; k++) { tm=dkoscvm[k]*ko; dkoscvm[k]=hf-tm; } 
tkuscv = new double[dmkooscv]; kuscv = new double[dmkooscv];
	if ((!tkuscv) || (!kuscv)) { cout << "No memory!" << endl; k = getchar(); exit(1); }
	k=0; tkuscv[k] = tnoscv; for (k = 1; k < dmkooscv; k++) tkuscv[k] = tkuscv[k - 1] + dtoscv;
kuscv=koefoslab(wmgo, wsio, walo, tkuscv, dmkooscv, kuscv); 
for (k = 0; k<cemv; k++) { g = epsisredver(tem[k], tkuscv, kuscv, dmkooscv, dkoscvt, dkoscvm, dkoscvl); stchsrver[k] = g; } 
for (k=0; k<cemv; k++) cout << "tem = " << tem[k] << "\tst_ch = " << stchsrver[k] << "\t"; cout << endl;	
initpervznach(cemv);
double *ktptsv=oprKTPTverKarkVermi(tem, laefm, por, wsio, walo, wmgo, vyfv, cemv, tena, dtos, dmkusc, kuscv, tkuscv, stchsrver, vysv, vpkf);
if (mu) delete[]mu; 
return ktptsv; }
double *oprKTPTverKarkVermi(double *tem, double *laefm, double por, double wsio, double walo, double wmgo, int vfv, int n, double tena, double dtos, int dmkusc, double *kuscv, double *tkuscv, double *stchv, int vysove, int vpkf)
{ 
int k, j=0, q=0, vv=0, u=n, w=0, na=0; 
double *srk=NULL, *lavo=new double[n], rmi=0.0, rma=0.0, ko=1e-3, srch=0.0;
for (k=0; k<n; k++) lavo[k]=opredTeploprovVozd(tem[k]); //for (j=0; j<n; j++) cout << "j = " << j << "\ttem = " << tem[j] << "\tlavo = " << lavo[j] << endl; cout << endl; 
if (!vfv) { rma=2e0; rmi=7e-1; }
else if (vfv==1) { rma=8e0; rmi=4e0; }
else if (vfv==2) { rma=1.6; rmi=0.35; } srch=(rmi+rma)*ko/2e0; 
srk=new double[n]; if (!srk) { cout << "No memory!"; j=getchar(); exit(1); } for (j=0; j<n; j++) srk[j]=0.0;
srk=DulnevZern(por, tem, srk, na, laefm, lavo, wsio, walo, wmgo, srch, n, tena, dtos, dmkusc, kuscv, tkuscv, stchv);	
for (j=0; j<n; j++) cout << "k = " << w+1 << "\ttem = " << tem[j] << "\tlam_tk = " << srk[j] << endl; cout << endl; 
k=getchar();
return srk; }
double opredDulnLam1(double po, double T, double eps, double lavo, double lae, double srch)
{ double ladb=1e5, lada=1e-5, ra=1e2, ep=tocras, ladc, fa, fb, fc; int h=0, kit=hko;
while ((ra>ep) && (h<kit)) {
ladc=(lada+ladb)/2.0;
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
{
double gam=7e0/5e0, H=101325.0, H0=1e5, kB=1.380658e-23, n=H/kB/T;
double d=(6e-1*21e-2+65e-2*79e-2)*1e-10, sig=pi*pow(d,2.0)/4.0;
double dli=1e0/pow(2e0,0.5)/n/sig, Pr=opredPrVozd(T), a=opredKoefAkkomodLandau(T); 
double cz=8.42e-3, Ty=113.0, la0=cz/H0/(1e0+Ty/T), Kn=la0*H0/H/de, B=4e0*gam/(gam+1e0)*(2e0-a)/a*Kn/Pr;
return lamg/(1e0+B); }
double opredPrVozd(double T)
{ int n=dmvtk; double Pr=opredKTPTKTochSha(Prvo, vtetk, T, n);
return Pr; }
double opredKoefAkkomodLandau(double T)
{ double PP=6.6260755e-34/2.0/pi, kB=1.380658e-23, NA=6.0221409e23, R=kB*NA, mu=29e-3;
double a=3e-9, gamv=7.0/5.0, m=mu/NA, ro=opredPlotnVozd(T-tem0);
double c=pow(gamv*R*T/mu,0.5), alpha=pow(a*c,2.0); alpha=(kB*T)/alpha; alpha=pow(alpha,3.0/2.0);
alpha=alpha/6.0/ro; alpha=alpha/pow(2.0*pi*m,0.5); return alpha; }
double *koefoslab(double wmg, double wsi, double wal, double *tere, int n, double *kuo)
{ int lt=n, k; 
double wo=wmg+wsi+wal, kmg=0.0, kal=0.0, ksi=0.0, ht=1e0, eps=1e-6;
double *mgo=new double[lt], *alo=new double[lt], *sio=new double[lt];
if ((!mgo) || (!alo) || (!sio)) { cout << "No memory" << endl; k=getchar(); exit(1); } //cout << "wmg = " << wmg << "\twsio = " << wsi << "\twal = " << wal << "\two = " << wo << "\tn = " << lt << endl;
sio=oslaintestepcher(lt,0,sio);
alo=oslaintestepcher(lt,1,alo);
mgo=oslaintestepcher(lt,2,mgo); //for (k=0; k<n; k++) cout << "sio = " << sio[k] << "\talo = " << alo[k] << "\tmgo = " << mgo[k] << endl;
for (k=0; k<n; k++) {
kmg=mgo[k]; kal=alo[k]; ksi=sio[k];
if (kmg<0.0) kmg=0.0; if (kmg>ht) kmg=ht; 
if (kal<0.0) kal=0.0; if (kal>ht) kal=ht; 
if (ksi<0.0) ksi=0.0; if (ksi>ht) ksi=ht;
if (fabs(wo)>eps) 
kuo[k]=(kmg*wmg+kal*wal+ksi*wsi)/wo; 
else kuo[k]=0.0; } //for (k=0; k<n; k++) cout << "kuo = " << kuo[k] << "\ttere = " << tere[k] << "\t"; cout << endl;
delete []mgo; delete []alo; delete []sio; return kuo; }
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
double *RasKorOV(double a11, double a12, double a21, double a22, double b1, double b2)
{ int l=2, k; double **mas=new double*[l], *st=new double[l], *bs=new double[l], *x, *kor=new double[l];
if ((!mas) || (!st) || (!bs)) { cout << "" << endl; k=getchar(); exit(1);}
st[0]=a11; st[1]=a12; mas[0]=st; st=new double[l];
if (st) { st[0]=a21; st[1]=a22; mas[1]=st; } else { cout << "No memory" << endl; k=getchar(); exit(1);}
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
double opredTeploprovVozd(double temp)
{ double ktptsp=opredKTPTKTochSha(ktpvotk, vtetk, temp, dmvtk);
return ktptsp; }
double opredKTPTKTochSha(double *ktptks, double *te, double temp, int ce)
{
	int n = ce, f = 1, p = 0, k=0; 
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
void NapMasVozdSha(double *ktpvoz, double *tvoz, int dmvo)
{
	int k = 0, dlma = dmvo; double te0=273.15;
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
}
double epsisredver(double T, double *tdkusctm, double *dkusctm, int dkusctl, double *dkoscet, double *dkoscem, int dkoscel)
{ int k=1; napstrdir(k,k);
double *dv=NULL; dv=dliny_voln_ver(dv, 0);
dv=new double[dlar]; 
if (!dv) { cout << "No memory" << endl; k=getchar(); exit(1); } 
for (k=0; k<dlar; k++) dv[k]=0.0;
dv=dliny_voln_ver(dv, 1);
double *npp=new double[dlar];
if (!npp) { cout << "No memory" << endl; k=getchar(); exit(1); }
for (k=0; k<dlar; k++) npp[k]=0.0; //cout << "ce = " << dlar << "\t";
npp=Kramers_Kronig_ver(dv, dlar, npp); //cout << "k = " << k << endl;
double dkusct=opredKTPTKTochSha(dkusctm, tdkusctm, T, dkusctl); 
if (dkusct>1e0) dkusct=1e0; if (dkusct<0.0) dkusct=0.0; 
double dkosce=opredKTPTKTochSha(dkoscem, dkoscet, T, dkoscel); 
if (dkosce>1e0) dkosce=1e0; if (dkosce<0.0) dkosce=0.0; 
double ume=dkusct*dkosce; 
double *epsil=new double[dlar];
if (!epsil) { cout << "No memory" << endl; k=getchar(); exit(1); } 
epsil=epsilnu(npp, dlar, epsil); 
double epssr=usredVelichPlank(dv, epsil, npp, T, dlar, ume); //cout << "T = " << T << "\teps_s = " << epssr << "\tume = " << ume << "\n";
if (dv) delete []dv; if (npp) delete []npp; if (epsil) delete []epsil;
k=1; osvpamov(k); return epssr; }
double *Kramers_Kronig_ver(double *dv, int lear, double *nnu)
{	int k=0; ofstream fo;
fo.open(sppovv,ios_base::out | ios_base::trunc); 
if (!fo.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
	double c0=299792458.0, *nu=new double[lear], *nus=new double[lear], *nnus=new double[lear];
	double *alsr=new double[lear], e=1e-20; //for (k=0; k<lear; k++) alsr[k]=0.0;
	if ((!nu) || (!nnus) || (!nus) || (!alsr)) { cout << "No memory" << endl; k=getchar(); exit(1); } 
	alsr=Koef_Pogl_ver(alsr); 
	for (k=0; k<lear; k++) if (dv[k]>e) nu[k]=2e0*pi*c0/dv[k]; else nu[k]=0.0; //for (k=0; k<lear; k++) cout << "k = " << k << "\tnu = " << nu[k] << endl; cout << "lear = " << lear << endl; k=getchar();
nus=izmMasChast(nu, enu, lear, nus); //for (k=0; k<lear; k++) cout << "k = " << k << "\tnu = " << nus[k] << endl; cout << "lear = " << lear << endl; k=getchar();
nnus=PokazPrelomAl(nu, nus, alsr, c0, lear, nnus); //for (k=0; k<lear; k++) cout << "k = " << k << "\tnu = " << nnus[k] << endl; cout << "lear = " << lear << endl; k=getchar();
nnu=PokazPrelomPribl(nu, nus, nnus, lear, nnu, enu); //for (k=0; k<lear; k++) cout << "k = " << k << "\tnu = " << nnu[k] << endl; cout << "lear = " << lear << endl; k=getchar();
if (nu) delete[]nu; if (nnus) delete[]nnus; if (nus) delete[]nus; if (alsr) delete[]alsr;
for (k=0; k<lear; k++) { nnu[k]=nnu[k]-1e0+1.543; fo << nnu[k] << endl; } fo.close(); 
return nnu; }
double *dliny_voln_ver(double *dv, int ide)
{ ifstream fin; double p; int k, q=100, lear; char *s=new char[q]; 
if (!s) { cout << "No memory" << endl; k=getchar(); exit(1); } for (k=0; k<q; k++) s[k]='\0';
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
double *DulnevZern(double por, double *te, double *srk, int na, double *laefm, double *lavo, double wsio, double walo, double wmgo, double srch, int n, double tena, double dtosc, int dm, double *kuscv, double *tkuscv, double *stchv)
{ int k, f; double uo, up=1e0, epsil, ts, la_voz, la_e, lam, m2, koscve=1e0, stchve;
uo=urovPod(por);
for (k=0; k<n; k++) {
ts=te[k];
stchve=opredKTPTKTochSha(stchv, te, ts, n);
epsil=stchve*koscve; 
la_voz=opredKTPTKTochSha(lavo, te, ts, n);
la_e=opredKTPTKTochSha(laefm, te, ts, n);
lam=opredDulnLam1(por, ts, epsil, la_voz, la_e, srch); 
srk[k]=lam; }
for (k=0; k<n; k++) cout << "tem = " << te[k] << "\tlaefm = " << laefm[k] << "\tlavo = " << lavo[k] << "\tlam = " << srk[k] << endl; 
srk=proverkakvi(srk, laefm, lavo, up, uo, n); 
return srk; }
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
double urovPodder(double pro, double ref, double urpo)
{ double fl; int f=1; 
if (pro<0) f=0; else if (pro<(urpo*ref)) f=0; else f=1;
if (f) fl=pro; else fl=0.0; return fl; }
double urovOtsechen(double pro, double ref, double urot)
{ double fl; int f=1; 
if (pro<0.0) f=-1; else if (pro>(urot*ref)) f=-1; else f=1;
if (f>0) fl=pro; else fl=0.0; return fl; }
void napstrdir(int vyve, int vm)
{	int k; 
ssdov=new char[dsov]; sfnoov=new char[dsov];
for (k=0; k<dsov; k++) { ssdov[k]='\0'; sfnoov[k]='\0'; }
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
sfnoov[k]='F';  k++; sfnoov[k]='i';  k++; sfnoov[k]='l';  k++; sfnoov[k]='e';  k++; sfnoov[k]='_';  k++; sfnoov[k]='i'; k++;
sfnoov[k]='s';  k++; sfnoov[k]='_';  k++; sfnoov[k]='n';  k++; sfnoov[k]='o';  k++; sfnoov[k]='t';  k++; sfnoov[k]='_'; k++;
sfnoov[k]='o';  k++; sfnoov[k]='p';  k++; sfnoov[k]='e';  k++; sfnoov[k]='n';  k++; sfnoov[k]='!';  k++; sfnoov[k]='\0'; 
//1 - вермикулит
if (vyve==1) { skpovvk=new char[dsov]; sdvovvk=new char[dsov]; sppovvk=new char[dsov]; 
skpovv=new char[2*dsov]; sdvovv=new char[2*dsov]; sppovv=new char[2*dsov]; 
if ((!skpovvk) || (!sdvovvk) || (!sppovvk) || (!skpovv) || (!sdvovv) || (!sppovv)) { cout << "No memory" << endl; k=getchar(); exit(1); } 
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
}
void osvpamov(int vyve)
{ delete []ssdov; delete []sfnoov; 
if (vyve==1) { 
delete []sppovv; delete []sppovvk; delete []sdvovv; 
delete []sdvovvk; delete []skpovv; delete []skpovvk; }
}
double *reshMetKram(double **mat, double *ssc, int rm)
{ int n=rm, j;
double d, *m=new double[n], mm, **mmat; if (!m) { cout << "No memory" << endl; j=getchar(); exit(1); } //PrintMatr(mat, rm);
d=vychopred(mat,n); for (j=0; j<n; j++) { mmat=polmat(mat,ssc,j,n); mm=vychopred(mmat,n); 
m[j]=mm/d; //cout << "x ( " << j << " ) = " << m[j] << endl;
mmat=osvpam(mmat,n); } 
return m; }
double **osvpam(double **a, int n)
{ int j; double *b; for (j=0; j<n; j++) { b=a[j]; delete []b; } return NULL; }
double **polmat(double **mat, double *b, int no, int n)
{ int k, j; double **ma=new double*[n], *m; 
if (!ma) { cout << "No memory" << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) {
	m=new double[n]; if (!m) { cout << "No memory" << endl; k=getchar(); exit(1); }
	for (j=0; j<n; j++) {
		if (j==no) m[j]=b[k]; else m[j]=mat[k][j];}
ma[k]=m; } 
return (ma); }
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
if (!fn) { cout << "No memory" << endl; k=getchar(); exit(1); } 
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
	if (!s) { cout << "No memory" << endl; k=getchar(); exit(1); }
fin.open(skpovv,ios_base::in); for (k=0; k<q; k++) s[k]='\0'; //cout << "k = " << k << endl;
if (!fin.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
k=0; while ((!fin.eof()) && (k<dlar)) { fin.getline(s,q,'\n'); p=atof(s); kp[k]=p; k++; }
fin.close(); delete []s; return kp; }
double usredVelichPlank(double *dv, double *uv, double *npp, double tem, int n, double umensh)
{ double PP=6.6260755e-34, PB=1.380658e-23, c0=299792458.0, vl, c1, c2, lambda; 
int k; double *Ib, *Ibn, *dl, nc, nz, e=1e-20;
Ib=new double[n]; Ibn=new double[n]; dl=new double[n];
if ((!Ib) || (!Ibn) || (!dl)) { cout << "No memory" << endl; k=getchar(); exit(1); } 
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
double vychopred(double  **mas, int m) {
  double d=Determinant(mas, m); //Вычисление определителя
  return d; }
double Determinant(double **mas, int m) { // Рекурсивное вычисление определителя
  int i;
  double **p, k=1., d=0.0, *po;
  p=new double*[m]; if (!p) { cout << "No memory" << endl; i=getchar(); exit(1); }
  for (i=0; i<m; i++) { po=new double[m]; if (!po) { cout << "No memory" << endl; i=getchar(); exit(1); } p[i]=po; }
  d=0;
  if (m<1) cout << "Определитель вычислить невозможно!";
  if (m==1) { d = mas[0][0]; return (d); } 
  if (m==2) { d = mas[0][0]*mas[1][1]-mas[1][0]*mas[0][1]; return(d); } 
  if (m>2) {
	for (i=0; i<m; i++) {
		p=GetMatr(mas, p, i, 0, m); 
		d=d+k*mas[i][0]*Determinant(p, m-1); //разложение по строкам по первому столбцу
		k=-k; } } //(-1) в степени i
for (i=0; i<m; i++) { po=p[i]; delete []po; } 
return (d); }
double trapz(double *arx, double *ary, int lenarr)
{ int k; double s=0.0;
for (k=1; k<lenarr; k++) s=s+(ary[k]+ary[k-1])*(arx[k]-arx[k-1])/2e0;
return s; }
double **GetMatr(double **mas, double **p, int i, int j, int m) { //Получение матрицы без i-й строки и j-го столбца
  int ki, kj, di, dj;
  di = 0;
  for (ki = 0; ki<(m-1); ki++) { // проверка индекса строки
    if (ki == i) di = 1;
    dj = 0;
    for (kj = 0; kj<(m-1); kj++) { // проверка индекса столбца
      if (kj == j) dj = 1;
      p[ki][kj] = mas[ki + di][kj + dj];
}  } return p; }
void initpervznach(int lentem)
{ int k, u=lentem; ktpvotk=new double[dmvtk]; vtetk=new double[dmvtk]; 
if ((ktpvotk) && (vtetk)) NapMasVozdSha(ktpvotk, vtetk, dmvtk); else { cout << "No memory!" << endl; k=getchar(); exit(1); }
Prvo=new double[dmvtk]; if (!Prvo) { cout << "No memory!" << endl; k=getchar(); exit(1); } else NapMasPr();
}
void NapMasPr()
{ int k=0; Prvo[k]=0.707; k++; Prvo[k]=0.705; k++; Prvo[k]=0.703; k++; Prvo[k]=0.701; k++; Prvo[k]=.699; k++;
Prvo[k]=0.698; k++; Prvo[k]=0.696; k++; Prvo[k]=0.694; k++; Prvo[k]=0.692; k++; Prvo[k]=0.69; k++;
Prvo[k]=0.688; k++; Prvo[k]=0.686; k++; Prvo[k]=0.684; k++; Prvo[k]=0.682; k++; Prvo[k]=0.681; k++;
Prvo[k]=0.68; k++; Prvo[k]=0.677; k++; Prvo[k]=0.674; k++; Prvo[k]=0.676; k++; Prvo[k]=0.678; k++;
Prvo[k]=0.687; k++; Prvo[k]=0.699; k++; Prvo[k]=0.706; k++; Prvo[k]=0.713; k++; Prvo[k]=0.717; k++;
Prvo[k]=0.719; k++; Prvo[k]=0.722; k++; Prvo[k]=0.724; }
double **napMasEKTPVerNac()
{
	int k=1, nk=0, j=0, jk=0, f=cemdu, qn=0, q=0, w=f-1, cemv=cem; 
	double hf=1e0, *po=NULL, s=0.0, *thol, *tgor, *qob, *ektpv, *tsred, *ete, por;
	double nf=0.0, t=0.0, g=0.0, e=1e-6;
	if (cemv>0) ete=new double[cemv]; if (!ete) { cout << "No memory" << endl; k = getchar(); exit(1); }
	k=0; ete[k]=tnac; 
	for (k=1; k<cemv; k++) ete[k]=ete[k-1]+dete; //cout << "por = " << porver << endl;
	double **mu=new double*[f], **muv=new double*[f]; if ((!mu) || (!muv)) { cout << "No memory" << endl; k=getchar(); exit(1); } 
	muv=napMasEKTPVer(vyfv, vysv, vmivmf, ete, k, ys0, vyuv, cemv, k, muv, f, dmko); //for (k=0; k<cem; k++) cout << "te = " << ete[k] << "\t"; cout << endl; //k=w; po=muv[k]; k=0; nf=po[k]; t=e; k=0; while (t<nf) { t=t+hf; k++; } cem=k; cout << "cem = " << cem << "\tnf = " << nf << "\n"; for (k=0; k<w; k++) { po=muv[k]; for (j=0; j<cem; j++) cout << " ( " << j << " ) = " << po[j] << "\t"; cout << endl; }
	k=0; j=1; mu=vydelPol(k, cemv, muv, mu, f, j);
	k=0; thol=mu[k]; k++; tgor=mu[k]; k++; qob=mu[k]; k++; ektpv=mu[k]; k++; tsred=mu[k]; etev=tsred; k++; po=mu[k]; 
	k=0; nf=po[k]; t=0.0; k=0; while (t<nf) { t=t+hf; k++; } cemv=k; cem=k; //cout << "cem = " << cem << "\tnf = " << nf << "\t";	
	for (k=0; k<cemv; k++) cout << "tem = " << tsred[k] << "\tktp = " << ektpv[k] << "\tth = " << thol[k] << "\ttg = " << tgor[k] << "\tqo = " << qob[k] << "\t"; cout << endl;
	//for (k=0; k<f; k++) { po=muv[k]; delete[]po; po=mu[k]; delete[]po; } 
	//if (muv) delete[]muv; if (mu) delete[]mu; 
	delete[]ete; return mu; 
}
double **napMasEKTPVer(int vyfrve, int vysove, int vymeizvemafr, double *te, int v, double h, int vyukve, int nt, int wa, double **mu, int rmu, int dmk)
{
	int k=0, nvm=dmk, j=0, vy=0, n84=0, n207=0, cedu=rmu, q=0, cedumi=cemdum, nn=0;
	double F=13.85*1e-4, nf84=0.0, nf207=0.0, *po=NULL, *tc84=NULL, s=0.0, t=0.0, e=1e-1; 
	double *th84=NULL, hf=1e0, *tp84=NULL, *th207=NULL, *tc207=NULL, *tp207=NULL;
	double *koeft=NULL, *koefh=NULL, *koefc=NULL, *koefq=NULL, *koefs=NULL, **muv=NULL, r=0.0; //cout << "vf = " << vyfrve << "\tvs = " << vysove << "\tvmi = " << vymeizvemafr << "\tvu = " << vyukve << endl; //for (k=0; k<nt; k++) cout << "te ( " << k << " ) = " << te[k] << "\t";
	double *ktp=NULL, *temvs=NULL, *temvh=NULL, *temvc=NULL, *tepv=NULL, *vm=NULL, *ma=NULL; 
	//-------
	if ((!vyfrve) || (vyfrve==2)) { //фракция 2-0,7 мм или фракция 1,6-0,35 мм
		th207=arrTemHigh207(); 
		k=cedumi; muv=new double*[k]; if (!muv) { cout << "No memory" << endl; k = getchar(); exit(1); }
		muv=arrTemCold207(muv); k=0; po=muv[k]; k++; tc207=muv[k]; 
		tp207=arrTepPot207();
		k=0; t=0.0; nf207=po[k]; while (t<nf207) { t=t+hf; k++; } n207=k; nn=k; //cout << "n207 = " << n207 << "\t";
		temvs=new double[n207]; temvh=new double[n207]; temvc=new double[n207]; 
		tepv=new double[n207]; ktp=new double[n207]; 
		if ((!temvs) || (!temvh) || (!temvc) || (!tepv) || (!ktp)) { cout << "No memory" << endl; k=getchar(); exit(1); } //for (k=0; k<n207; k++) cout << "th = " << th207[k] << "\ttc = " << tc207[k] << "\tQ = " << tp207[k] << endl;
		for (k=0; k<n207; k++) {
		temvh[k]=th207[k]; temvc[k]=tc207[k]; temvs[k]=(temvh[k]+temvc[k])/2e0; 
		tepv[k]=tp207[k]/F; ktp[k]=fabs(temvh[k]-temvc[k])/h; ktp[k]=tepv[k]/ktp[k]; }
		if (tc207) delete[]tc207; if (th207) delete[]th207; if (tp207) delete[]tp207; if (po) delete[]po; //for (k=0; k<n207; k++) cout << "ktp = " << ktp[k] << "\tqob = " << tepv[k] << endl;
		if (!vysove) {
		if (!vymeizvemafr) { //установка Netzsch
		if (muv) delete[]muv; muv=new double*[cedu]; if (!muv) { cout << "No memory" << endl; k=getchar(); exit(1); }
		muv=oprkoefKTPiskhchao(vymeizvemafr, k, te, h, n207, muv, cedu, nt, dmk); 
		k=0; koefc=muv[k]; k++; koefh=muv[k]; k++; koefq=muv[k]; k++;
		koeft=muv[k]; k++; koefs=muv[k]; k++; po=muv[k]; }
		else if (vymeizvemafr==1) { //данные 2020 года - ГОСТ 12170 - стационарный метод
		if (ktp) delete[]ktp; ktp=arrKTP_2020(); if (temvh) delete[]temvh; temvh=arrTem1_2020(); 
		k=cedumi; muv=new double*[k]; if (!muv) { cout << "No memory" << endl; k=getchar(); exit(1); }
		muv=arrTem2_2020(muv);
		k=0; po=muv[k]; nf207=po[k]; k++; if (temvc) delete[]temvc; temvc=muv[k]; 
		if (temvs) delete[]temvs; temvs=arrTem3_2020(); 
		k=0; s=e; while (s<nf207) { s=s+hf; k++; } n207=k; k=0; q=3; //cout << "nk 2020 = " << n207 << endl; //for (k=0; k<n207; k++) cout << "ts ( " << k << " ) = " << temvs[k] << "\tktp = " << ktp[k] << "\tth = " << temvh[k] << "\ttc = " << temvc[k] << endl;
		double *ts=new double[n207], *tst=new double[n207]; if ((!ts) || (!tst)) { cout << "No memory" << endl; k=getchar(); exit(1); }
		for (k=0; k<n207; k++) { tst[k]=(temvh[k]+temvc[k])/2e0; ts[k]=tst[k]; temvs[k]=ts[k]; }
		ts=danIskh207(ts, tst, k, q, q); ktp=danIskh207(ktp, tst, k, q, q); temvh=danIskh207(temvh, tst, k, q, q); 
		temvc=danIskh207(temvc, tst, k, q, q); temvs=danIskh207(temvs, tst, k, q, q); if (tst) delete[]tst; //for (k=0; k<q; k++) cout << "ts ( " << k << " ) = " << ts[k] << "\tktp = " << ktp[k] << "\tth = " << temvh[k] << "\ttc = " << temvc[k] << endl;
		koeft=new double[dmk]; koefh=new double[dmk];
		koefc=new double[dmk]; koefs=new double[dmk]; koefq=new double[dmk];
		if ((!koeft) || (!koefh) || (!koefc) || (!koefs) || (!koefq)) { cout << "No memory" << endl; k=getchar(); exit(1); }
		for (k=0; k<dmk; k++) { koefh[k]=0.0; koefc[k]=0.0; koefq[k]=0.0; koeft[k]=0.0; koefs[k]=0.0; } 
		if (ktp) koeft=koefPribSha(ktp, ts, q, koeft, "No memory"); 
		if (temvh) koefh=koefPribSha(temvh, ts, q, koefh, "No memory");
		if (temvc) koefc=koefPribSha(temvc, ts, q, koefc, "No memory");
		if (temvs) koefs=koefPribSha(temvs, ts, q, koefs, "No memory");
		tepv=new double[q]; if (!tepv) { cout << "No memory" << endl; k=getchar(); exit(1); }
		for (k=0; k<q; k++) {
			s=0.0; r=0.0; for (j=0; j<dmk; j++) { s=s+koeft[j]*pow(ts[k], r); r=r+hf; }
			tepv[k]=s*fabs(temvh[k]-temvc[k])/h; }
			koefq=koefPribSha(tepv, ts, q, koefq, "No memory"); }
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
		k=0; t=0.0; nf84=po[k]; while (t<nf84) { t=t+hf; k++; } n84=k; nn=k; 
		temvs=new double[n84]; temvh=new double[n84]; temvc=new double[n84]; tepv=new double[n84]; ktp=new double[n84];
		if ((!temvs) || (!temvh) || (!temvc) || (!tepv) || (!ktp)) { cout << "No memory" << endl; k=getchar(); exit(1); }
		for (k=0; k<n84; k++) { 
		temvs[k]=(th84[k]+tc84[k])/2e0; 
		temvh[k]=th84[k]; temvc[k]=tc84[k]; tepv[k]=tp84[k]; 
		tepv[k]=tepv[k]/F; ktp[k]=fabs(th84[k]-tc84[k])/h; 
		ktp[k]=tepv[k]/ktp[k]; } //for (k=0; k<n84; k++) cout << "th = " << temvh[k] << "\ttc = " << temvc[k] << "\ttp = " << tepv[k] << endl;
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
	k=1; muv=new double*[cedu]; po=new double[k]; if ((!muv) || (!po)) { cout << "No memory" << endl; k=getchar(); exit(1); }
	k=0; po[k]=t; muv[k]=koefc; k++; muv[k]=koefh; k++; 
	muv[k]=koefq; k++; muv[k]=koeft; k++; muv[k]=koefs; k++; muv[k]=po; //vy=cedu-1; for (k=0; k<vy; k++) { ma=muv[k]; for (j=0; j<dmk; j++) cout << "k = " << k << "\tkoef ( " << j << " ) = " << ma[j] << "\t"; cout << endl; }
	q=cedu-1; for (vy=0; vy<q; vy++) { 
	vm=muv[vy]; ma=new double[nt]; if (!ma) { cout << "No memory" << endl; k=getchar(); exit(1); }
	for (k=0; k<nt; k++) { t=te[k]; s=0.0; r=s;
	for (j=0; j<nvm; j++) { s=s+vm[j]*pow(t, r); r=r+hf; } 
	ma[k]=s; //if (!vy) cout << "te_s = " << t << "\tt_c = " << s << "\t"; if (vy==1) cout << "te_s = " << t << "\tt_h = " << s << "\t"; if (vy==2) cout << "te_s = " << t << "\tq = " << s << "\t"; if (vy==3) cout << "te_s = " << t << "\tktp = " << s << "\t";if (vy==4) cout << "te_s = " << t << "\tts = " << s << "\t";
	} //cout << endl;
	mu[vy]=ma; if (vm) delete[]vm; } po=muv[q]; if (po) delete[]po; if (muv) delete[]muv;
	k=1; po=new double[k]; if (!po) { cout << "No memory" << endl; k=getchar(); exit(1); }
	s=0.0; for (k=0; k<nt; k++) s=s+hf; k=0; po[k]=s; mu[q]=po; //cout << "nt = " << nt << "\t";
	muv=new double*[cedu]; if (!muv) { cout << "No memory" << endl; k=getchar(); exit(1); }
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
		k=cedumi; muv=new double*[k]; if (!muv) { cout << "No memory" << endl; k = getchar(); exit(1); }
		muv=arrTem_Netzsch(muv);
		k=0; po=muv[k]; k++; mt=muv[k];
		k=0; s=e; nf207=po[k]; while (s<nf207) { s=s+hf; k++; } nk=k; //cout << "nk = " << nk << "\t"; //nk=8 - длина массива ktpv
		ktpv=arrKTP_Netzsch(); 
		if (muv) delete[]muv; if (po) delete[]po;
		k=rmu; muv=new double*[k]; k=0; if (!muv) { cout << "No memory" << endl; k=getchar(); exit(1); }
		muv=opredTempHolGor(ktpv, mt, nk, nn, h, k, muv, rmu, cem, efte, dmko); //cem - длина массива efte //for (k=0; k<f; k++) { po=muv[k]; delete[]po; }
		if (ktpv) delete[]ktpv; if (mt) delete[]mt;
		k=0; thv=muv[k]; k++; tgv=muv[k]; k++; qov=muv[k]; k++; ktpv=muv[k]; k++; 
		tsv=muv[k]; k++; po=muv[k]; if (muv) delete[]muv;
		q=0; nf207=po[q]; s=0.0; while (s<nf207) { s=s+hf; q++; } qn=q; 
		koefh=new double[dmk]; koefc=new double[dmk]; koefq=new double[dmk]; koeft=new double[dmk]; koefs=new double[dmk];
		if ((!koefh) || (!koefc) || (!koefq) || (!koeft) || (!koefs)) { cout << "No memory" << endl; k=getchar(); exit(1); }
		for (k=0; k<dmk; k++) { koefh[k]=0.0; koefc[k]=0.0; koefq[k]=0.0; koeft[k]=0.0; koefs[k]=0.0; }
		koefc=koefPribSha(thv, tsv, qn, koefc, "No memory"); koefh=koefPribSha(tgv, tsv, qn, koefh, "No memory"); 
		koefq=koefPribSha(qov, tsv, qn, koefq, "No memory"); koeft=koefPribSha(ktpv, tsv, qn, koeft, "No memory"); 
		koefs=koefPribSha(tsv, tsv, qn, koefs, "No memory"); 
		if (thv) delete[]thv; if (tgv) delete[]tgv; if (qov) delete[]qov; 
		if (ktpv) delete[]ktpv; if (tsv) delete[]tsv; if (po) delete[]po; 
	}
	k=1; po=new double[k]; if (!po) { cout << "No memory" << endl; k=getchar(); exit(1); }
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
		if (!m) { cout << "No memory" << endl; k=getchar(); exit(1); }
		for (k=0; k<np; k++) { s=s+ma[k]*x[k]; t=t+x[k]; } s=s/t;
		k=0; m[k]=s; k++;
		p=3; m[k]=ma[p]; k++;
		p=4; q=5; m[k]=(ma[p]*x[p]+ma[q]*x[q])/(x[p]+x[q]); //600, 800 и 1000 °C
		if (ma) delete[]ma; return m;
	}
	else return ma;
}
double **arrTem_Netzsch(double **mu)
{
	int k = 0, le = 8; double *tem = new double[le], hf = 1e0, nf = 0.0, *po = new double[1];
	if ((!tem) || (!po)) { cout << "No memory" << endl; k = getchar(); exit(1); }
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
	if (!ktp) { cout << "No memory" << endl; k = getchar(); exit(1); }
	ktp[k] = 7.0*1e-2;   k++; ktp[k] = 8.0*1e-2;   k++; ktp[k] = 103.0*1e-3; k++; 
	ktp[k] = 15.0*1e-2;  k++; ktp[k] = 207.0*1e-3; k++; ktp[k] = 283.0*1e-3; k++; 
	ktp[k] = 373.0*1e-3; k++; ktp[k] = 477.0*1e-3; 
	return ktp;
}
double *arrKTP_2020()
{
	int k=0, n=6; double *a = new double[n]; if (!a) { cout << "No memory" << endl; k = getchar(); exit(1); }
	a[k] = 0.175566644058715; k++; a[k] = 0.176801537812368; k++; a[k] = 0.179324717653617; k++;
	a[k] = 0.211768953068592; k++; a[k] = 0.237194543621728; k++; a[k] = 0.237231989760775;
	return a;
}
double *arrTem1_2020()
{
	int k=0, n=6; double *a=new double[n]; if (!a) { cout << "No memory" << endl; k = getchar(); exit(1); }
	a[k] = 585.0; k++; a[k] = 600.0;  k++; a[k] = 585.0; k++; 
	a[k] = 800.0; k++; a[k] = 1000.0; k++; a[k] = 1000.0;
	for (k = 0; k < n; k++) a[k] = a[k] + te0;
	return a;
}
double **arrTem2_2020(double **mu)
{
	int k=1, n=6; double *a=new double[n], nf=0.0, hf=1e0, *po=new double[k];
	if ((!a) || (!po)) { cout << "No memory" << endl; k = getchar(); exit(1); }
	for (k = 0; k < n; k++) nf = nf + hf; k=0; po[k]=nf;
	a[k] = 119.75; k++; a[k] = 138.0; k++; a[k] = 129.5; k++; 
	a[k] = 200.0;  k++; a[k] = 273.0; k++; a[k] = 261.0;
	for (k=0; k<n; k++) a[k]=a[k]+te0;
	k=0; mu[k]=po; k++; mu[k]=a;
	return mu;
}
double *arrTem3_2020()
{
	int k=0, n=6; double *a=new double[n]; if (!a) { cout << "No memory" << endl; k = getchar(); exit(1); }
	a[k] = 377.0; k++; a[k] = 396.0; k++; a[k] = 383.5; k++; 
	a[k] = 548.0; k++; a[k] = 703.0; k++; a[k] = 697.25;
	for (k = 0; k < n; k++) a[k] = a[k] + te0;
	return a;
}
double *danPoTemTepl840(double *temvs, double *temvh, int n) //Засыпка плоско-параллельная, исходный
{
	double tho1=0.0, tho2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n]; 
	int k=0, p=0, q=0;
	if (!isdan) { cout << "No memory" << endl; k=getchar(); exit(1); }
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
	if (!isdan) { cout << "No memory" << endl; k=getchar(); exit(1); }
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
	if (!isdan) { cout << "No memory" << endl; k=getchar(); exit(1); }
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
	if (!isdan) { cout << "No memory" << endl; k=getchar(); exit(1); }
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
	if (!isdan) { cout << "No memory" << endl; k = getchar(); exit(1); }
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
	if (!isdan) { cout << "No memory" << endl; k=getchar(); exit(1); }
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
double *danPoTemTepl2071(double *temvs, double *tepv, int n) //Засыпка исходная, фракция 2-0,7 мм
{
	double tepv1=0.0, tepv2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n], hf=1e0; int k=0;
	if (!isdan) { cout << "No memory" << endl; k = getchar(); exit(1); }
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
	if (!isdan) { cout << "No memory" << endl; k = getchar(); exit(1); }
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
	if (!isdan) { cout << "No memory" << endl; k = getchar(); exit(1); }
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
	if (!isdan) { cout << "No memory" << endl; k=getchar(); exit(1); }
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
double **arrTemCold84(double **mu)
{
	int k=0, n=26; double *temcol=new double[n], *po=new double[1], nf=0.0, hf=1e0;
	if ((!temcol) || (!po)) { cout << "No memory" << endl; k = getchar(); exit(1); }
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
	if (!temh) { cout << "No memory" << endl; k=getchar(); exit(1); }
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
	if (!tep) { cout << "No memory" << endl; k=getchar(); exit(1); }
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
	if ((!temcol) || (!po)) { cout << "No memory" << endl; k=getchar(); exit(1); } k=0;
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
	if (!temh) { cout << "No memory" << endl; k=getchar(); exit(1); }
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
	if (!tepot) { cout << "No memory" << endl; k = getchar(); exit(1); }
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
	t=e; while (t<nf) { k++; t=t+hf; } nk=k; mui=new int[nk]; if (!mui) { cout << "No memory" << endl; k=getchar(); exit(1); }
	q=0; for (k=0; k<nk; k++) {
		x=1; for (m=0; m<w; m++) { 
			po=mu[m]; if ((po[k]<e) && (x>0)) { x=-1; break; } } 
		if (x>0) { mui[k]=k; q++; } else mui[k]=-1; } qn=q;
		if (fl>0) { m=1; po=mu[m]; for (k=0; k<nk; k++) if (po[k]>templa) { mui[k]=-1; qn--; } } //if (fl>0) { cout << "qn = " << qn << "\tnf = " << nf << "\tw = " << w << "\tn = " << n << "\t"; } 
	for (m=0; m<w; m++) {
	vm=new double[qn]; if (!vm) { cout << "No memory" << endl; k=getchar(); exit(1); } //for (k=0; k<qn; k++) vm[k]=0.0;
	q=0; po=mu[m]; for (k=0; k<nk; k++) {
		x=mui[k]; if (x>=0) { vm[q]=po[x]; q++; } } muv[m]=vm; } qn=q; qk=q; if (mui) delete[]mui; 
	nf=0.0; for (k=0; k<qn; k++) nf=nf+hf; 
	k=1; po=new double[k]; if (!po) { cout << "No memory" << endl; k=getchar(); exit(1); } k=0; po[k]=nf; muv[w]=po; //if (fl>0) { q=w; po=muv[q]; q=0; t=po[q]; cout << "r = " << t << endl; s=e; q=0; while (s<t) { s=s+hf; q++; } qk=q; q=0; po=muv[q]; for (q=0; q<qk; q++) cout << "temc = " << po[q] << "\t"; cout << endl; q=1; po=muv[q]; for (q=0; q<qk; q++) cout << "temh = " << po[q] << "\t"; cout << endl; q=2; po=muv[q]; for (q=0; q<qk; q++) cout << "qo = " << po[q] << "\t"; cout << endl; q=3; po=muv[q]; for (q=0; q<qk; q++) cout << "ktp = " << po[q] << "\t"; cout << endl; q=4; po=muv[q]; for (q=0; q<qk; q++) cout << "tems = " << po[q] << "\t"; cout << endl; q=w; po=mu[q]; q=0; t=po[q]; cout << "r = " << t << endl; s=e; q=0; while (s<t) { s=s+hf; q++; } qk=q; q=0; po=mu[q]; for (q=0; q<qk; q++) cout << "temc = " << po[q] << "\t"; cout << endl; q=1; po=mu[q]; for (q=0; q<qk; q++) cout << "temh = " << po[q] << "\t"; cout << endl; q=2; po=mu[q]; for (q=0; q<qk; q++) cout << "qo = " << po[q] << "\t"; cout << endl; q=3; po=mu[q]; for (q=0; q<qk; q++) cout << "ktp = " << po[q] << "\t"; cout << endl; q=4; po=mu[q]; for (q=0; q<qk; q++) cout << "tems = " << po[q] << "\t"; cout << endl;  cout << mu << "\t" << muv << "\t"; }
	for (k=0; k<f; k++) { vm=mu[k]; //if ((fl>0) && (vm) && (k<w)) { for (q=0; q<nk; q++) cout << "k = "<< k << "\tq = " << q << "\tp = " << vm[q] << "\t"; cout << endl; } 
	if ((fl<0) && (vm)) delete[]vm; } //if (fl>0) { cout << "qn = " << qn << "\tn = " << nf << "\t"; } 
	return muv;
}
double *koefPribSha(double *ktp, double *te, int le, double *ko, char *snome)
{
	int k, kem=3; double **A=new double*[kem], *AA, *b=new double[kem];
	double yx2=0.0, yx=0.0, p=0.0, hf=1e0, x4=0.0, x3=0.0, x2=0.0;
	double x=0.0, y=0.0, de=0.0, de1=0.0, de2=0.0, de3=0.0;
	if (A) 
	{ 
	for (k=0; k<kem; k++) 
	{ 
	AA=new double[kem]; 
	if (AA) A[k]=AA; else { cout << snome << endl; k=getchar(); exit(1); } 
	} 
	}
	else { cout << snome << endl; k=getchar(); exit(1); }
	if ((!ko) || (!b)) { cout << snome << endl; k=getchar(); exit(1); }
	for (k=0; k<le; k++) 
	{
		yx2=yx2+ktp[k]*pow(te[k], 2e0); 
		yx=yx+ktp[k]*te[k]; 
		y=y+ktp[k];
		x4=x4+pow(te[k], 4e0); 
		x3=x3+pow(te[k], 3e0); 
		x2=x2+pow(te[k], 2e0); 
		x=x+te[k];
		p=p+hf;
	} //применение метода наименьших квадратов
	k=0; b[k]=yx2; k++; b[k]=yx; k++; b[k]=y; 
	A[0][0] = x4; A[0][1] = x3; A[0][2] = x2; A[1][0] = x3; A[1][1] = x2; A[1][2] = x; A[2][0] = x2; A[2][1] = x; A[2][2] = p;
	de = A[0][0] * (A[2][2] * A[1][1] - A[2][1] * A[1][2]) - A[0][1] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]) + A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
	de1 = b[0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]) - A[0][1] * (b[1] * A[2][2] - b[2] * A[1][2]) + A[0][2] * (b[1] * A[2][1] - b[2] * A[1][1]);
	de2 = A[0][0] * (b[1] * A[2][2] - b[2] * A[1][2]) - b[0] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]) + A[0][2] * (A[1][0] * b[2] - A[2][0] * b[1]);
	de3 = A[0][0] * (A[1][1] * b[2] - A[2][1] * b[1]) - A[0][1] * (A[1][0] * b[2] - A[2][0] * b[1]) + b[0] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
	k=2; ko[k] = de1 / de; k--; ko[k] = de2 / de; k--; ko[k] = de3 / de;
	if (b) delete[]b; 
	for (k = 0; k < kem; k++) { AA = A[k]; if (AA) delete[]AA; } 
	return ko;
}
double novNapMas(int vyve, int vymave, int vmimfv, int vyfrve, int vyukve, int vysove, int vpmf, int vpkf, int cemv)
{	
	int k=0; double por;
	if (vyve==1) 
	{ //если выбран вермикулит
	double ko=1e-2, por207=66.35*ko, poris84=55.75*ko, porin84=81.53*ko;
	double poro84=86.61*ko, poro16035=84.36*ko, por16035=83.97*ko;
	if ((!vysv) || (vysv==1)) 
	{  //исходный или после повторных измерений
	if (!vyfv) 
	{ 
	if (!vpmf) por = por207; 
	else if (vpmf==1) por=por16035; 
	} //выбор пористости мелкой фракции
	else if (vyfrve==1) 
	{ 
	if (!vpkf) por = poris84; 
	else if (vpkf==1) por=porin84; 
	} //выбор пористости крупной фракции
	else if (vyfrve==2) por = poro16035; 
	}
	if (vysove==2) 
	{ //после обжига
	if ((!vyfrve) || (vyfrve==2)) por=poro16035;
	else if (vyfrve==1) por=poro84; 
	} 
	}
	return por;
}
double **opredTempHolGor(double *ektp, double *ete, int n, int l, double h0, int v, double **mu, int rmu, int ni, double *efte, int dmkoef) //моделирование процесса теплообмена в образце
{ //n - длина массива ektp, l - длина массивов temvs, qob, ktp, temvh, temvc, ni - длина efte, qob - плотность теплового потока, которую может создать лабораторная установка 
	int cemf=rmu, k=0, j=1, w=0, rm=k, nk=n; 
	double g=0.0, p=0.0, *qon=NULL, *tena=NULL, nit=1e10, hf=1e0, *po=new double[j], del=hf;
	double ep=1e-3, d=1e-4, Thna=tn, Thnac=0.0, *koeq=new double[dmkoef], ts=0.0, laef=0.0;
	double tgor=0.0, thol=0.0, tsred=0.0, etem=0.0, dt=hf, **muv=new double*[cemf], nf=0.0, s=0.0;
	if ((!koeq) || (!po) || (!muv)) { cout << "No memory" << endl; k=getchar(); exit(1); }
	muv=PoiskZavVelTem(k, muv, cemf, ni, efte, h0);
	k=2; qon=muv[k]; k=4; tena=muv[k]; k++; po=muv[k];
	k=0; nf=po[k]; w=k; g=ep; while (g<nf) { g=g+hf; w++; } //cout << "w = " << w << endl; for (k=0; k<w; k++) cout << "q ( " << k << " ) = " << qon[k] << "\tt = " << tena[k] << endl;
	for (k=0; k<dmkoef; k++) koeq[k]=0.0; koeq=koefPribSha(qon, tena, w, koeq, "No memory");
	double *qob=new double[n], *temvs=new double[n], *temvc=new double[n], *temvh=new double[n], *ktp=new double[n];
	if ((!qob) || (!temvs) || (!temvc) || (!temvh) || (!ktp)) { cout << "No memory" << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) { qob[k]=0.0; temvs[k]=0.0; temvc[k]=0.0; temvh[k]=0.0; ktp[k]=0; }
		for (k=0; k<ni; k++) {
			ts=ete[k]; g=0.0; p=g;
			for (j=0; j<dmkoef; j++) { g=g+pow(ts, p)*koeq[j]; p=p+hf; } 
			qob[k]=g; } //cout << "n = " << n << endl; for (k=0; k<n; k++) cout << "qob = " << qob[k] << "\tektp = " << ektp[k] << "\tte = " << ete[k] << endl; 
		if (koeq) delete[]koeq; for (k=0; k<cemf; k++) { po=muv[k]; delete[]po; } delete[]muv; 
		g=0.0; for (k=0; k<n; k++) {
			s=qob[k];
			if (s>ep) {
				laef=ektp[k]; p=0.0; Thnac=Thna+g*dt; del=hf; etem=ete[k]; ktp[k]=laef;
				while ((del>ep) && (p<nit)) {
					thol=Thnac+p*d; //Tg - массив температур горячих стенок, Th - массив температур холодных стенок
					tgor=thol+qob[k]*h0/laef;
					del=fabs(2e0*etem-(thol+tgor));
					p=p+hf; } 
				g=g+hf; }
			else { thol=0.0; tgor=0.0; qob[k]=0.0; ktp[k]=0.0; } tsred=(tgor+thol)/2e0; 
				temvc[k]=thol; temvh[k]=tgor; temvs[k]=tsred; }
	g=0.0; for (k=0; k<n; k++) g=g+hf; //cout << "n = " << g << endl; for (k=0; k<n; k++) cout << "qob = " << qob[k] << "\tektp = " << ektp[k] << "\tte = " << ete[k] << "\tth = " << temvh[k] << "\ttc = " << temvc[k] << endl;
	k=1; po=new double[k]; muv=new double*[cemf]; if ((!po) || (!muv)) { cout << "No memory" << endl; k=getchar(); exit(1); } k=0; po[k]=g;
	muv[k]=temvc; k++; muv[k]=temvh; k++; muv[k]=qob; k++; muv[k]=ktp; k++; muv[k]=temvs; k++; muv[k]=po; for (k=0; k<cemf; k++) mu[k]=muv[k]; //for (k=0; k<n; k++) cout << "tc = " << temvc[k] << "\ttg = " << temvh[k] << "\tts = " << temvs[k] << endl; 
	k=0; j=1; mu=vydelPol(k, n, muv, mu, rmu, j); //cout << muv << "\t" << mu << "\t"; //temvc=mu[k]; k++; temvh=mu[k]; k++; qob=mu[k]; k++; ktp=mu[k]; k++; temvs=mu[k]; k++; po=mu[k]; k=0; g=po[k]; nf=ep; while (nf<g) { nf=nf+hf; k++; } nk=k; for (k=0; k<nk; k++) cout << "tc = " << temvc[k] << "\ttg = " << temvh[k] << "\tts = " << temvs[k] << endl; for (k=0; k<cemf; k++) { po=muv[k]; if (po) delete[]po; }
	if (temvc) delete[]temvc; if (temvh) delete[]temvh; //if (qob) delete[]qob; 
	if (ktp) delete[]ktp; if (temvs) delete[]temvs; if (po) delete[]po; if (muv) delete[]muv; 
	return mu;
}
double **PoiskZavVelTem(int v, double **mu, int rmu, int n, double *efte, double h)
{
	int k=1, f=rmu, j=0, q=j, w=j, b=10, l=0; double *poin=new double[k], cf=0.0, hf=1e0, nf=cf, tf=cf, e=1e-1, cfp=cf, cft=cf;
	double *temvht=NULL, *temvct=NULL, *tepvt=NULL, *po=NULL, *temhq=NULL, *ts=NULL, ***muu=NULL;
	double *temcq=NULL, *temvs=NULL, *tepv=NULL, *cemt=NULL, *ktpq=NULL, *ktp=NULL, **muv=NULL;
	int vyfrve=0, vysove=0, vymivmf=1, vyukve=1, c=0, nnvyfv=2, nnvysv=3, nnvmivmf=2;
	int nnvyuv=2, kvf=0, jvsv=0, qvmi=0, qvuv=0, d=0;
	double *nvyfv=new double[nnvyfv], *nvysv=new double[nnvysv], *nvmivmf=new double[nnvmivmf];
	double *nvyuv=new double[nnvyuv]; temhq=new double[n]; temcq=new double[n]; temvs=new double[n]; 
	tepv=new double[n]; cemt=new double[n]; ktpq=new double[n]; ktp=new double[n]; muu=new double**[b]; 
	if ((!nvyfv) || (!nvysv) || (!nvmivmf) || (!nvyuv) || (!temhq) || (!temcq) || (!temvs) || (!tepv) || (!cemt) || (!ktpq) || (!ktp) || (!muu)) 
	{ cout << "No memory" << endl; k=getchar(); exit(1); } 
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
					muv=new double*[f]; if (!muv) { cout << "No memory" << endl; k=getchar(); exit(1); }
					muv=napMasEKTPVer(vyfrve, vysove, vymivmf, efte, k, h, vyukve, n, k, muv, f, dmko); muu[b]=muv; b++; if (vysove>0) break; } }
			else if (vyfrve==1) {
				for (qvuv=0; qvuv<nnvyuv; qvuv++) {
					vyukve=nvyuv[qvuv];
					muv=new double*[f]; if (!muv) { cout << "No memory" << endl; k=getchar(); exit(1); }
					muv=napMasEKTPVer(vyfrve, vysove, vymivmf, efte, k, h, vyukve, n, k, muv, f, dmko); muu[b]=muv; b++; } } } } 
	for (j=0; j<b; j++) { muv=muu[j];
	d=0; temvct=muv[d]; d++; temvht=muv[d]; d++; tepvt=muv[d]; d++; ktp=muv[d]; d++; ts=muv[d]; d++; po=muv[d]; k=0; nf=po[k]; 
					cf=e; while (cf<nf) { cf=cf+hf; k++; } q=k; 
					d=1; k=0; cfp=tepvt[k]; for (k=d; k<q; k++) { cft=tepvt[k]; if ((cft<=cfp) && (d>0)) { d=-1; break; } cfp=cft; } //{ for (k=0; k<q; k++) cout << "c = " << c << "\tqo = " << tepvt[k] << "\tts = " << ts[k] << "\t"; cout << endl; } //vyvodfile(ts, q, k, nf, "C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\tmp.txt");
						if (d>0) { for (w=0; w<q; w++) { 
							tf=ts[w]; cf=tepvt[w]*tf;
							for (k=0; k<n; k++) {
							if (fabs(tf-efte[k])<=hf) {
								temvs[k]=temvs[k]+tf;
								temhq[k]=temhq[k]+tf*temvht[w]; 
								temcq[k]=temcq[k]+tf*temvct[w]; 
								ktpq[k]=ktpq[k]+tf*ktp[w];
								tepv[k]=tepv[k]+cf; 
								cemt[k]=cemt[k]+hf; //cout << "tf = " << tf << "\tptp = " << tepvt[k] << "\t";
								break; } } } } c++; }
	b=c; for (k=0; k<b; k++) { muv=muu[k]; for (j=0; j<f; j++) { po=muv[j]; delete[]po; } if (muv) delete[]muv; } //for (k=0; k<n; k++) cout << "cemt = " << cemt[k] << "\t";
		for (k=0; k<n; k++) {
			cf=temvs[k]; if (cf>e) {
				temhq[k]=temhq[k]/cf; temcq[k]=temcq[k]/cf; ktpq[k]=ktpq[k]/cf; tepv[k]=tepv[k]/cf; }
			cf=cemt[k]; if (cf>e) temvs[k]=temvs[k]/cf; else temvs[k]=0.0; }
	cf=0.0; for (k=0; k<n; k++) cf=cf+hf;
	k=0; poin[k]=cf;
	mu[k]=temcq; k++; mu[k]=temhq; k++; mu[k]=tepv; k++; 
	mu[k]=ktpq; k++; mu[k]=temvs; k++; mu[k]=poin; //for (k=0; k<n; k++) cout << "cemt = " << cemt[k] << "\tth = " << temhq[k] << "\ttc = " << temcq[k] << "\tq = " << tepv[k] << "\tktp = " << ktpq[k] << "\tts = " << temvs[k] << endl;
	if (nvyfv) delete[]nvyfv; if (nvysv) delete[]nvysv;
	if (nvmivmf) delete[]nvmivmf; if (nvyuv) delete[]nvyuv;
	if (muu) delete[]muu; if (cemt) delete[]cemt;
	return mu; 
}