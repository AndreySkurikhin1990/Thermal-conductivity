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
#define vmivmf 0 //выбор метода измерений для фракции 2-0,7 мм
#define vyfv 1 //выбор фракции: 0 - фракция 2-0,7 мм, 1 - фракция 8-4 мм
#define vyuv 2 //выбор укладки: 1 - плоскопараллельная, 2 - вертикальная
#define vysv 2 //выбор состояния: 0 - исходное, 1 - после повторных измерений, 2 - после прокаливания при 1000 °С
using namespace std;
const double pi = acos(-1.0), enu=1e-3;
const int dsov=60, dmkoscrt=14;
const char snmnr[]="No memory!";
const int dsv=50, dmkooscv=14;
const double ssiv=368e-3, salv=132e-3, smgv=217e-3, detev=1e2, tev0=273.15;
const double tnoscv=3e2; //dkoscv - максимальное отклонение от экспериментальных данных в литературных источниках
const double dtoscv=1e2, tnacv=2e2;
const double epsi=1e-15, pksvv=0.0, dkoalv=0.444905; //dkoalv - определение КП из КО
//----------
int dlar=1;
char *sndov, *ssdov, *skpovsk, *sdvovsk, *szfovk, *snmov;
char *sfnoov, *sppovv, *sppovvk, *sppovs, *sppovsk;
char *skpovs, *sdvovs, *szfov, *sdvovv, *sdvovvk, *skpovv;
char *skpovvk, *sdvovi, *sppovi, *skpovi, *sdvovik, *sppovik, *skpovik;
int dkoscvl=6, cemv=11; /*выбор температуры вермикулита*/
char *snv=NULL, *ssv=NULL, *skptv=NULL, *svsv=NULL, *snmv=NULL, *sfnov=NULL;
char *sfatv=NULL, *sfov=NULL, *svfdv=NULL, *svfdvu=NULL, nfv='D';
double *ktrv=NULL, *dkoscvm=NULL, *dkoscvt=NULL, *etev;
double dkospv=1.3002913762; //dkospv - дополнительный коэффициент для КО с учетом пористой структуры;
double *tkuscv, *kuscv, tnav=tnacv+tev0;
//------
double *koefoslab(double, double, double, double *, int);
double epsisred(double);
double trapz(double *, double *, int);
double *Kramers_Kronig(double *, int);
double *dliny_voln();
double *epsilnu(double *, int);
double *Koef_Pogl();
double *PokazPrelomAl(double *, double *, double *, double, int);
double *PokazPrelomPribl(double *, double *, double *, int);
void zapisvfile(double *, int, char *);
double *dliny_voln_sham();
double *Koef_Pogl_sha();
double *ronusha(double *, int);
void napstrdir();
double *reshMetKram(double **, double *, int);
double *oslaintestepcher(int, int);
double *RasKorOV(double, double, double, double, double, double);
double UsredMasOV(double **, int, int);
double KoefPoglRossel(double, double *, double *, double *, int);
double nsredPlank2(double *, double *, double, int);
void osvpamov();
double usredVelichPlank(double *, double *, double *, double, int);
double *KoefPoglRosselNac(double *, int, int, double, double, double, double *, double *, int, double, double, double *, double *, int);
double PoiskReflPhi(double, double, double, double);
double ReflSredVer(double);
double PoiskTetaShtrNach(double, double, double, double);
double PraCha10(double, double, double, double);
double *PoiSpeKoeOtr(double *, int);
double provUgla(double);
double **polmat(double **, double *, int, int);
double vychopred(double  **, int);
double **osvpam(double **, int);
double opredKTPTKTochSha(double *, double *, double, int);
double *izmMasChast(double *, double, int);
void initarrver(double, double, double);
//------------
void vyzovKoefPogl()
{ int idve=1, k;
initarrver(smgv, ssiv-pksvv, salv+pksvv); 
ktrv=KoefPoglRosselNac(etev, idve, cemv, smgv, ssiv-pksvv, salv+pksvv, dkoscvm, dkoscvt, dkoscvl, dkospv, dkoalv, kuscv, tkuscv, dmkooscv);
//k=getchar();
}	
double *KoefPoglRosselNac(double *tem, int ident, int arle, double wmgo, double wsio, double walo, double *dkoscem, double *dkoscet, int dkoscel, double dkosups, double dkokpiko, double *dkusctm, double *tdkusctm, int dkusctl)
{	napstrdir(); int k, j, dm;
	double *dl=NULL, *alfs=NULL, *npp=NULL, *kpr=new double[arle], krpk, temk;
	double *ktr=new double[arle], dkusct, dkosce, dko2=dkosups, dko4=dkokpiko; //КТП по Росселанду
	if (ident==1) { //для вермикулита
		dl=dliny_voln(); alfs=Koef_Pogl(); npp=Kramers_Kronig(dl, dlar); dm=dlar; } 
	for (k=0; k<dm; k++) alfs[k]=alfs[k]*dko2*dko4; cout << "dm = " << dm << endl;
	double *GnpT=new double[arle], sigma=5.67e-8, *npT=new double[arle]; //производная dn/dT, функция <n>(T)
	if ((!kpr) || (!GnpT) || (!ktr)) { cout << snmov << endl; k=getchar(); exit(1); }
	for (k=0; k<arle; k++) {
	npT[k]=usredVelichPlank(dl, npp, npp, tem[k], dm); 
	cout << "al_sred = " << usredVelichPlank(dl, alfs, npp, tem[k], dm) << "\t"; 
	temk=tem[k];
	dkusct=opredKTPTKTochSha(dkusctm, tdkusctm, temk, dkusctl);
	dkosce=opredKTPTKTochSha(dkoscem, dkoscet, temk, dkoscel);
	if ((dkosce<=0.0) || (dkosce>1e0)) dkosce=1e0; if ((dkusct<=0.0) || (dkusct>1e0)) dkusct=1e0; 
	for (j=0; j<dm; j++) alfs[j]=alfs[j]*dkusct*dkosce;
	krpk=KoefPoglRossel(temk, dl, alfs, npp, dm); kpr[k]=krpk; 
	for (j=0; j<dm; j++) alfs[j]=alfs[j]/dkusct/dkosce; 
	cout << "\ttemp = " << temk << "\tKoef_pogl_Ross = " << krpk << "\tdkod = " << dkusct << endl; } 
	for (k=1; k<arle; k++) GnpT[k-1]=(npT[k]-npT[k-1])/(tem[k]-tem[k-1]); k=arle-1; GnpT[k]=2e0*GnpT[k-1]-GnpT[k-2];
	for (k=0; k<arle; k++) { 
	temk=tem[k]; krpk=kpr[k];
    krpk=fabs(8e0*npT[k]*sigma*pow(temk,3e0)/(3e0*krpk));
    krpk=krpk*fabs(2e0*npT[k]+temk*GnpT[k]);
	ktr[k]=krpk; 
	cout << "ktrv = " << krpk << endl; } 
	zapisvfile(ktr, cemv, svfdv);
	delete []dl; delete []npp; delete []alfs; delete []kpr; delete []GnpT; delete []npT;
	osvpamov(); return ktr; }
double *koefoslab(double wmg, double wsi, double wal, double *tere, int n)
{ int lt=n, k; double wo=wmg+wsi+wal, kmg=0.0, kal=0.0, ksi=0.0, *kuo=new double[lt];
if (!kuo) { cout << snmov << endl; k=getchar(); exit(1); }
double *mgo, *alo, *sio; 
sio=oslaintestepcher(lt,0);
alo=oslaintestepcher(lt,1);
mgo=oslaintestepcher(lt,2);
for (k=0; k<n; k++) {
kmg=sio[k]; kal=alo[k]; ksi=sio[k];
if ((kmg<0.0) || (kmg>1.0)) kmg=1e0; if ((kal<0.0) || (kal>1.0)) kal=1e0; if ((ksi<0.0) || (ksi>1.0)) ksi=1e0;
if (fabs(wo)>epsi) kuo[k]=(kmg*wmg+kal*wal+ksi*wsi)/wo; else kuo[k]=0.0; }
delete []mgo; delete []alo; delete []sio; return kuo; }
double epsisred(double T)
{  napstrdir();
double *dv=dliny_voln(), *npp=Kramers_Kronig(dv,dlar);
double *epsil=epsilnu(npp,dlar), epssr=usredVelichPlank(dv, epsil, npp, T, dlar); //cout << "eps = " << epssr << "\t";
delete []dv; delete []npp; delete []epsil;
osvpamov(); return epssr; }
double reflver(double T)
{ napstrdir(); 
double *dv=dliny_voln(), *npp=Kramers_Kronig(dv,dlar), *ronu=ronusha(npp,dlar), t;
t=usredVelichPlank(dv, ronu, npp, T, dlar);
delete []dv; delete []npp; delete []ronu;
osvpamov(); return t; }
double *Kramers_Kronig(double *dv, int lear)
{	int k=0; ofstream fo; fo.open(sppovv,ios_base::out | ios_base::trunc); 
if (!fo.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
	double c0=299792458.0, *nu=new double[lear], *nus, *nnus, *nnu;
	double *alsr=Koef_Pogl(); if (!nu) { cout << snmov << endl; k=getchar(); exit(1); } 
	for (k=0; k<lear; k++) nu[k]=2e0*pi*c0/dv[k];
nus=izmMasChast(nu, enu, lear);
nnus=PokazPrelomAl(nu, nus, alsr, c0, lear);
nnu=PokazPrelomPribl(nu, nus, nnus, lear); 
for (k=0; k<lear; k++) { nnu[k]=nnu[k]-1e0+1.543; fo << nnu[k] << endl; } fo.close();
delete []nu; delete []nnus; delete []alsr;
return nnu; }
double *ronusha(double *npp, int lear)
{ int k; double *ronu=new double[lear], n; if (!ronu) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<lear; k++) {
	n=fabs(npp[k]);
    ronu[k]=0.5+((n-1.0)*(3.0*n+1.0))/(6.0*pow(n+1.0,2.0))-(2.0*pow(n,3.0)*(pow(n,2.0)+2.0*n-1.0))/((pow(n,2.0)+1.0)*(pow(n,4.0)-1.0)); 
    ronu[k]=ronu[k]+(8.0*pow(n,4.0)*(pow(n,4.0)+1.0)*log(n))/((pow(n,2.0)+1.0)*pow((pow(n,4.0)-1.0),2.0));
    ronu[k]=ronu[k]+pow(n,2.0)*pow((pow(n,2.0)-1.0),2.0)*log(fabs((n-1.0)/(n+1.0)))/pow((pow(n,2.0)+1.0),3.0); }
return ronu; }
double *epsilnu(double *npp, int lear)
{ int k; double eps, n, *epsarr=new double[lear]; 
if (!epsarr) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<lear; k++) {
	n=fabs(npp[k]);
	eps=(4.0*n+2.)/3.0/pow((n+1.0),2.0);
	eps=eps+2.0*pow(n,3.0)*(pow(n,2.0)+2.0*n-1.0)/(pow(n,2.0)+1.0)/(pow(n,4.0)-1.0);
	eps=eps-8.0*pow(n,4.0)*(pow(n,4.0)+1.0)*log(n)/(pow(n,2.0)+1.0)/pow((pow(n,4.0)-1.0),2.0);
	eps=eps-pow(n,2.0)*log((n-1.0)/(n+1.0))*pow((pow(n,2.0)-1.0),2.0)/pow((pow(n,2.0)+1.0),3.0); 
	epsarr[k]=eps; }
return epsarr; }
double *PokazPrelomAl(double *nu, double *nus, double *ka, double vl, int arle)
{ int k, j, p=arle-1, q; double dkpdo, podln, *fn=new double[arle], *np=new double[arle]; 
if ((!fn) || (!np)) { cout << snmov << endl; k=getchar(); exit(1); } 
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
double *izmMasChast(double *nu, double epnu, int chel)
{ double *nus=new double[chel]; int k;
if (!nus) { cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<(chel-1); k++) nus[k]=epnu*(nu[k+1]-nu[k])+nu[k];
nus[chel-1]=2e0*nus[chel-2]-nus[chel-3];
return nus; }
double *PokazPrelomPribl(double *nu, double *nus, double *na, int arle)
{ double *nnu=new double[arle], de; int k, p=arle-1;
if (!nnu) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<p; k++) {
    de=(nu[k]-nus[k])*(na[k+1]-na[k])/(nus[k+1]-nus[k]);
    nnu[k]=de+na[k]; }
nnu[p]=(nu[p]-nus[p])*(na[p]-na[p-1])/(nus[p]-nus[p-1])+na[p];
return nnu; }
double *Koef_Pogl()
{	ifstream fin; double *kp=new double[dlar], p; int k=0, q=100; char *s=new char[q]; 
	if ((!s) || (!kp)) { cout << snmov << endl; k=getchar(); exit(1); }
fin.open(skpovv,ios_base::in); for (k=0; k<q; k++) s[k]='\0';
if (!fin.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
k=0; while ((!fin.eof()) && (k<dlar)) { fin.getline(s,q,'\n'); p=atof(s); kp[k]=p; k++; }
fin.close(); delete []s; return kp; }
double *dliny_voln()
{ ifstream fin; double p, *dv; int k, q=100, lear; char *s=new char[q]; 
if (!s) { cout << snmov << endl; k=getchar(); exit(1); } for (k=0; k<q; k++) s[k]='\0';
fin.open(sdvovv,ios_base::in); 
if (!fin.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
k=0; while (!fin.eof()) { fin.getline(s,q,'\n'); p=atof(s); k++; } 
fin.clear(); fin.seekg(0);
lear=k; dlar=lear; dv=new double[lear]; 
if (!dv) { cout << snmov << endl; k=getchar(); exit(1); } 
k=0; while (!fin.eof()) { fin.getline(s,q,'\n'); p=atof(s); dv[k]=p;  k++; }
fin.close();
for (k=0; k<lear; k++) dv[k]=(1e-2)/dv[k];
delete []s; return dv; }
void zapisvfile(double *vyvod, int dlina, char *nazf)
{ 	if (dlina>0) { int k; FILE *fo = fopen(nazf, "a+"); 
	if (!fo) { cout << sfnoov << endl; k=getchar(); exit(1); } 
	fprintf(fo,"vmivmf = %d\n",vmivmf);
	fprintf(fo,"vyfv = %d\n",vyfv);
	fprintf(fo,"vyuv = %d\n",vyuv);
	fprintf(fo,"vysv = %d\n",vysv);
	for (k=0; k<dlina; k++) fprintf(fo,"%0.20lf\n",vyvod[k]); fprintf(fo,"%\n"); 
	fclose(fo); } }
double trapz(double *arx, double *ary, int lenarr)
{ int k; double s=0.0;
for (k=1; k<lenarr; k++) s=s+(ary[k]+ary[k-1])*(arx[k]-arx[k-1])/2e0;
return s; }
double *oslaintestepcher(int dm, int vy) //рассчитывает ослабление интегральной степени черноты для SiO2, Al2O3, MgO
{	int n=2, k, j, w=0, p=6, q=10, m, l; double **scv=new double*[p], *sc=new double[dm], **hsv=new double*[p];
	double *hs=new double[n], **kor=new double*[q], t0, t1, *vscs=new double[dm], *vscm=new double[dm], *vsca=new double[dm], *po;
if ((!scv) || (!sc) || (!hsv) || (!hs) || (!kor)) { cout << snmov << endl; k=getchar(); exit(1);}
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
	if ((!sc) || (!hs)) { cout << snmov << endl; k=getchar(); exit(1); }
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
if ((!sc) || (!hs)) { cout << snmov << endl; k=getchar(); exit(1);}
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
if ((!sc) || (!hs)) { cout << snmov << endl; k=getchar(); exit(1);}
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
if ((!sc) || (!hs)) { cout << snmov << endl; k=getchar(); exit(1); }
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
if ((!sc) || (!hs)) { cout << snmov << endl; k=getchar(); exit(1);}
k=0;
sc[k]=1.0000; k++; sc[k]=984e-3; k++; sc[k]=953e-3; k++; sc[k]=917e-3; k++; sc[k]=854e-3; k++; 
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
delete []Ibc; delete []Ibz; delete []dv;
return (1e0/rs); }
double nsredPlank2(double *dv, double *npp, double tem, int n)
{ int k; double *npp2=new double[n], t; 
if (!npp2) {cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) npp2[k]=pow(npp[k],2e0); 
t=usredVelichPlank(dv, npp2, npp, tem, n);
delete []npp2; return t; }
double usredVelichPlank(double *dv, double *uv, double *npp, double tem, int n)
{ double PP=6.6260755e-34, PB=1.380658e-23, c0=299792458.0, vl, c1, c2, lambda; 
int k; double *Ib, *Ibn, *dl, nc, nz;
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
delete []Ibn; delete []Ib; delete []dl;
return (nc/nz); }
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
double Determinant(double **mas, int m) { // Рекурсивное вычисление определителя
  int i;
  double **p, k=1., d=0.0, *po;
  p=new double*[m]; if (!p) { cout << snmnr << endl; i=getchar(); exit(1); }
  for (i=0; i<m; i++) { po=new double[m]; if (!po) { cout << snmnr << endl; i=getchar(); exit(1); } p[i]=po; }
  d=0;
  if (m<1) cout << "Определитель вычислить невозможно!";
  if (m==1) { d = mas[0][0]; return (d); } 
  if (m==2) { d = mas[0][0]*mas[1][1]-mas[1][0]*mas[0][1]; return(d); } 
  if (m>2) {
	for (i=0; i<m; i++) {
		p=GetMatr(mas, p, i, 0, m); 
		d=d+k*mas[i][0]*Determinant(p, m-1); //разложение по строкам по первому столбцу
		k=-k; } } //(-1) в степени i
for (i=0; i<m; i++) { po=p[i]; delete []po; } return (d); }
double *reshMetKram(double **mat, double *ssc, int rm)
{ int n=rm, j;
double d, *m=new double[n], mm, **mmat; if (!m) { cout << snmnr << endl; j=getchar(); exit(1); } //PrintMatr(mat, rm);
d=vychopred(mat,n); for (j=0; j<n; j++) { mmat=polmat(mat,ssc,j,n); mm=vychopred(mmat,n); 
m[j]=mm/d; /*cout << "x ( " << j << " ) = " << m[j] << endl;*/
mmat=osvpam(mmat,n); } 
return m; }
double **osvpam(double **a, int n)
{ int j; double *b; 
for (j=0; j<n; j++) 
{ b=a[j]; delete []b; } 
return NULL; }
double **polmat(double **mat, double *b, int no, int n)
{ int k, j; double **ma=new double*[n], *m; 
if (!ma) { cout << snmnr << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) {
	m=new double[n]; if (!m) { cout << snmnr << endl; k=getchar(); exit(1); }
	for (j=0; j<n; j++) {
		if (j==no) m[j]=b[k]; else m[j]=mat[k][j];}
ma[k]=m; } 
return (ma); }
double provresh(double **A, double *x, double *xp, double *b, int n, double toch)
{ int k, j, q=0; double s, *pr=new double[n], m, t, mp;
if (!pr) { cout << snmnr << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) { s=0.0; 
for (j=0; j<n; j++) s=s+A[k][j]*x[j]; pr[k]=s-b[k]; //cout << "pr = " << pr[k] << "\t"; 
} //cout << endl; 
m=0.0; for (j=0; j<n; j++) { t=fabs(pr[k]); if (t>m) { m=t; q=k; } } 
mp=0.0; for (k=0; k<n; k++) { t=fabs(xp[k]-x[k]); if (t>mp) { mp=t; } }
delete []pr; //cout << "pr = " << m << endl; //return mp; 
return m; } 
double vychopred(double  **mas, int m) {
  double d=Determinant(mas, m); //Вычисление определителя
  return d; }
void vyvodmatr(double **a, double *b, int n) // Вывод системы уравнений
{ int i, j; 
for (i=0; i<n; i++) 
{ for (j=0; j<n; j++) 
{ cout << a[i][j] << "\t";} 
cout << " --- " << b[i] << endl; } }
double *gauss(double **a, double *y, int n, double *xs) 
{ double maxi, temp=0.0, *x=new double[n]; int k=0, index=0, i=0, j=0, q;
if (!x) { cout << snmnr << endl; i=getchar(); exit(1); } 
for (j=0; j<n; j++) x[j]=xs[j];
k=0; while (k<n) { // Поиск строки с максимальным a[i][k]
    maxi=fabs(a[k][k]); index=k;
    for (i=k+1; i<n; i++) if (fabs(a[i][k])>maxi) { maxi=fabs(a[i][k]); index=i; } 
    if (!(fabs(maxi))) { for (q=0; q<n; q++) x[q]=0.0; return x; //cout << "Resheniye poluchit nevozmozhno iz-za nulevogo stolbtsa " << index << " matritsy A" << endl; k=getchar(); exit(1); 
	} //нет ненулевых диагональных элементов
    for (j=0; j<n; j++) { temp=a[k][j]; a[k][j]=a[index][j]; a[index][j]=temp; } temp=y[k]; y[k]=y[index]; y[index]=temp; //Перестановка строк
    for (i=k; i<n; i++) {
      temp=a[i][k]; if (!(fabs(temp))) continue; // для нулевого коэффициента пропустить
      for (j=0; j<n; j++) if (fabs(temp)>0.0) a[i][j]=a[i][j]/temp; else a[i][j]=0.0;
	  if (fabs(temp)>0.0) y[i]=y[i]/temp; else y[i]=0.0;
      if (i==k) continue; // уравнение не вычитать само из себя
      for (j=0; j<n; j++) a[i][j]=a[i][j]-a[k][j]; y[i]=y[i]-y[k]; } k++; }
for (k=n-1; k>=0; k--) { // обратная подстановка
	x[k]=y[k]; for (i=0; i<k; i++) y[i]=y[i]-a[i][k]*x[k]; } 
return x; }
double *reshMetGau(double **a, double *y, int n, double *x) 
{ int i, j; double **mt=new double*[n], *b, *xs;
if (!mt) { cout << snmnr << endl; i=getchar(); exit(1); } 
for (j=0; j<n; j++) y[j]=-y[j];
for (i=0; i<n; i++) { b=new double[n]; if (!b) { cout << snmnr << endl; i=getchar(); exit(1); } 
for (j=0; j<n; j++) b[j]=a[i][j]; mt[i]=b; }  
b=new double[n]; if (!b) { cout << snmnr << endl; i=getchar(); exit(1); } 
for (i=0; i<n; i++) b[i]=y[i]; //PrintMatr(a,n); 
xs=gauss(a,y,n,x); 
for (j=0; j<n; j++) x[j]=xs[j];
for (i=0; i<n; i++) { for (j=0; j<n; j++) a[i][j]=mt[i][j]; y[i]=b[i]; }
delete []b; delete []xs; for (i=0; i<n; i++) { b=mt[i]; delete []b; } //for (i=0; i<n; i++) cout << "x ( " << i << " ) = " << x[i] << endl; //i=getchar();
return x; }
void napstrdir()
{	int k; sndov=new char[dsov]; ssdov=new char[dsov]; skpovsk=new char[dsov]; skpovvk=new char[dsov];
	sdvovsk=new char[dsov]; sdvovvk=new char[dsov]; sppovvk=new char[dsov]; sppovsk=new char[dsov]; 
	szfovk=new char[dsov]; snmov=new char[dsov]; sfnoov=new char[dsov];
	skpovs=new char[2*dsov]; skpovv=new char[2*dsov]; sdvovs=new char[2*dsov]; sdvovv=new char[2*dsov];
	sppovs=new char[2*dsov]; sppovv=new char[2*dsov]; szfov=new char[2*dsov];
	sdvovik=new char[dsov]; sppovik=new char[dsov]; skpovik=new char[dsov];
	sdvovi=new char[2*dsov]; sppovi=new char[2*dsov]; skpovi=new char[2*dsov];
if ((!sndov) || (!ssdov) || (!skpovsk) || (!skpovvk) || (!sdvovsk) || (!sdvovvk) || (!sppovvk) || (!sppovsk) || (!szfovk) || (!snmov) || (!sfnoov)) 
	{ cout << "No_memory!" << endl; k=getchar(); exit(1); }
if ((!skpovs) || (!skpovv) || (!sdvovs) || (!sdvovv) || (!sppovs) || (!sppovv) || (!szfov) || (!sdvovi) || (!sppovi) || (!skpovi)) 
{ cout << "No_memory!" << endl; k=getchar(); exit(1);}
for (k=0; k<(2*dsov); k++) { skpovs[k]='\0'; skpovv[k]='\0'; sdvovs[k]='\0'; sdvovv[k]='\0'; sppovs[k]='\0'; sppovv[k]='\0'; szfov[k]='\0'; 
sdvovi[k]='\0'; sppovi[k]='\0'; skpovi[k]='\0'; } 
for (k=0; k<dsov; k++) { sndov[k]='\0'; ssdov[k]='\0'; skpovsk[k]='\0'; skpovvk[k]='\0'; sdvovsk[k]='\0'; sdvovvk[k]='\0'; sppovvk[k]='\0'; 
sppovsk[k]='\0'; szfovk[k]='\0'; snmov[k]='\0'; sfnoov[k]='\0'; sdvovik[k]='\0'; sppovik[k]='\0'; skpovik[k]='\0'; }
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
skpovsk[k]='K'; k++; skpovsk[k]='o'; k++; skpovsk[k]='e'; k++; skpovsk[k]='f'; k++; skpovsk[k]='f'; k++; skpovsk[k]='i'; k++;
skpovsk[k]='c'; k++; skpovsk[k]='i'; k++; skpovsk[k]='e'; k++; skpovsk[k]='n'; k++; skpovsk[k]='t'; k++; skpovsk[k]='_'; k++;
skpovsk[k]='p'; k++; skpovsk[k]='o'; k++; skpovsk[k]='g'; k++; skpovsk[k]='l'; k++; skpovsk[k]='o'; k++; skpovsk[k]='s'; k++;
skpovsk[k]='c'; k++; skpovsk[k]='h'; k++; skpovsk[k]='e'; k++; skpovsk[k]='n'; k++; skpovsk[k]='i'; k++; skpovsk[k]='y'; k++;
skpovsk[k]='a'; k++; skpovsk[k]='_'; k++; skpovsk[k]='s'; k++; skpovsk[k]='h'; k++; skpovsk[k]='a'; k++; skpovsk[k]='.'; k++;
skpovsk[k]='t'; k++; skpovsk[k]='x'; k++; skpovsk[k]='t'; k++; skpovsk[k]='\0';
k=0;
skpovvk[k]='K'; k++; skpovvk[k]='o'; k++; skpovvk[k]='e'; k++; skpovvk[k]='f'; k++; skpovvk[k]='f'; k++; skpovvk[k]='i'; k++;
skpovvk[k]='c'; k++; skpovvk[k]='i'; k++; skpovvk[k]='e'; k++; skpovvk[k]='n'; k++; skpovvk[k]='t'; k++; skpovvk[k]='_'; k++;
skpovvk[k]='p'; k++; skpovvk[k]='o'; k++; skpovvk[k]='g'; k++; skpovvk[k]='l'; k++; skpovvk[k]='o'; k++; skpovvk[k]='s'; k++;
skpovvk[k]='c'; k++; skpovvk[k]='h'; k++; skpovvk[k]='e'; k++; skpovvk[k]='n'; k++; skpovvk[k]='i'; k++; skpovvk[k]='y'; k++;
skpovvk[k]='a'; k++; skpovvk[k]='_'; k++; skpovvk[k]='v'; k++; skpovvk[k]='e'; k++; skpovvk[k]='r'; k++; skpovvk[k]='.'; k++;
skpovvk[k]='t'; k++; skpovvk[k]='x'; k++; skpovvk[k]='t'; k++; skpovvk[k]='\0';
k=0;
sdvovsk[k]='D'; k++; sdvovsk[k]='l'; k++; sdvovsk[k]='i'; k++; sdvovsk[k]='n'; k++; sdvovsk[k]='a'; k++; sdvovsk[k]='_'; k++; 
sdvovsk[k]='v'; k++; sdvovsk[k]='o'; k++; sdvovsk[k]='l'; k++; sdvovsk[k]='n'; k++; sdvovsk[k]='y'; k++; sdvovsk[k]='_'; k++; 
sdvovsk[k]='s'; k++; sdvovsk[k]='h'; k++; sdvovsk[k]='a'; k++; sdvovsk[k]='.'; k++; sdvovsk[k]='t'; k++; sdvovsk[k]='x'; k++; 
sdvovsk[k]='t'; k++; sdvovsk[k]='\0';
k=0;
sdvovvk[k]='D'; k++; sdvovvk[k]='l'; k++; sdvovvk[k]='i'; k++; sdvovvk[k]='n'; k++; sdvovvk[k]='y'; k++; sdvovvk[k]='_'; k++; 
sdvovvk[k]='v'; k++; sdvovvk[k]='o'; k++; sdvovvk[k]='l'; k++; sdvovvk[k]='n'; k++; sdvovvk[k]='_'; k++; sdvovvk[k]='v'; k++; 
sdvovvk[k]='e'; k++; sdvovvk[k]='r'; k++; sdvovvk[k]='.'; k++; sdvovvk[k]='t'; k++; sdvovvk[k]='x'; k++; sdvovvk[k]='t'; k++;
sdvovvk[k]='\0'; 
k=0;
szfovk[k]='V'; k++; szfovk[k]='y'; k++; szfovk[k]='v'; k++; szfovk[k]='o'; k++; szfovk[k]='d'; k++; szfovk[k]='v'; k++; 
szfovk[k]='F'; k++; szfovk[k]='i'; k++; szfovk[k]='l'; k++; szfovk[k]='e'; k++; szfovk[k]='.'; k++; szfovk[k]='t'; k++;
szfovk[k]='x'; k++; szfovk[k]='t'; k++; szfovk[k]='\0';
k=0;
snmov[k]='N'; k++; snmov[k]='o'; k++; snmov[k]='_'; k++; snmov[k]='m'; k++; snmov[k]='e'; k++; snmov[k]='m'; k++; 
snmov[k]='o'; k++; snmov[k]='r'; k++; snmov[k]='y'; k++; snmov[k]='!'; k++; snmov[k]='\0';
k=0;
sfnoov[k]='F'; k++; sfnoov[k]='i'; k++; sfnoov[k]='l'; k++; sfnoov[k]='e'; k++; sfnoov[k]='_'; k++; sfnoov[k]='i'; k++;
sfnoov[k]='s'; k++; sfnoov[k]='_'; k++; sfnoov[k]='n'; k++; sfnoov[k]='o'; k++; sfnoov[k]='t'; k++; sfnoov[k]='_'; k++;
sfnoov[k]='o'; k++; sfnoov[k]='p'; k++; sfnoov[k]='e'; k++; sfnoov[k]='n'; k++; sfnoov[k]='!'; k++; sfnoov[k]='\0';
k=0;
sppovsk[k]='P'; k++; sppovsk[k]='o'; k++; sppovsk[k]='k'; k++; sppovsk[k]='a'; k++; sppovsk[k]='z'; k++; sppovsk[k]='a'; k++; 
sppovsk[k]='t'; k++; sppovsk[k]='_'; k++; sppovsk[k]='p'; k++; sppovsk[k]='r'; k++; sppovsk[k]='e'; k++; sppovsk[k]='l'; k++;
sppovsk[k]='o'; k++; sppovsk[k]='m'; k++; sppovsk[k]='l'; k++; sppovsk[k]='e'; k++; sppovsk[k]='n'; k++; sppovsk[k]='_'; k++;
sppovsk[k]='s'; k++; sppovsk[k]='h'; k++; sppovsk[k]='a'; k++; sppovsk[k]='.'; k++; sppovsk[k]='t'; k++; sppovsk[k]='x'; k++; 
sppovsk[k]='t'; k++; sppovsk[k]='\0';
k=0;
sppovvk[k]='P'; k++; sppovvk[k]='o'; k++; sppovvk[k]='k'; k++; sppovvk[k]='a'; k++; sppovvk[k]='z'; k++; sppovvk[k]='a'; k++;
sppovvk[k]='t'; k++; sppovvk[k]='_'; k++; sppovvk[k]='p'; k++; sppovvk[k]='r'; k++; sppovvk[k]='e'; k++; sppovvk[k]='l'; k++;
sppovvk[k]='o'; k++; sppovvk[k]='m'; k++; sppovvk[k]='l'; k++; sppovvk[k]='e'; k++; sppovvk[k]='n'; k++; sppovvk[k]='_'; k++;
sppovvk[k]='v'; k++; sppovvk[k]='e'; k++; sppovvk[k]='r'; k++; sppovvk[k]='.'; k++; sppovvk[k]='t'; k++; sppovvk[k]='x'; k++;
sppovvk[k]='t'; k++; sppovvk[k]='\0';
k=0;
sdvovik[k]='D'; k++; sdvovik[k]='l'; k++; sdvovik[k]='i'; k++; sdvovik[k]='n'; k++; sdvovik[k]='y'; k++; sdvovik[k]='_'; k++;
sdvovik[k]='v'; k++; sdvovik[k]='o'; k++; sdvovik[k]='l'; k++; sdvovik[k]='n'; k++; sdvovik[k]='_'; k++; sdvovik[k]='i'; k++;
sdvovik[k]='t'; k++; sdvovik[k]='o'; k++; sdvovik[k]='m'; k++; sdvovik[k]='.'; k++; sdvovik[k]='t'; k++; sdvovik[k]='x'; k++;
sdvovik[k]='t'; k++; sdvovik[k]='\0';
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
skpovik[k]='.'; k++; skpovik[k]='t'; k++; skpovik[k]='x'; k++; skpovik[k]='t'; k++; skpovik[k]='\0'; k++;
strcpy(sppovv,ssdov); strcat(sppovv,sppovvk); sppovv[strlen(sppovv)+1]='\0'; 
strcpy(sppovs,ssdov); strcat(sppovs,sppovsk); sppovv[strlen(sppovs)+1]='\0';
strcpy(skpovs,ssdov); strcat(skpovs,skpovsk); skpovs[strlen(skpovs)+1]='\0'; 
strcpy(skpovv,ssdov); strcat(skpovv,skpovvk); skpovs[strlen(skpovv)+1]='\0'; 
strcpy(sdvovs,ssdov); strcat(sdvovs,sdvovsk); sdvovs[strlen(sdvovs)+1]='\0'; 
strcpy(sdvovv,ssdov); strcat(sdvovv,sdvovvk); sdvovs[strlen(sdvovs)+1]='\0';
strcpy(szfov,ssdov); strcat(szfov,szfovk); szfov[strlen(szfov)+1]='\0'; 
strcpy(sdvovi,ssdov); strcat(sdvovi,sdvovik); szfov[strlen(sdvovi)+1]='\0'; 
strcpy(sppovi,ssdov); strcat(sppovi,sppovik); szfov[strlen(sppovi)+1]='\0'; 
strcpy(skpovi,ssdov); strcat(skpovi,skpovik); szfov[strlen(skpovi)+1]='\0'; }
void osvpamov()
{ delete []sndov; delete []ssdov; delete []skpovsk; 
delete []sdvovsk; delete []szfovk; delete []snmov; 
delete []sfnoov; delete []sppovv; delete []sppovvk; 
delete []sppovs; delete []sppovsk; delete []skpovs; 
delete []sdvovs; delete []szfov; delete []sdvovv;
delete []sdvovvk; delete []skpovv; delete []skpovvk; 
delete []sdvovi; delete []sppovi; delete []skpovi; 
delete []sdvovik; delete []sppovik; delete []skpovik; }
void napstrver()
{ 	if ((!sfatv) || (!sfov) || (!snv) || (!ssv) || (!skptv) || (!svsv) || (!snmv) || (!sfnov)) { cout << "No_memory!" << endl; getchar(); exit(1); } 
int k=0; 
snv[k]='D'; k++; snv[k]=':';  k++; snv[k]='\\'; k++; snv[k]='\\'; k++; snv[k]='_'; k++; 
snv[k]='А'; k++; snv[k]='с';  k++; snv[k]='п';  k++; snv[k]='и';  k++; snv[k]='р'; k++; 
snv[k]='а'; k++; snv[k]='н';  k++; snv[k]='т';  k++; snv[k]='у';  k++; snv[k]='р'; k++; 
snv[k]='а'; k++; snv[k]='\\'; k++; snv[k]='\\'; k++; snv[k]='t';  k++; snv[k]='m'; k++; 
snv[k]='p'; k++; snv[k]='\\'; k++; snv[k]='\\'; k++; snv[k]='\0';
k=0;
ssv[k]='C';  k++; ssv[k]=':'; k++; ssv[k]='\\'; k++; ssv[k]='\\'; k++; ssv[k]='U';  k++; 
ssv[k]='s';  k++; ssv[k]='e'; k++; ssv[k]='r';  k++; ssv[k]='s';  k++; ssv[k]='\\'; k++; 
ssv[k]='\\'; k++; ssv[k]='А'; k++; ssv[k]='н';  k++; ssv[k]='д';  k++; ssv[k]='р';  k++; 
ssv[k]='е';  k++; ssv[k]='й'; k++; ssv[k]='\\'; k++; ssv[k]='\\'; k++; ssv[k]='D';  k++; 
ssv[k]='o';  k++; ssv[k]='c'; k++; ssv[k]='u';  k++; ssv[k]='m';  k++; ssv[k]='e';  k++; 
ssv[k]='n';  k++; ssv[k]='t'; k++; ssv[k]='s';  k++; ssv[k]='\\'; k++; ssv[k]='\\'; k++;
ssv[k]='_';  k++; ssv[k]='А'; k++; ssv[k]='с';  k++; ssv[k]='п';  k++; ssv[k]='и';  k++; 
ssv[k]='р';  k++; ssv[k]='а'; k++; ssv[k]='н';  k++; ssv[k]='т';  k++; ssv[k]='у';  k++; 
ssv[k]='р';  k++; ssv[k]='а'; k++; ssv[k]='\\'; k++; ssv[k]='\\'; k++; ssv[k]='t';  k++; 
ssv[k]='m';  k++; ssv[k]='p'; k++; ssv[k]='\\'; k++; ssv[k]='\\'; k++; ssv[k]='\0';
k=0;
skptv[k]='K'; k++; skptv[k]='o'; k++; skptv[k]='e'; k++; skptv[k]='f'; k++; skptv[k]='f'; k++;
skptv[k]='i'; k++; skptv[k]='c'; k++; skptv[k]='i'; k++; skptv[k]='e'; k++; skptv[k]='n'; k++; 
skptv[k]='t'; k++; skptv[k]='_'; k++; skptv[k]='p'; k++; skptv[k]='o'; k++; skptv[k]='g'; k++; 
skptv[k]='l'; k++; skptv[k]='o'; k++; skptv[k]='s'; k++; skptv[k]='c'; k++; skptv[k]='h'; k++;
skptv[k]='e'; k++; skptv[k]='n'; k++; skptv[k]='i'; k++; skptv[k]='y'; k++; skptv[k]='a'; k++; 
skptv[k]='_'; k++; skptv[k]='v'; k++; skptv[k]='e'; k++; skptv[k]='r'; k++; skptv[k]='_'; k++;
skptv[k]='T'; k++; skptv[k]='.'; k++; skptv[k]='t'; k++; skptv[k]='x'; k++; skptv[k]='t'; k++; 
skptv[k]='\0';
k=0;
svsv[k]='D'; k++; svsv[k]='o'; k++; svsv[k]='l'; k++; svsv[k]='i'; k++; svsv[k]='_'; k++; 
svsv[k]='p'; k++; svsv[k]='r'; k++; svsv[k]='o'; k++; svsv[k]='p'; k++; svsv[k]='_'; k++;
svsv[k]='V'; k++; svsv[k]='e'; k++; svsv[k]='r'; k++; svsv[k]='m'; k++; svsv[k]='i'; k++; 
svsv[k]='k'; k++; svsv[k]='-'; k++; svsv[k]=nfv; k++; svsv[k]='.'; k++; svsv[k]='t'; k++;
svsv[k]='x'; k++; svsv[k]='t'; k++; svsv[k]='\0';
k=0;
snmv[k]='N'; k++; snmv[k]='o'; k++; snmv[k]='_'; k++; snmv[k]='m'; k++; snmv[k]='e'; k++;
snmv[k]='m'; k++; snmv[k]='o'; k++; snmv[k]='r'; k++; snmv[k]='y'; k++; snmv[k]='!'; k++; snmv[k]='\0';
k=0;
sfnov[k]='F'; k++; sfnov[k]='i'; k++; sfnov[k]='l'; k++; sfnov[k]='e'; k++; sfnov[k]='_'; k++;
sfnov[k]='i'; k++; sfnov[k]='s'; k++; sfnov[k]='_'; k++; sfnov[k]='n'; k++; sfnov[k]='o'; k++;
sfnov[k]='t'; k++; sfnov[k]='_'; k++; sfnov[k]='o'; k++; sfnov[k]='p'; k++; sfnov[k]='e'; k++;
sfnov[k]='n'; k++; sfnov[k]='!'; k++; sfnov[k]='\0';
k=0;
svfdvu[k]='V'; k++; svfdvu[k]='y'; k++; svfdvu[k]='v'; k++; svfdvu[k]='o'; k++; svfdvu[k]='d'; k++; 
svfdvu[k]='v'; k++; svfdvu[k]='F'; k++; svfdvu[k]='i'; k++; svfdvu[k]='l'; k++; svfdvu[k]='e'; k++;
svfdvu[k]='.'; k++; svfdvu[k]='t'; k++; svfdvu[k]='x'; k++; svfdvu[k]='t'; k++; svfdvu[k]='\0';
for (k=0; k<(2*dsv); k++) { sfatv[k]='\0'; sfov[k]='\0'; }
strcpy(sfatv,ssv); strcat(sfatv,skptv); sfatv[strlen(sfatv)+1]='\0';
strcpy(sfov,ssv); strcat(sfov,svsv); sfov[strlen(sfov)+1]='\0'; 
strcpy(svfdv,ssv); strcat(svfdv,svfdvu); svfdv[strlen(svfdv)+1]='\0';}
void initarrver(double wmg, double wsi, double wal)
{ snv=new char[dsv]; ssv=new char[dsv]; skptv=new char[dsv]; svsv=new char[dsv]; 
snmv=new char[dsv]; sfnov=new char[dsv]; svfdvu=new char[dsv];
int j; if ((!snv) || (!ssv) || (!skptv) || (!svsv) || (!snmv) || (!sfnov) || (!svfdvu)) { cout << "No  memory!" << endl; j=getchar(); exit(1); }
for (j=0; j<dsv; j++) { snv[j]='\0'; ssv[j]='\0'; skptv[j]='\0'; svsv[j]='\0'; snmv[j]='\0'; sfnov[j]='\0'; svfdvu[j]='\0'; }
j=2*dsv; sfatv=new char[j]; sfov=new char[j]; svfdv=new char[j];
if ((!sfatv) || (!sfov) || (!svfdv)) { cout << "No memory!" << endl; j=getchar(); exit(1); }
for (j=0; j<(2*dsv); j++) { sfatv[j]='\0'; sfov[j]='\0'; svfdv[j]='\0'; }
napstrver();
dkoscvm=new double[dkoscvl]; dkoscvt=new double[dkoscvl]; int k; 
if ((!dkoscvm) || (!dkoscvt)) { cout << snmv << endl; k=getchar(); exit(1); }
double tnd=6e2, dtd=2e2, tm; dkoscvm[0]=tnd; for (k=1; k<dkoscvl; k++) dkoscvt[k]=dkoscvt[k-1]+dtd;
k=0; dkoscvt[k]=4.68; k++; dkoscvt[k]=4.69; k++; dkoscvt[k]=5.65; k++; dkoscvt[k]=13.17; k++; dkoscvt[k]=20.2; k++; dkoscvt[k]=27.81;
for (k=0; k<dkoscvl; k++) { tm=dkoscvt[k]/1e2; dkoscvt[k]=1e0-tm; cout << "dkosc = " << dkoscvt[k] << "\t";
} //cout << endl; 
etev=new double[cemv];
if (!etev) { cout << snmv << endl; j=getchar(); exit(1); }
for (j=0; j<cemv; j++) { etev[j]=0.0; }
tkuscv=new double[dmkooscv]; kuscv=new double[dmkooscv]; 
if ((!tkuscv) || (!kuscv)) { cout << snmv << endl; k=getchar(); exit(1); }
tkuscv[0]=tnoscv; for (k=1; k<dmkooscv; k++) tkuscv[k]=tkuscv[k-1]+dtoscv;
kuscv=koefoslab(wmg, wsi, wal, tkuscv, dmkooscv); 
etev[0]=tnav; for (j=1; j<cemv; j++) etev[j]=etev[j-1]+detev;
for (k=0; k<cemv; k++) cout << "te = " << etev[k] << endl; }
void osvpamver()
{ 	delete []snv; delete []ssv; delete []skptv; delete []svsv; delete []snmv; 
	delete []sfnov; delete []svfdvu; delete []sfatv; delete []sfov; 
	delete []svfdv; delete []etev; }
double opredKTPTKTochSha(double *ktptks, double *te, double temp, int ce)
{ int n=ce, f=1, p=0, k;
for (k=0; k<n; k++)
	if ((te[k]>=temp) && (f>0)) { p=k; f=0; }
	double ktp=0.0, ko=0.0; if (!p) if (!f) p=1; else { p=n-1; f=0; }
	if ((!f) && (p>0)) { ko=(ktptks[p]-ktptks[p-1])/(te[p]-te[p-1]);
	ktp=ktptks[p-1]+ko*(temp-te[p-1]); } else ktp=0.0;
return ktp; }