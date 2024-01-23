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
class KoePog { public: double alp, tem; KoePog *nex; };
//-------
double **RaschRTA(int, double, double, double, int, double, int, double, int, int, int, int, double *, double, int, double, double, 
	double, double, double, char *, double, double *, double *, int, double *, double *, int, double, double *, double *, int, double, double);
double **KoefPoglSten(int, double, int, double, double, double, double, double, int, double, double, double, int, int);
double **izmRTA(int, int, double *, double *, double *, double *, double *, double *, double *, double *, double *, char *, double **, double *, 
	double *, double *, int, double *, double *, int);
char **napolStrok(int, int);
char **napNazvFile(int, int, char *);
double F0_lamT(double, double, char *);
double *EffectTols(double *, double *, double *, double, double, int, char *);
double RaschAlphaTvKar(int, int, int, int, int, int, int, double, char *, int, double, int, double *, double *, double *, int, int, 
	int, double *, int, double *, double *, double *, double *, int, double, double, double, double, double, double, double, int, 
	double, double, double, double *, double, double, int, double *, double *, int, int, int);
double *KoefPoglRosselNac(double *, int, double, double, int, int);
double **RasPorPoRazVer(double, int, int, int, int, char *, int);
double **opredTempLPSten(double, int, double *, double *, double *, double *, int, int, double, double, double, int, char *, double);
double UtochKoefPogl(double, double, double, int, int, double);
double opredPokPreSr(double, int, int, double);
double **oprRasTemNach(int, double, double, double, double, double, double, char *);
double bbfn(double, double);
double opredKTPTKToch(double *, double *, double, int);
double ***zadrkt(int, char *, int, double *, double *, double, double, double, double, int);
double BolTochRasAlpha(int, double, double, double, double *, double *, double *, int, int, char *, double *, int, int, int, double *, 
	double *, double *, double *, int, double, double, double, double, double, double, double, int, double, double, double, double, 
	double *, double, double, int, double *, double *, int, int, int);
double opredTolTverChasObFigObPor(double, double, double, double, double, int, double, double);
double opredTempStenFragm(double, int, double *, double *, double *, double *, int, int, double, double, double, double);
double *RasIzlSerStenNacSok(double *, double *, double *, double ***, double ***, double, double, double, double, double, double, double, int, double *, 
	double *, int, double *, double *, char *, double *, int, double, double, double, double, double, double, double, int, int, int,
	double *, int, double, int, int, int, int, int);
double koefPoglSred(double, int, int, double);
double integroexponen(int, double);
double VIEW(int, int, double *, int);
double **vybFunRasPorPoRazSha(double, int, int, char *, int);
double **RasPorPoRazitom(int, char *);
double **RasPorPoRazkvi(int, char *);
double **promPreobMas(double **);
//-------
double **RaschRTA(int kost, double htk, double kta, double ktb, int izm, double hvo, int v, double ti, int c, int w, int vybves, 
	int cem, double *ete, double hk, int vybmar, double dkusce, double dkoscee, double thol, double tgor, double y0, char *snm, 
	double dkosp, double *kusc, double *tkusc, int dmkoosc, double *dkoscm, double *dkosct, int dkoscl, double dkoal, double *Rmp, 
	double *rts, int dmRef, double a, double b)
{ 
	int k=kost, j=0, q=0, x=0, kk=9; //cout << "htk = " << htk << "\t";
	double *te=NULL, **mauk=NULL, hf=1e0, *alpha=NULL, t=0.0, r=0.0, *po=NULL;
	double dko=dkoal*dkosp, *arg=NULL, hmi=1e-3, ko=1e-6;
	if (dko<0.0) dko=0.0; if (dko>hf) dko=hf;
	if (hvo>hmi) hvo=hvo*ko; if (htk>hmi) htk=htk*ko;
	if (!v) { if ((vybves>=0) && (vybves<4)) { //определение коэффициентов поглощения при температуре стенок
		mauk=KoefPoglSten(v, ti, kost, htk, hvo, dkoal, dkosp, hk, cem, thol, tgor, y0, vybves, vybmar);
		q=0; alpha=mauk[q]; q++; te=mauk[q]; }
		if (!c) { k=kost; //задание коэффициентов отражения, поглощения и пропускания: v=0, c=0
double *Ra=NULL, *R=new double[kost], *Rt=new double[kost]; 
Ra=new double[kost]; for (k=0; k<kost; k++) { t=opredKTPTKToch(Rmp, rts, te[k], dmRef);
if (t<0.0) t=0.0; if (t>hf) t=hf; Ra[k]=t; }
for (k=0; k<kost; k++) { Rt[k]=Ra[k]; R[k]=Ra[k]; }
double *Ta=new double[kost], *Tt=new double[kost], *T=new double[kost];
double *Aa=new double[kost], *At=new double[kost], *A=new double[kost];
k=3; x=k; j=38; arg=new double[k]; 
k=0; arg[k]=a; k++; arg[k]=b; k++; arg[k]=hvo;
if ((!A) || (!Aa) || (!At) || (!Ta) || (!Tt) || (!T) || (!Rt) || (!R) || (!Ra) || (!arg)) { cout << snm << endl; k=getchar(); exit(1); }
	for (k=0; k<kost; k++) {
	Ta[k]=-alpha[k]*htk; Ta[k]=exp(Ta[k]); 
	T[k]=Ta[k];
	Aa[k]=hf-Ta[k]-Ra[k]; 
	if (Aa[k]<0.0) { r=alpha[k]*htk; Aa[k]=r; Ra[k]=hf-Ta[k]-Aa[k]; 
	if (Ra[k]<0.0) Ra[k]=0.0; if (Ra[k]>hf) Ra[k]=hf; } 
	A[k]=Aa[k];
	q=3; 
	r=alpha[k]*htk; r=integroexponen(q, r); if (r<0.0) r=0.0; if (r>hf) r=hf;
	t=VIEW(j, x, arg, k); if (t<0.0) t=0.0; if (t>hf) t=hf;
	Tt[k]=r*t;
	At[k]=hf-Rt[k]-Tt[k];
	if (At[k]<0.0) { r=alpha[k]*htk; At[k]=r; Rt[k]=hf-Tt[k]-At[k]; 
	if (Rt[k]<0.0) Rt[k]=0.0; if (Rt[k]>hf) Rt[k]=hf; } }
		if (mauk) delete[]mauk; if (arg) delete[]arg;
		mauk=new double*[kk]; if (!mauk) { cout << snm << endl; k=getchar(); exit(1); }
		k=kost; q=1;
		mauk=izmRTA(k, q, R, T, A, Ra, Ta, Aa, Rt, Tt, At, snm, mauk, te, kusc, tkusc, dmkoosc, dkoscm, dkosct, dkoscl); //for (k=0; k<kost; k++) cout << "T ( " << k << " ) = " << T[k] << "\tR = " << R[k] << "\tA = " << A[k] << "\tal = " << alpha[k] << "\t";
		if (alpha) delete[]alpha; if (te) delete[]te; } } //отправляется массив температур
	if (v==1) mauk=KoefPoglSten(v, ti, kost, htk, hvo,dkoal, dkosp, hk, cem, thol, tgor, y0, vybves, vybmar); //поиск коэффициента поглощения при температуре ti
	return mauk;
}
//-----
double **KoefPoglSten(int vyb, double tl, int kost, double htk, double hvp, double dkoal, double dkosp, double hk, int cem, 
	double thol, double tgor, double tol, int vybves, int vybmar)
{
	int k=0, j=0, y=j+1, l=y+1, m=l+1, i=y, q=100, r=0, f=0;
	char **unau=napolStrok(vybves, vybmar), *snm=unau[j], *sfno=unau[y];
	char **mauk=napNazvFile(vybves, vybmar, snm);
	char *sppt=mauk[j], *skpt=mauk[y], *sdvt=mauk[l], *skpovt=mauk[m]; 
	double hf=1e0, te0=273.15, t0=22.0, te=t0+te0, dte=hf, t=0.0, *p=new double[i], dko=t;
	double a1=0.0, a2=a1, t1=a1, t2=a1, **mu=NULL, **muv=NULL, ko=1e-6, hmi=1e-3;
	if (hvp>hmi) hvp=hvp*ko; if (htk>hmi) htk=htk*ko;
	char *s=new char[q]; 
	KoePog *kp=new KoePog, *ne=NULL, *roo=kp, *pre=NULL;
	if ((!kp) || (!s)) { cout << snm << endl; k=getchar(); exit(1); } for (k=0; k<q; k++) s[k]='\0';
	ifstream fin; //cout << "skpovt = " << skpovt << "\ttl = " << tl << "\t";
	fin.open(skpovt); if (!fin.is_open()) { cout << sfno << endl; k=getchar(); exit(1); }
	if (dkosp>hf) dkosp=hf; if (dkosp<0.0) dkosp=0.0; 
	if (dkoal>hf) dkoal=hf; if (dkoal<0.0) dkoal=0.0;
	dko=dkoal*dkosp;
	roo=kp; k=0; 
	while (!fin.eof()) {
		fin.getline(s, q, '\n');
		ne = new KoePog; if (!ne) { cout << snm << endl; k=getchar(); exit(1); }
		kp->alp=atof(s); kp->tem=te; 
		kp->nex=ne; pre=kp; kp=ne; k++; te=te+dte; kp->nex=NULL; }
	if (ne) delete[]ne; kp=NULL; pre->nex=kp; if (s) delete[]s; 
	r=k; fin.close(); 
	if (!vyb) { //поиск массива коэффициентов поглощения при температурах стенок
		f=2; mu=new double*[f];
		double knat=0.0, *al=new double[kost], *teks=NULL, *xi=NULL;
		if ((!al) || (!mu)) { cout << snm << endl; k=getchar(); exit(1); } //cout << "th = " << thol << "\ttg = " << tgor << "\ttol = " << tol << "\thtk = " << htk << "\thvp = " << hvp << "\thk = " << hk << "\tks = " << kost << "\t";
		muv=oprRasTemNach(kost, htk, hvp, hk, thol, tgor, tol, snm); 
		k=0; teks=muv[k]; k++; xi=muv[k]; 
		for (k=0; k<kost; k++) {
			kp=roo->nex; pre=roo; j=0; t=teks[k]; //cout << "alp = " << pre->alp << "\ttem = " << pre->tem << "\t";
			while ((kp) && (j<r)) {
				if (kp->tem>=t) {
					a1=pre->alp; a2=kp->alp; 
					t1=pre->tem; t2=kp->tem; 
					knat=(a2-a1)/(t2-t1); //cout << "alp = " << a2 << "\ttem = " << t2 << "\tdko = " << dko << "\tknat = " << knat << "\tt = " << t << "\t";
					al[k]=a1+knat*(t-t1); al[k]=dko*al[k]; break; }
				else { pre=kp; kp=kp->nex; } j++; } } //for (k=0; k<kost; k++) cout << "al ( " << k << " ) = " << al[k] << "\tte = " << teks[k] << "\t"; k=getchar();
		if (xi) delete[]xi; if (muv) delete[]muv; if (p) delete[]p;
		k=0; mu[k]=al; k++; mu[k]=teks; }
	if (vyb==1) { //поиск конкретного коэффициента поглощения при данной температуре tl
		kp=roo->nex; j=0; pre=roo; t=0.0;
		while ((kp) && (j<r)) {
			if (kp->tem>tl) {
				a1=pre->alp; a2=kp->alp; 
				t1=pre->tem; t2=kp->tem;
				t=(a2-a1)/(t2-t1);
				t=a1+t*(tl-t1); break; }
			else { pre=kp; kp=kp->nex; } j++; }
	f=1; mu=new double*[f]; if ((!p) || (!mu)) { cout << snm << endl; j=getchar(); exit(1); } 
	j=0; p[j]=t; mu[j]=p; }
	kp=roo; while (kp) { ne=kp->nex; delete kp; kp=ne; } //удаление списка
	if (snm) delete[]snm; if (sfno) delete[]sfno; if (unau) delete[]unau;
	if (sppt) delete[]sppt; if (skpt) delete[]skpt; if (sdvt) delete[]sdvt; 
	if (skpovt) delete[]skpovt; if (mauk) delete[]mauk; //удаление имен файлов
	return mu;
}
//-----
double **izmRTA(int kost, int izm, double *Ra, double *Ta, double *Aa, double *Rb, double *Tb, double *Ab, double *Rc, 
	double *Tc, double *Ac, char *snm, double **mu, double *temp, double *kusc, double *tkusc, int dmkoosc, double *dkoscm, 
	double *dkosct, int dkoscl)
{ //поиск изменения степени черноты или безразмерного коэффициента поглощения
	double hf=1e0, dkusce=hf, dkoscee=hf, dko=hf; 
	int k=0, v=0;
	if (izm) { for (k=0; k<kost; k++) {
	dkusce=opredKTPTKToch(kusc, tkusc, temp[k], dmkoosc); 
	if (dkusce>hf) dkusce=hf; if (dkusce<0.0) dkusce=0.0;
	dkoscee=opredKTPTKToch(dkoscm, dkosct, temp[k], dkoscl); 
	if (dkoscee>hf) dkoscee=hf; if (dkoscee<0.0) dkoscee=0.0; //izm = 0 - нет изменений, izm=1 - учитываются изменения
	dko=dkusce*dkoscee;
	Aa[k]=Aa[k]*dko; Ta[k]=hf-Aa[k]-Ra[k]; Tb[k]=Ta[k]; Ab[k]=Aa[k]; } } //Ra - КО с учетом уменьшения, Rb - с учетом переотражений для внешнего излучения, Rc - для собственного излучения
		for (k=0; k<kost; k++) {
			Rb[k]=pow((hf-Ra[k])*Ta[k], 2e0)*Ra[k]/(hf-pow(Ra[k]*Ta[k], 2e0)*(hf-Ra[k]))+Ra[k]; 
			Rc[k]=pow((hf-Ra[k])*Tc[k], 2e0)*Ra[k]/(hf-pow(Ra[k]*Tc[k], 2e0)*(hf-Ra[k]))+Ra[k]; 
			Tb[k]=pow(hf-Ra[k], 2e0)*Ta[k]/(hf-pow(Ra[k]*Ta[k], 2e0)*(hf-Ra[k]));
			Tc[k]=pow(hf-Ra[k], 2e0)*Tc[k]/(hf-pow(Ra[k]*Tc[k], 2e0)*(hf-Ra[k]));
			Ab[k]=hf-Tb[k]-Rb[k]; 
			Ab[k]=hf-Tc[k]-Rc[k]; } 
	k=0; mu[k]=Ra; k++; mu[k]=Ta; k++; mu[k]=Aa; k++; mu[k]=Rb; k++; mu[k]=Tb; k++; mu[k]=Ab;
	k++; mu[k]=Rc; k++; mu[k]=Tc; k++; mu[k]=Ac;
	return mu;
}
//-----
double BolTochRasAlpha(int kost, double hvp, double por, double tc, double *ktpvo, double *vte, double *ete, int cem, int dmkvoz, 
	char *snm, double *ktptk, int vybves, int vybmar, int dmkoosc, double *kusc, double *tkusc, double *dkoscm, double *dkosct, 
	int dkoscl, double thol, double tgor, double y0, double dkoal, double dkosp, double hk, double qobsh, int kolPok, double htk, 
	double dkusce, double dkoscee, double dpctpe, double *stch, double p1, double p2, int vp, double *Refm, double *rts, int dmrts, 
	int vfqi, int vogf)
{
	int j=0, nk=j, l=0, k=0, zf=0, iz=1, v=0, c=0, w=0; //w - расчет температуры, iz - изменение
	double ***protv=NULL, ***prots=NULL, hf=1e0, hlr0=hf, hrl0=0.0, vyzn=0.0, ra=fabs(hlr0-hrl0), tol=0.0, *ktr=NULL, *reiz=NULL;
	double *Ts=NULL, *Tss=NULL, **mu=NULL, **mauk=NULL, *po=NULL, *sislr=NULL, *sisrl=NULL;
	double d=0.0, ka=d, kb=d, e=1e-9, mo=-hf, tx=tgor-fabs(tgor-thol)*hk/y0, r=0.0, ko=1e-6, hmi=1e-2;
	if (hvp>hmi) hvp=hvp*ko; if (htk>hmi) htk=htk*ko;
	for (j=0; j<kost; j++) d=d+hf; tol = d*htk + (d - hf)*hvp; 
	if ((vybves>=4) || (vybves<0)) { cout << "Net takogo veschestva!" << endl; k = getchar(); exit(1); }
	mauk=RaschRTA(kost, htk, ka, kb, iz, hvp, v, tc, c, w, vybves, cem, ete, hk, vybmar, dkusce, dkoscee, thol, tgor, y0, snm, dkosp, 
		kusc, tkusc, dmkoosc, dkoscm, dkosct, dkoscl, dkoal, Refm, rts, dmrts, p1*hvp, p2*hvp);
	double *Refa=NULL, *Traa=NULL, *Aba=NULL, *Ref=NULL, *Tra=NULL, *Ab=NULL, *Reft=NULL, *Trat=NULL, *Abt=NULL;
	k=0; Refa=mauk[k]; k++; Traa=mauk[k]; k++; Aba=mauk[k];  
	k++; Ref=mauk[k];  k++; Tra=mauk[k];  k++; Ab=mauk[k];   
	k++; Reft=mauk[k]; k++; Trat=mauk[k]; k++; Abt=mauk[k];  k++; nk=k; 
	j=0; protv=zadrkt(kost, snm, kolPok, Ref, Tra, htk, hvp, p1*hvp, p2*hvp, j); //vyv - выбор вещества: 0 - шамот, 1 - вермикулит
	j++; prots=zadrkt(kost, snm, kolPok, Ref, Tra, htk, hvp, p1*hvp, p2*hvp, j); 
	mu=opredTempLPSten(tc, 2*kost, ktpvo, vte, ete, ktptk, cem, dmkvoz, htk, hvp, qobsh, kost, snm, mo);
	k=0; Ts=mu[k]; k++; Tss=mu[k]; k++; if (mu) delete[]mu;
	ktr=KoefPoglRosselNac(ete, cem, dkosp, dkoal, vybves, vybmar); 
	hlr0=opredKTPTKToch(ktr, ete, tc, cem); if (hlr0<0.0) hlr0=0.0; 
	hlr0=fabs(hlr0*opredTempStenFragm(tc, 2*kost, ktpvo, vte, ete, ktptk, cem, dmkvoz, htk, hvp, qobsh, mo)/tol); 
	k=0; r=0.0; 
	po=RasIzlSerStenNacSok(Ref, Tra, Ab, protv, prots, tc, hlr0, hrl0, por, dpctpe, htk, hvp, kost, ete, ktptk, cem, ktpvo, vte, snm, stch, dmkvoz, 
		qobsh, hk, p1, p2, tx, dkosp, dkoal, vybves, vybmar, vp, ktr, cem, r, k, -1, vfqi, kolPok, vogf); 
	k=0; vyzn=po[k]; if (po) delete[]po;
	l=4; for (k=0; k<l; k++) { 
		mu=protv[k]; 
		for (j=0; j<kost; j++) { 
			po=mu[j]; if (po) delete[]po; } 
		mu=prots[k]; 
		for (j=0; j<kost; j++) { 
			po=mu[j]; if (po) delete[]po; } } 
	if (protv) delete[]protv; if (prots) delete[]prots; if (Ts) delete[]Ts; if (Tss) delete[]Tss; 
	for (k=0; k<nk; k++) { po=mauk[k]; if (po) delete[]po; } if (mauk) delete[]mauk; if (ktr) delete[]ktr;
	return vyzn;
}
//-----
double RasFracXeff(int fr, int v, int vyve, int vyma, char *snm) //fr - фракция, v - концентрация
{
	int l=2, k=0, j=0, i=0;
	if (vyve==1) { 
		if ((!fr) || (fr==3)) l=6; 
		if ((fr==1) || (fr==2)) l=9; } //если выбран  вермикулит
	double *mkbr=new double[l], *mv=new double[l], *tol=new double[l];
	double rokbr=2.75, *xv=NULL, xsr=0.0, hf=1e0, r=0.0, s=r, rov=hf;
	if (!vyve) { //если выбран шамот
	k=0; mkbr[k]=238.57; k++; mkbr[k]=227.973; k++;
	k=0; mv[k]=1.706;    k++; mv[k]=1.1;       k++;
	k=0; tol[k]=0.7;     k++; tol[k]=0.68;     k++;
	rov = 2.69; xv=EffectTols(mkbr, mv, tol, rov, rokbr, l, snm);
	k=0; j=1; if (!v) xsr = xv[k]; else if (v == 1) xsr = xv[j]; }
	//-----
	if (vyve==1) { //если выбран  вермикулит
	if (!fr) {
	k=0; mkbr[k]=250.239; k++; mkbr[k]=249.740; k++; mkbr[k]=249.223; k++;
		 mkbr[k]=250.395; k++; mkbr[k]=250.336; k++; mkbr[k]=249.55;
	k=0; mv[k]=0.283;     k++; mv[k]=0.464;     k++; mv[k]=0.812; k++;
		 mv[k]=0.22;      k++; mv[k]=0.547;     k++; mv[k]=0.777;
	k=0; tol[k]=0.73;     k++; tol[k]=0.72;     k++; tol[k]=0.72; k++;
		 tol[k]=0.72;     k++; tol[k]=0.71;     k++; tol[k]=0.7;
	rov=0.49; xv=EffectTols(mkbr, mv, tol, rov, rokbr, l, snm);
		if (!v) { k=0; j=3; } //выбор концентрации вермикулита
		if (v == 1) { k=1; j=4; }
		if (v == 2) { k=2; j=5; } 
		xsr=(xv[k]+xv[j])/2e0; }
	if (fr==1) {	
	k=0; mkbr[k]=250.0;  k++; mkbr[k]=249.629; k++; mkbr[k]=249.294; k++; mkbr[k]=249.706; k++;
	     mkbr[k]=249.51; k++; mkbr[k]=249.307; k++; mkbr[k]=250.328; k++; mkbr[k]=249.604; k++; mkbr[k]=249.206;
	k=0; mv[k]=0.255;    k++; mv[k]=0.539;     k++; mv[k]=0.809;     k++; mv[k]=0.295;     k++;
	     mv[k]=0.517;    k++; mv[k]=0.756;     k++; mv[k]=0.36;      k++; mv[k]=0.534;	   k++; mv[k]=0.843;
	k=0; tol[k]=0.72;    k++; tol[k]=0.71;     k++; tol[k]=0.7;      k++; tol[k]=0.7;      k++;
		 tol[k]=0.73;    k++; tol[k]=0.72;     k++; tol[k]=0.74;     k++; tol[k]=0.7;      k++; tol[k]=0.76;
	rov=0.52; xv=EffectTols(mkbr, mv, tol, rov, rokbr, l, snm);
	if (!v) { k=0; j=3; i=6; } 
	if (v == 1) { k=1; j=4; i=7; } 
	if (v == 2) { k=2; j=5; i=8; } 
	xsr=(xv[k]+xv[j]+xv[i])/3e0; }
	if (fr==2) {
	k=0; mkbr[k]=249.913; k++; mkbr[k]=249.607; k++; mkbr[k]=249.218; k++; mkbr[k]=249.929; k++;
		 mkbr[k]=249.695; k++; mkbr[k]=249.306; k++; mkbr[k]=250.405; k++; mkbr[k]=249.625; k++; mkbr[k]=249.348;
	k=0; mv[k]=0.315;     k++; mv[k]=0.473;     k++; mv[k]=0.709;     k++; mv[k]=0.293;     k++;
		 mv[k]=0.528;     k++; mv[k]=0.83;      k++; mv[k]=0.27;      k++; mv[k]=0.493;     k++; mv[k]=0.764;
	k=0; tol[k]=0.74;     k++; tol[k]=0.74;     k++; tol[k]=0.72;     k++; tol[k]=0.72;     k++;
		 tol[k]=0.71;     k++; tol[k]=0.7;      k++; tol[k]=0.78;     k++; tol[k]=0.73;     k++; tol[k]=0.76;
	rov = 0.53; xv = EffectTols(mkbr, mv, tol, rov, rokbr, l, snm);
	if (!v) { i=0; j=3; k=6; } 
	if (v == 1) { i=1; j=4; k=7; } 
	if (v == 2) { i=2; j=5; k=8; } 
	xsr=(xv[i]+xv[j]+xv[k])/3e0; }
	if (fr==3) {
	k=0; mkbr[k]=250.882; k++; mkbr[k]=249.590; k++; mkbr[k]=249.213; k++;
		 mkbr[k]=250.299; k++; mkbr[k]=249.441; k++; mkbr[k]=249.365;
	k=0; mv[k]=0.320;     k++; mv[k]=0.533;     k++; mv[k]=0.849; k++;
		 mv[k]=0.223;     k++; mv[k]=0.502;     k++; mv[k]=0.797;
	k=0; tol[k]=0.76;     k++; tol[k]=0.72;     k++; tol[k]=0.69; k++;
		 tol[k]=0.73;     k++; tol[k]=0.73;     k++; tol[k]=0.73;
	rov = 0.56; xv=EffectTols(mkbr, mv, tol, rov, rokbr, l, snm);
	if (!v) { j=0; k=3; } 
	if (v == 1) { j=1; k=4; } 
	if (v == 2) { j=2; k=5; } 
	xsr=(xv[j]+xv[k])/2e0; } }
	//-----
	if (vyve==2) { //если выбран ИТОМ
	k=0; mkbr[k]=230.078;
	k++; mkbr[k]=231.006;
	k=0; mv[k]=0.95;
	k++; mv[k]=1.19;
	k=0; tol[k]=0.64;
	k++; tol[k]=0.64;
	if (!vyma) rov=0.44; 
	if (vyma==1) rov=0.62; 
	if (vyma==2) rov=0.86;
	if (vyma==3) rov=hf;
	xv = EffectTols(mkbr, mv, tol, rov, rokbr, l, snm);
	if ((!v) || (v==1)) xsr = xv[v]; }
	if (vyve==3) { //если выбран КВИ
	r=230.0; k=0; mkbr[k]=r; k++; mkbr[k]=r; 
	r=0.23; k=0; mv[k]=r; k++; mv[k]=2e0*r;  
	r=0.64; k=0; tol[k]=r; k++; tol[k]=r; 
	s=1e1; r=0.0; for (k=0; k<vyma; k++) r=r+hf; r=r/s; rov=r;
	xv = EffectTols(mkbr, mv, tol, rov, rokbr, l, snm); xsr = 0.0;
	if ((!v) || (v==1)) xsr=xv[v]; }
	if (mkbr) delete[]mkbr; if (mv) delete[]mv; if (tol) delete[]tol; if (xv) delete[]xv; 
	return xsr;
}
double *EffectTols(double *mkbr, double *mv, double *tol, double rov, double rokbr, int n, char *snm)
{
	int k=0; 
	double *vkbr=new double[n], *vv=new double[n], *xvo=new double[n], hf=1e0, t=hf*1e3;
	if ((!vkbr) || (!vv) || (!xvo)) { cout << snm << endl; k=getchar(); exit(1); }
	for (k = 0; k < n; k++) {
		vkbr[k] = mkbr[k] / (t*rokbr);
		vv[k] = mv[k] / (t*rov);
		xvo[k] = (vv[k] / (vv[k] + vkbr[k]))*tol[k] * t; }
	if (vkbr) delete[]vkbr; if (vv) delete[]vv; 
	return xvo;
}
double RaschAlphaTvKar(int vybves, int vybmar, int vystsha, int vfv, int vsv, int isrp, int vpkf, double por, char *snm, int d, double dpctpe, 
	int vrsh, double *ktpvo, double *vte, double *ete, int cem, int dmkvoz, int kost, double *ktptk, int dmkoosc, double *kusc, double *tkusc, 
	double *dkoscm, double *dkosct, int dkoscl, double thol, double tgor, double y0, double dkoal, double dkosp, double hk, double qobsh, 
	int kolPok, double dkusce, double dkoscee, double tc, double *stch, double p1, double p2, int vp, double *Rmp, double *rts, int dmRef, 
	int vfqi, int vogf)
{
	int l=3, k=0, fr=0, nfr=4, q=0, cfk=l*nfr, cei=0, u=0, cm=5; 
	if (vybves!=1) { nfr=1; l=2; cfk=l*nfr; }
	long j=0, jk=100*100*100, pj=0, h=0, w=0;
	double x=0.0, y=x, p=x, xt=x, yp=x, hf=1e0, xa=x, xb=x, pore=x;
	double *srra=NULL, *porm=new double[nfr], srp=0.0, e=1e-8, koefc=1e-6, s=srp, r=s;
	double gfv1=6e1, gfv2=1e2, gfv3=15e1, gfv4=2e2, x1=0.0, x2=x1, x3=x1, x4=x1, gv0=63.0;
	double *pn=new double[d], *alx=new double[d], *rcf=new double[cfk];
	double *rapo=NULL, **unau=NULL, **mu=NULL, ce=0.0, *pp=NULL, lamtem=0.0; 
	double dvpn=1448.0/tc, df=0.0, tmss=df, htc=df, asp=0.0; //меньше 1 %
	double pr=0.0, rmf=0.0, prf=0.0, po=0.0, *alsf=new double[cfk]; //cfk - число различных вариантов таблеток
	double *prgr=NULL, *legr=NULL, *mez=NULL;
	if (!vybves) unau=vybFunRasPorPoRazSha(por, vybmar, vystsha, snm, vybves); 
	else if (vybves==1) unau=RasPorPoRazVer(por, vfv, vsv, isrp, vpkf, snm, vybves); 
	else if (vybves==2) unau=RasPorPoRazitom(vybmar, snm); 
	else if (vybves==3) unau=RasPorPoRazkvi(vybmar, snm); 
	else { cout << "Nepravilno vybran nomer!" << endl; k=getchar(); exit(1); }
	k=0; rapo=unau[k]; k++; srra=unau[k]; k++; prgr=unau[k]; 
	k++; legr=unau[k]; k++; mez=unau[k];
	k=0; srp=mez[k]; k++; ce=mez[k]; 
	unau=promPreobMas(unau);
	k=0; rapo=unau[k]; k++; srra=unau[k]; k++; prgr=unau[k]; 
	k++; legr=unau[k]; k++; mez=unau[k];
	k=0; srp=mez[k]; k++; ce=mez[k]; 
	for (j=0; j<d; j++) df=df+hf;
	xa=p1*srp; xb=p2*srp; x=srp; y=ce; 
	htc=opredTolTverChasObFigObPor(dpctpe, df, x, por, y, vogf, xa, xb);
	tmss=srp*(df-hf)+htc*df;
	j=0; s=x; x=e; while (x<ce) { x=x+hf; j++; } cei=j; //cout << "srp = " << s << "\tmrp = " << y << "\tcei = " << cei << "\thtc = " << htc << "\ttmss = " << tmss << "\t";
	if ((!pn) || (!alx) || (!rcf) || (!alsf) || (!porm)) { cout << snm << endl; j = getchar(); exit(1); }
		x1=0.0; yp=x1; x2=x1; x3=x1; x4=x1; 
		for (j=0; j<cei; j++) {
			x=srra[j]; p=rapo[j]; 
			if (vybves==1) { 
			if (yp<=gfv1) x1=x1+x*p; //размер пор должен быть ограничен размерами частиц
			if (yp<=gfv2) x2=x2+x*p; 
			if (yp<=gfv3) x3=x3+x*p; 
			if (yp<=gfv4) x4=x4+x*p; 
			yp=yp+hf; }
			else if (yp<gv0) { x1=x1+x*p; //cout << "x_i = " << x << "\tri = " << p << "\t"; 
			} }
		if (vybves==1) {
		k=0; porm[k]=x1*por/srp; k++; porm[k]=x2*por/srp; 
		k++; porm[k]=x3*por/srp; k++; porm[k]=x4*por/srp; } //пористость для фракций
		else { k=0; porm[k]=x1*por/srp; pore=porm[k]; //cout << "pore = " << pore << "\tx1 = " << x1 << "\tsrp = " << srp << "\tpor = " << por << "\t"; 
		}
	rmf=0.0; po=0.0; 
	for (j=0; j<RAND_MAX; j++) rmf=rmf+hf; 
	for (j=0; j<jk; j++) po=po+hf; 
	srand(time(0));
	for (j=0; j<d; j++) { pn[j]=0.0; alx[j]=0.0; } 
	for (j=0; j<cfk; j++) { rcf[j]=0.0; alsf[j]=0.0; } 
	u=0; 
	for (fr=0; fr<nfr; fr++) { //по фракциям
		for (k=0; k<l; k++) { //концентрация от 0,1 до 0,3 %
			x = RasFracXeff(fr, k, vybves, vybmar, snm); //эффективная толщина образца в таблетке
			rcf[u]=x; 
			p=porm[fr]; //cout << "por = " << p << "\tfr = " << fr << endl; //cout << "x = " << x << "\tk = " << k << "\tvv = " << vybves << "\tvm = " << vybmar << "\t"; 
			y=x*p; //размер поры 
			j=0; while (j<jk) {
				pj=rand();
				prf=0.0; for (h=0; h<pj; h++) prf=prf+hf; 
				pr=prf/rmf; //формируем случайное число от 0 до 1 с равномерной функцией распределения
				yp=y*pr; //выстрел по размеру пор
				p=porm[fr];
				pr=yp*koefc;
				xa=p1*pr;
				xb=p2*pr;
				if (fabs(p)>0.0) {
					xt=opredTolTverChasObFigObPor(dpctpe, df, pr, p, ce, vogf, xa, xb); if (xt>hf) xt=xt/koefc; }
				else continue; //определяем толщину стенки твердого каркаса
				xt=xt+yp; //пора с перегородкой
				if (fabs(xt)>0.0) p=x/xt; else continue; //определяем число стенок, укладывающихся в одной частице
				h=0; pr=hf+e; while ((pr<p) && (h<d)) { pr=pr+hf; h++; } 
				if (h<d) pn[h]=pn[h]+hf; 
				j++; } 
pr=0.0; for (j=0; j<d; j++) { pn[j]=pn[j]/po; pr=pr+pn[j]; } 
for (j = 0; j < d; j++) pn[j]=pn[j]/pr; //нормировка на единицу, получение вероятности для числа стенок, укладывающихся в тонком слое образца //pr=0.0; for (j=0; j<d; j++) pr=pr+pn[j]; //for (j=0; j<d; j++) cout << "pn ( " << j << " ) = " << pn[j] << "\t"; //cout << "sum = " << pr << endl; 
	for (j=0; j<d; j++) alx[j]=0.0;
	for (j=2; j<d; j++) {
		p=0.0; for (h=0; h<j; h++) p=p+hf;
		pore=porm[fr];
			if (p>hf) yp=pore*x/(p-hf); else yp=0.0;
			xt=opredTolTverChasObFigObPor(dpctpe, p, yp, pore, ce, vogf, xa, xb)/koefc; //cout << "xt = " << xt << "\typ = " << yp << "\t";
				if (yp>dvpn) { 
			lamtem=bbfn(yp, tc);
			if (lamtem<=0.0) lamtem=0.0; if (lamtem>hf) lamtem=hf; //cout << "lamtem = " << lamtem << "\t" << F0_lamT(tc, yp, snm) << "\t";
			alx[j]=BolTochRasAlpha(j, yp, pore, tc, ktpvo, vte, ete, cem, dmkvoz, snm, ktptk, vybves, vybmar, dmkoosc, 
				kusc, tkusc, dkoscm, dkosct, dkoscl, thol, tgor, y0, dkoal, dkosp, hk, qobsh, kolPok, xt, dkusce, dkoscee, 
				dpctpe, stch, p1, p2, vp, Rmp, rts, dmRef, vfqi, vogf); alx[j]=alx[j]*lamtem; } } //kolPok - число поколений
	for (j=0; j<d; j++) alsf[u]=alsf[u]+alx[j]*pn[j]; alsf[u]=alsf[u]*x;
	u++; } } //cout << "u = " << u << "\t";
	x=0.0; yp=0.0; for (j = 0; j < cfk; j++) { x = x + alsf[j]; yp = yp + rcf[j]; } //for (j=0; j<cfk; j++) cout << "alsf ( " << j << " ) = " << alsf[j] << "\trcf = " << rcf[j] << "\t"; cout << endl;
	asp=koefPoglSred(tc, vybves, vybmar, dkoal*dkosp); cout << "asp = " << asp << "\t";
	if ((fabs(yp)>0.0) && (fabs(asp)>0.0)) { x=x/yp; x=(x+asp)/asp; } else x=0.0; 
	for (k=0; k<cm; k++) { pp=unau[k]; if (pp) delete[]pp; } if (unau) delete[]unau; cout << "x = " << x << "\t";
	if (rcf) delete[]rcf; if (pn) delete[]pn; if (alx) delete[]alx; 
	if (alsf) delete[]alsf; if (porm) delete[]porm; printf("sr_okp = %0.10lf\talpha_sr = %0.10lf\n", x, asp); //ослабление КП за счет пористой структуры вермикулита //k=getchar();
	return x;
}
//-----
double UtochKoefPogl(double alp, double htc, double T, int vybves, int vybmar, double dko)
{
	double ns=opredPokPreSr(T, vybves, vybmar, dko), hf=1e0;
	double Refl=pow(ns-hf,2e0)/pow(ns+hf,2e0), Dnom=exp(-alp*htc), al=0.0;
	double Ap=(hf-Refl)*(hf-Dnom*(hf+Refl*Dnom-Refl))/(hf-pow(Refl*Dnom,2e0));
	al=fabs(log(Ap)/htc);
	return al;
}
double integroexponen(int n, double x)
{ 
	int Ni=100, k=0, m=0, j=0, i=0;
	double E=0.0, hf=1e0, s=0.0, e=0.577215665057043, nf=s, t=s, r=s, d=s, eps=1e-6, psi=s, ksi=s;
	if (((n==1) && (x<eps)) || (x<0.0)) { cout << "Beskonechnost!\tA = " << x << "\t"; E=hf; }
	else {
	for (k=0; k<n; k++) nf=nf+hf;
if ((x<eps) && (n>1))
    E=hf/(nf-hf);
if (x>eps) {
	s=0.0;
for (m=0; m<=Ni; m++) { 
    if (m==(n-1))
        continue; 
        t=hf;
        for (k=1; k<=m; k++) {
			r=0.0; for (j=0; j<k; j++) r=r+hf;
            t=t*r; }
		k=m-(n-1); d=hf; if (k<0) { d=-d; k=-k; }
		r=0.0; for (j=0; j<k; j++) r=r+d;
		ksi=0.0; for (j=0; j<m; j++) ksi=ksi+hf;
		if (fabs(r)>eps)
			s=s+pow(-x, ksi)/(r*t); } 
t=hf;
for (k=1; k<=(n-1); k++) {
	r=0.0; for (j=0; j<k; j++) r=r+hf;
	t=t*r; } //(n-1)!
k=n-1; r=0.0; for (j=0; j<k; j++) r=r+hf;
E=pow(-x, r); //(-x)^(n-1)
E=E/t; 
psi=-e;
if (n>1) { 
t=0.0;
for (k=1; k<=(n-1); k++) {
	r=0.0; for (j=0; j<k; j++) r=r+hf;
	t=t+hf/r; }
psi=psi+t; } 
if (x>eps) E=E*(-log(x)+psi); else E=0.0;
E=E-s; } }
return E; }
double **promPreobMas(double **mu)
{
	double hmi=1e-3, hf=1e0, s=0.0, w=0.0, koefc=1e6, srp=0.0, mrp=0.0, r=1e1, ko=1e-2;
	int k=0, qg=0;
	double *prgr01=NULL, *srra=NULL, *legr01=NULL, *ras02=NULL, *mez=NULL, x=0.0, e=1e-1;
	k=0; ras02=mu[k];  k++; srra=mu[k]; k++; prgr01=mu[k]; 
	k++; legr01=mu[k]; k++; mez=mu[k];
	k=0; srp=mez[k]; k++; mrp=mez[k];
	if (mrp>hmi) koefc=hf; 
	w=koefc; s=w*mrp;
	k=0; x=e; while (x<s) { x=x+hf; k++; } qg=k; //cout << "ce = " << qg << "\t"; //for (k=0; k<w; k++) cout << "rpr ( " << k << " ) = " << rapo[k] << "\tsrp = " << srra[k] << "\t"; k=getchar();
for (k=0; k<qg; k++) { srra[k]=srra[k]*w; legr01[k]=legr01[k]*w; prgr01[k]=prgr01[k]*w; }
srp=srp*w; mrp=mrp*w;
s=0.0; for (k=0; k<qg; k++) s=s+ras02[k];
if (s>r) w=ko; else w=hf; for (k=0; k<qg; k++) ras02[k]=ras02[k]*w;
k=0; mez[k]=srp; k++; mez[k]=mrp;
k=0; mu[k]=ras02;  k++; mu[k]=srra; k++; mu[k]=prgr01; 
k++; mu[k]=legr01; k++; mu[k]=mez;
return mu;
}