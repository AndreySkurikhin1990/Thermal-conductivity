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
double **KoefPoglSten(double, double, int, double, int, double, double, int, double, double, double, int, double, double, double, int, int);
double **izmRTA(int, int, double *, double *, double *, double *, double *, double *, double *, double *, double *, char *, double **, double *, 
	double *, double *, int, double *, double *, int);
char **napolStrok(int, int, char);
char **napNazvFile(int, int, char *, char);
double F0_lamT(double, char *);
double *EffectTols(double *, double *, double *, double, double, int, char *);
double RaschAlphaTvKar(int, int, int, int, int, int, int, double, char *, int, double, int, double *, double *, double *, int, int, 
	int, double *, int, double *, double *, double *, double *, int, double, double, double, double, double, double, double, int, 
	double, double, double, double *, double, double, int, double *, double *, int, int);
double *KoefPoglRosselNac(double *, int, double, double, int, int);
double **RasPorPoRazVer(double, int, int, int, int, char *);
double **opredTempLPSten(double, int, double *, double *, double *, double *, int, int, double, double, double, int, char *, double);
double UtochKoefPogl(double, double, double, int, int, double);
double **RaschSobLuchPlotTepPot(int, double *, double *, double *, double *, double *, double *, double, double, double *, int, char *); //double **VybFunRasPorpoRazSha(double, int, int); double **RasPorPoRazitom(int); double **RasPorpoRazkvi(int);
double opredPokPreSr(double, int, int, double);
double **oprRasTemNach(int, double, double, double, double, double, int, double, double, double, char *);
double bbfn(double);
double *opredKoefOtr(double *, int, int, int, double *, int, char *, double);
double opredKTPTKToch(double *, double *, double, int);
double ***zadrkt(int, char *, int, double *, double *, double, double, double, double, int);
double BolTochRasAlpha(int, double, double, double, double *, double *, double *, int, int, char *, double *, int, int, int, double *, 
	double *, double *, double *, int, double, double, double, double, double, double, double, int, double, double, double, double, 
	double *, double, double, int, double *, double *, int, int);
double oprSrRazPor(int, int, double, int, int, int, int, int, int, char *);
double oprMaxRazPor(int, int, double, int, int, int, int, int, int, char *);
double opredTolTverChasObFig(double, double, double, double, double);
double opredTolTverChasObFigObPor(double, double, double, double, double);
double opredTempStenFragm(double, int, double *, double *, double *, double *, int, int, double, double, double, double);
double *RasIzlSerStenNacSok(double *, double *, double *, double ***, double ***, double, double, double, double, double, double, double, int, double *, 
	double *, int, double *, double *, char *, double *, int, double, double, double, double, double, double, double, int, int, int,
	double *, int, double, int, int, int, int);
double koefPoglSred(double, int, int, double);
double integroexponen(int, double);
double VIEW(int, int, double *, int);
//-------
double **RaschRTA(int kost, double htk, double kta, double ktb, int izm, double hvo, int v, double ti, int c, int w, int vybves, int cem, 
	double *ete, double hk, int vybmar, double dkusce, double dkoscee, double thol, double tgor, double y0, char *snm, double dkosp, 
	double *kusc, double *tkusc, int dmkoosc, double *dkoscm, double *dkosct, int dkoscl, double dkoal, double *Rmp, double *rts, int dmRef,
	double a, double b)
{ 
	int k=kost, j=0, q=0, x=0; //cout << "ht = " << htk << "\thv = " << hvo << endl;
	double *te=NULL, **mauk=NULL, hf=1e0, *alpha=NULL, t=0.0, r=0.0, *po=NULL, dko=dkoal*dkosp, *arg=NULL; 
	if (dko<0.0) dko=0.0; if (dko>hf) dko=hf;
	if (!v) { if ((vybves>=0) && (vybves<4)) { //определение коэффициентов поглощения при температуре стенок
		mauk=KoefPoglSten(kta, ktb, v, ti, kost, htk, hvo, x, dkoal, dkosp, hk, cem, thol, tgor, y0, vybves, vybmar);
		q=0; alpha=mauk[q]; q++; te=mauk[q]; /*for (k=0; k<kost; k++) cout << "al = " << alpha[k] << "\tte = " << te[k] << endl;*/ }
		if (!c) { k=kost; //задание коэффициентов отражения, поглощения и пропускания: v=0, c=0
double *Ra=NULL, *R=new double[kost], *Rt=new double[kost]; //Ra=opredKoefOtr(te, k, vybves, cem, ete, vybmar, snm, dko); //for (k=0; k<kost; k++) cout << "R = " << Ra[k] << "\t"; cout << endl;
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
	A[k]=Aa[k];
	q=3; Tt[k]=integroexponen(q, alpha[k]*htk)*VIEW(j, x, arg, k);
	At[k]=hf-Rt[k]-Tt[k]; }
		if (mauk) delete[]mauk; if (arg) delete[]arg;
		k=9; mauk=new double*[k]; if (!mauk) { cout << snm << endl; k=getchar(); exit(1); }
		k=kost; q=1;
		mauk=izmRTA(k, q, R, T, A, Ra, Ta, Aa, Rt, Tt, At, snm, mauk, te, kusc, tkusc, dmkoosc, dkoscm, dkosct, dkoscl); 
		if (alpha) delete[]alpha; if (te) delete[]te; } } //отправляется массив температур
	else if (v == 1) mauk=KoefPoglSten(t, r, v, ti, kost, htk, hvo, w, dkoal, dkosp, hk, cem, thol, tgor, y0, vybves, vybmar); //поиск коэффициента поглощения при температуре ti
	return mauk;
}
//-----
double **KoefPoglSten(double kta, double ktb, int vyb, double tl, int kost, double htk, double hvp, int w, double dkoal, 
	double dkosp, double hk, int cem, double thol, double tgor, double tol, int vybves, int vybmar)
{
	int k=0, q=100, j=0, r=0, f=0, i=2, u=1;
	char ch='\0', **unau=napolStrok(vybves, vybmar, ch), *snm=unau[k], *sfno=unau[u];
	char **mauk=napNazvFile(vybves, vybmar, snm, ch), *skpt=mauk[i];
	double te0=273.15, t0=22.0, te=t0+te0, dte=1e0, t=0.0, *p=NULL, dko=0.0, hf=1e0;
	double a1=0.0, a2=0.0, t1=0.0, t2=0.0, **mu=NULL, **muv=NULL;
	char *s=new char[q]; KoePog *kp=new KoePog, *ne=NULL, *roo=NULL, *pre=NULL;
	if ((!kp) || (!s)) { cout << snm << endl; j=getchar(); exit(1); } for (j=0; j<q; j++) s[j]=ch;
	ifstream fin; fin.open(skpt); if (!fin.is_open()) { cout << sfno << endl; j=getchar(); exit(1); }
	if (dkosp>hf) dkosp=hf; if (dkosp<0.0) dkosp=0.0; if (dkoal>hf) dkoal=hf; if (dkoal<0.0) dkoal=0.0;
	dko=dkoal*dkosp;
	roo=kp; k=0; while (!fin.eof()) {
		fin.getline(s, q, '\n');
		ne = new KoePog; if (!ne) { cout << snm << endl; j = getchar(); exit(1); }
		kp->alp=atof(s)*dko; kp->tem=te; 
		kp->nex=ne; pre=kp; kp=ne; k++; te=te+dte; kp->nex=NULL; }
	if (ne) delete[]ne; kp=NULL; pre->nex=kp; if (s) delete[]s; r=k; fin.close(); 
	if (!vyb) { //поиск массива коэффициентов поглощения при температурах стенок
		f=2; mu=new double*[i]; 
		double knat=0.0, *al=new double[kost], *teks=NULL, *xi=NULL;
		if ((!al) || (!mu)) { cout << snm << endl; j=getchar(); exit(1); }
		muv=oprRasTemNach(kost, htk, hvp, hk, kta, ktb, w, thol, tgor, tol, snm); //cout << "w = " << w << "\tvyb = " << vyb << endl;
		k=0; teks=muv[k]; k++; xi=muv[k]; 
		for (k=0; k<kost; k++) {
			kp=roo->nex; pre=roo; j=0; t=teks[k];
			while ((kp) && (j<r)) {
				if (kp->tem>=t) {
					a1=pre->alp; a2=kp->alp; 
					t1=pre->tem; t2=kp->tem; 
					knat=(a2-a1)/(t2-t1); 
					al[k]=a1+knat*(t-t1); break; }
				else { pre=kp; kp=kp->nex; } j++; } } 
		if (xi) delete[]xi; if (muv) delete[]muv;
		k=0; mu[k]=al; k++; mu[k]=teks; }
	if (vyb==1) { //поиск конкретного коэффициента поглощения при данной температуре tl
		kp=roo->nex; j=0; pre=roo; t=0.0;
		while ((kp) && (j<r)) {
			if (kp->tem>tl) {
				a1=pre->alp; a2=kp->alp; 
				t1=pre->tem; t2=kp->tem; 
				t=a1+(a2-a1)*(tl-t1)/(t2-t1); break; }
			else { pre=kp; kp=kp->nex; } j++; } j=1; p=new double[j]; j=0; w=1; mu=new double*[w];
			if ((p) && (mu)) { p[j]=t; mu[j]=p; } else { cout << snm << endl; j=getchar(); exit(1); } }
	kp=roo; while (kp) { ne=kp->nex; delete kp; kp=ne; } //удаление списка
	j=4; for (k=0; k<j; k++) { skpt=mauk[k]; if (skpt) delete[]skpt; } if (mauk) delete[]mauk; //удаление имен файлов
	j=2; for (k=0; k<j; k++) { skpt=unau[k]; if (skpt) delete[]skpt; } if (unau) delete[]unau; 
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
	dkusce=opredKTPTKToch(kusc, tkusc, temp[k], dmkoosc); if (dkusce>hf) dkusce=hf; if (dkusce<0.0) dkusce=0.0;
	dkoscee=opredKTPTKToch(dkoscm, dkosct, temp[k], dkoscl); if (dkoscee>hf) dkoscee=hf; if (dkoscee<0.0) dkoscee=0.0; //izm = 0 - нет изменений, izm=1 - учитываются изменения
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
	int vfqi)
{
	int j=0, nk=j, l=0, k=0, zf=0, iz=1, v=0, c=0, w=0; //w - расчет температуры, iz - изменение
	double ***protv=NULL, ***prots=NULL, hf=1e0, hlr0=hf, hrl0=0.0, vyzn=0.0, ra=fabs(hlr0-hrl0), tol=0.0, *ktr=NULL, *reiz=NULL;
	double *Ts=NULL, *Tss=NULL, **mu=NULL, **mauk=NULL, *po=NULL, *sislr=NULL, *sisrl=NULL;
	double d=0.0, ka=d, kb=d, e=1e-9, mo=-hf, tx=tgor-fabs(tgor-thol)*hk/y0, r=0.0;
	for (j=0; j<kost; j++) d=d+hf; tol = d*htk + (d - hf)*hvp;
	if (vybves>=4) { cout << "Net takogo veschestva!" << endl; k = getchar(); exit(1); }
	mauk=RaschRTA(kost, htk, ka, kb, iz, hvp, v, tc, c, w, vybves, cem, ete, hk, vybmar, dkusce, dkoscee, thol, tgor, y0, snm, dkosp, 
		kusc, tkusc, dmkoosc, dkoscm, dkosct, dkoscl, dkoal, Refm, rts, dmrts, p1*hvp, p2*hvp);
	double *Refa=NULL, *Traa=NULL, *Aba=NULL, *Ref=NULL, *Tra=NULL, *Ab=NULL, *Reft=NULL, *Trat=NULL, *Abt=NULL;
	k=0; Refa=mauk[k]; k++; Traa=mauk[k]; k++; Aba=mauk[k];  k++; Ref=mauk[k];
	k++; Tra=mauk[k];  k++; Ab=mauk[k];   k++; Reft=mauk[k]; k++; Trat=mauk[k]; k++; Abt=mauk[k]; k++; nk=k; //for (k=0; k<kost; k++) cout << "R = " << Ref[k] << "\tT = " << Tra[k] << "\tA = " << Ab[k] << endl;
	j=0; protv=zadrkt(kost, snm, kolPok, Ref, Tra, htk, hvp, p1*hvp, p2*hvp, j); //vyv - выбор вещества: 0 - шамот, 1 - вермикулит
	j++; prots=zadrkt(kost, snm, kolPok, Ref, Tra, htk, hvp, p1*hvp, p2*hvp, j); 
	mu=opredTempLPSten(tc, 2*kost, ktpvo, vte, ete, ktptk, cem, dmkvoz, htk, hvp, qobsh, kost, snm, mo);
	k=0; Ts=mu[k]; k++; Tss=mu[k]; k++; if (mu) delete[]mu;
	ktr=KoefPoglRosselNac(ete, cem, dkosp, dkoal, vybves, vybmar); //for (k=0; k<cem; k++) cout << "ktr ( " << k << " ) = " << ktr[k] << "\t"; cout << endl;
	hlr0=opredKTPTKToch(ktr, ete, tc, cem); if (hlr0<0.0) hlr0=0.0; //cout << "KTP_Ros = " << hlr0 << endl;
	hlr0=fabs(hlr0*opredTempStenFragm(tc, 2*kost, ktpvo, vte, ete, ktptk, cem, dmkvoz, htk, hvp, qobsh, mo)/tol); //cout << "hlr = " << hlr0 << endl;
	k=0; r=0.0; 
	po=RasIzlSerStenNacSok(Ref, Tra, Ab, protv, prots, tc, hlr0, hrl0, por, dpctpe, htk, hvp, kost, ete, ktptk, cem, ktpvo, vte, snm, stch, dmkvoz, 
		qobsh, hk, p1, p2, tx, dkosp, dkoal, vybves, vybmar, vp, ktr, cem, r, k, -1, vfqi, kolPok); //for (k=0; k<kost; k++) cout << "Ts = " << Ts[k] << "\tTss = " << Tss[k] << "\t\t"; cout << endl; 
	k=0; vyzn=po[k]; if (po) delete[]po;
	if (protv) delete[]protv; if (prots) delete[]prots; if (Ts) delete[]Ts; if (Tss) delete[]Tss; 
	for (k=0; k<nk; k++) { po=mauk[k]; if (po) delete[]po; } if (mauk) delete[]mauk; if (ktr) delete[]ktr;
	return vyzn;
}
//-----
double RasFracXeff(int fr, int v, int vyve, int vyma, char *snm)
{
	int l=6, k=0, j=0, i=0; 
	if (!vyve) l=2; //если выбран шамот
	else if (vyve==1) { if ((!fr) || (fr==3)) l=6; else if ((fr==1) || (fr==2)) l=9; } //если выбран  вермикулит
	double *mkbr=new double[l], *mv=new double[l], *tol=new double[l];
	double rokbr=2.75, rov=0.0, *xv=NULL, xsr=0.0, hf=1e0;
	if (!vyve) { //если выбран шамот
	k=0; mkbr[k]=238.57; k++; mkbr[k]=227.973; k++;
	k=0; mv[k]=1.706;    k++; mv[k]=1.1;       k++;
	k=0; tol[k]=0.7;     k++; tol[k]=0.68;     k++;
	rov = 2.69; xv=EffectTols(mkbr, mv, tol, rov, rokbr, l, snm);
	k=0; j=1; if (!v) xsr = xv[k]; else if (v == 1) xsr = xv[j]; }
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
		if (v == 2) { k=2; j=5; } xsr=(xv[k]+xv[j])/2e0; }
	else if (fr==1) {	
	k=0; mkbr[k]=250.0;  k++; mkbr[k]=249.629; k++; mkbr[k]=249.294; k++; mkbr[k]=249.706; k++;
	     mkbr[k]=249.51; k++; mkbr[k]=249.307; k++; mkbr[k]=250.328; k++; mkbr[k]=249.604; k++; mkbr[k]=249.206;
	k=0; mv[k]=0.255;    k++; mv[k]=0.539;     k++; mv[k]=0.809;     k++; mv[k]=0.295;     k++;
	     mv[k]=0.517;    k++; mv[k]=0.756;     k++; mv[k]=0.36;      k++; mv[k]=0.534;	   k++; mv[k]=0.843;
	k=0; tol[k]=0.72;    k++; tol[k]=0.71;     k++; tol[k]=0.7;      k++; tol[k]=0.7;      k++;
		 tol[k]=0.73;    k++; tol[k]=0.72;     k++; tol[k]=0.74;     k++; tol[k]=0.7;      k++; tol[k]=0.76;
	rov=0.52; xv=EffectTols(mkbr, mv, tol, rov, rokbr, l, snm);
	if (!v) { k=0; j=3; i=6; } 
	if (v == 1) { k=1; j=4; i=7; } 
	if (v == 2) { k=2; j=5; i=8; } xsr=(xv[k]+xv[j]+xv[i])/3e0; }
	else if (fr==2) {
	k=0; mkbr[k]=249.913; k++; mkbr[k]=249.607; k++; mkbr[k]=249.218; k++; mkbr[k]=249.929; k++;
		 mkbr[k]=249.695; k++; mkbr[k]=249.306; k++; mkbr[k]=250.405; k++; mkbr[k]=249.625; k++; mkbr[k]=249.348;
	k=0; mv[k]=0.315;     k++; mv[k]=0.473;     k++; mv[k]=0.709;     k++; mv[k]=0.293;     k++;
		 mv[k]=0.528;     k++; mv[k]=0.83;      k++; mv[k]=0.27;      k++; mv[k]=0.493;     k++; mv[k]=0.764;
	k=0; tol[k]=0.74;     k++; tol[k]=0.74;     k++; tol[k]=0.72;     k++; tol[k]=0.72;     k++;
		 tol[k]=0.71;     k++; tol[k]=0.7;      k++; tol[k]=0.78;     k++; tol[k]=0.73;     k++; tol[k]=0.76;
	rov = 0.53; xv = EffectTols(mkbr, mv, tol, rov, rokbr, l, snm);
	if (!v) { i=0; j=3; k=6; } 
	if (v == 1) { i=1; j=4; k=7; } 
	if (v == 2) { i=2; j=5; k=8; } xsr=(xv[i]+xv[j]+xv[k])/3e0; }
	else if (fr==3) {
	k=0; mkbr[k]=250.882; k++; mkbr[k]=249.590; k++; mkbr[k]=249.213; k++;
		 mkbr[k]=250.299; k++; mkbr[k]=249.441; k++; mkbr[k]=249.365;
	k=0; mv[k]=0.320;     k++; mv[k]=0.533;     k++; mv[k]=0.849; k++;
		 mv[k]=0.223;     k++; mv[k]=0.502;     k++; mv[k]=0.797;
	k=0; tol[k]=0.76;     k++; tol[k]=0.72;     k++; tol[k]=0.69; k++;
		 tol[k]=0.73;     k++; tol[k]=0.73;     k++; tol[k]=0.73;
	rov = 0.56; xv=EffectTols(mkbr, mv, tol, rov, rokbr, l, snm);
	if (!v) { j=0; k=3; } 
	if (v == 1) { j=1; k=4; } 
	if (v == 2) { j=2; k=5; } xsr=(xv[j]+xv[k])/2e0; } }
	else if (vyve==2) { //если выбран ИТОМ
	k=0; mkbr[k]=230.078; k++; mkbr[k]=231.006; k++;
	k=0; mv[k]=0.95;      k++; mv[k]=1.19;      k++;
	k=0; tol[k]=0.64;     k++; tol[k]=0.64;     k++;
	if (!vyma) rov=0.44; else if (vyma==1) rov=0.62; 
	else if (vyma==2) rov=0.86; else if (vyma==4) rov=1e0;
	xv = EffectTols(mkbr, mv, tol, rov, rokbr, l, snm);
	if ((!v) || (v==1)) xsr = xv[v]; }
	else if (vyve==3) { //если выбран КВИ
	k=0; mkbr[k]=230.078; k++; mkbr[k]=231.006; k++;
	k=0; mv[k]=0.95;  k++; mv[k]=1.19;  k++;
	k=0; tol[k]=0.64; k++; tol[k]=0.64; k++;
	rov=0.4; xv = EffectTols(mkbr, mv, tol, rov, rokbr, l, snm); xsr = 0.0;
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
	int vfqi)
{
	int l=3, k=0, fr=0, nfr=4, q=0, cfk=l*nfr, cei=0, u=0;
	long j=0, jk=100*100*100, pj=0, h=0, w=0;
	double x=0.0, y=0.0, p=0.0, xt=0.0, yp=0.0, hf=1e0;
	double *srra=NULL, *porm=new double[nfr], srp=0.0, e=1e-10, koefc=1e-6, s=srp, r=s;
	double gfv1=6e1, gfv2=1e2, gfv3=15e1, gfv4=2e2, x1=0.0, x2=x1, x3=x1, x4=x1;
	double *pn=new double[d], *alx=new double[d], *rcf=new double[cfk];
	double *rapo=NULL, **unau=NULL, **mu=NULL, ce=0.0, *pp=NULL, lamtem=0.0; 
	double dvpn=1448.0/tc, df=0.0, tmss=df, htc=df, asp=0.0; //меньше 1 %
	double pr=0.0, rmf=0.0, prf=0.0, po=0.0, *alsf=new double[cfk]; //cfk - число различных вариантов таблеток
	if (vybves==1) unau=RasPorPoRazVer(por, vfv, vsv, isrp, vpkf, snm); //if (!vybves) unau=VybFunRasPorpoRazSha(por, vybmar, vystsha); if (vybves==2) unau=RasPorPoRazitom(vybmar); if (vybves==3) unau=RasPorpoRazkvi(vybmar); else { cout << "Nepravilno vybran nomer!" << endl; k=getchar(); exit(1); }
	k=0; rapo=unau[k]; k++; srra=unau[k]; 
	srp=oprSrRazPor(vybves, vybmar, por, vfv, vsv, isrp, vpkf, vrsh, vystsha, snm)/koefc; //в микронах
	ce=oprMaxRazPor(vybves, vybmar, por, vfv, vsv, isrp, vpkf, vrsh, vystsha, snm)/koefc; //в микронах
	for (j=0; j<d; j++) df=df+hf; 
	htc=opredTolTverChasObFig(dpctpe, df, srp*koefc, por, ce*koefc);
	if (htc<e) htc=opredTolTverChasObFigObPor(dpctpe, df, srp*koefc, por, ce*koefc);
	tmss=srp*(df-hf)+htc*df;
	j=0; x=e; while (x<ce) { x=x+hf; j++; } cei=j; cout << endl;
	if ((!pn) || (!alx) || (!rcf) || (!alsf) || (!porm)) { cout << snm << endl; j = getchar(); exit(1); }
		x1=0.0; yp=x1; x2=x1; x3=x1; x4=x1;
		cout << "cei = " << cei << "\tce = " << ce << "\tsrp = " << srp << "\ttmss = " << tmss << "\tpor = " << por << "\tdp = " << dpctpe << "\thtc = " << htc << "\tdv_min = " << dvpn << endl;
		for (j=0; j<cei; j++) {
			x=srra[j]; p=rapo[j];
			if (yp<=gfv1) x1=x1+x*p; //размер пор должен быть ограничен размерами частиц
			if (yp<=gfv2) x2=x2+x*p; 
			if (yp<=gfv3) x3=x3+x*p; 
			if (yp<=gfv4) x4=x4+x*p; 
			yp=yp+hf; }
		k=0; porm[k]=x1*por/srp; k++; porm[k]=x2*por/srp; 
		k++; porm[k]=x3*por/srp; k++; porm[k]=x4*por/srp; //пористость для фракций
	rmf=0.0; po=rmf; for (j=0; j<RAND_MAX; j++) rmf=rmf+hf; for (j=0; j<jk; j++) po=po+hf; //cout << "x1 = " << x1 << "\tx2 = " << x2 << "\tx3 = " << x3 << "\tx4 = " << x4; cout << endl; //for (j=0; j<nfr; j++) cout << "j = " << j << "\tporm = " << porm[j] << "\t"; cout << endl;
	srand(time(0));
	for (j=0; j<d; j++) { pn[j]=0.0; alx[j]=0.0; } for (j=0; j<cfk; j++) { rcf[j]=0.0; alsf[j]=0.0; } u=0; //long lt=0; unsigned int st=0; lt = time(NULL); cout << "t = " << lt << endl; st = (unsigned int)(lt - (lt % 2)) / 2; cout << "st = " << st << endl;
	for (fr=0; fr<nfr; fr++) {
		for (k=0; k<l; k++) { //концентрация от 0,1 до 0,3 %
			x = RasFracXeff(fr, k, vybves, vybmar, snm); //эффективная толщина образца в таблетке
			rcf[u]=x; 
			y=x*porm[fr]; //размер поры 
			cout << "Tol Obr = " << x << "\ty = " << y << "\tp = " << porm[fr] << "\tfr = " << fr << "\tk = " << k << "\tlTm = " << dvpn << endl;
			j=0; while (j<jk) {
				pj=rand();
				prf=0.0; for (h=0; h<pj; h++) prf=prf+hf; 
				pr=prf/rmf; //формируем случайное число от 0 до 1 с равномерной функцией распределения
				yp=y*pr; //выстрел по размеру пор
				p=porm[fr]; 
				if (fabs(p)>0.0) xt=opredTolTverChasObFig(dpctpe, df, yp*koefc, p, ce*koefc)/koefc; //xt=yp*(hf-p)/p; //определяем толщину стенки твердого каркаса
				xt=xt+yp;
				if (fabs(xt)>0.0) p=x/xt; else continue; //определяем число стенок, укладывающихся в одной частице
				h=0; pr=hf+e; while (pr<p) { pr=pr+hf; h++; } 
				if (h<d) pn[h]=pn[h]+hf; 
				j++; } 
pr=0.0; for (j=0; j<d; j++) { pn[j]=pn[j]/po; pr=pr+pn[j]; } for (j = 0; j < d; j++) pn[j] = pn[j]/pr; //нормировка на единицу, получение вероятности для числа стенок, укладывающихся в тонком слое образца //pr=0.0; for (j=0; j<d; j++) pr=pr+pn[j];
	for (j=0; j<d; j++) cout << "pn ( " << j << " ) = " << pn[j] << "\t"; //cout << "sum = " << pr << endl; 
	for (j=0; j<d; j++) alx[j]=0.0;
	for (j=2; j<d; j++) {
		p = 0.0; for (h = 0; h<j; h++) p = p + hf;
			yp = porm[fr]*x/(p-hf); xt=opredTolTverChasObFig(dpctpe, p, yp*koefc, porm[fr], ce*koefc)/koefc; //xt = (hf - porm[fr])*x/p;
				if (yp>dvpn) {
					cout << "Tol voz Pros = " << yp << "\ttc = " << tc << "\tlam*T = " << yp*tc << "\txt = " << xt << "\tks = " << j << endl;
			lamtem=bbfn(yp*tc); if (lamtem<=0.0) lamtem=0.0; if (lamtem>hf) lamtem=hf; cout << "F(lam*T) = " << lamtem << "\t"; //lamtem=F0_lamT(yp*tc*koefc, snm); if (lamtem<=0.0) lamtem=0.0; if (lamtem>hf) lamtem=hf; cout << "F(lam*T) = " << lamtem << "\n";
			alx[j]=BolTochRasAlpha(j, yp*koefc, porm[fr], tc, ktpvo, vte, ete, cem, dmkvoz, snm, ktptk, vybves, vybmar, dmkoosc, 
				kusc, tkusc, dkoscm, dkosct, dkoscl, thol, tgor, y0, dkoal, dkosp, hk, qobsh, kolPok, xt*koefc, dkusce, dkoscee, 
				dpctpe, stch, p1, p2, vp, Rmp, rts, dmRef, vfqi); alx[j]=alx[j]*lamtem; } } //kolPok - число поколений
	for (j=0; j<d; j++) alsf[u]=alsf[u]+alx[j]*pn[j]; alsf[u]=alsf[u]*x; u++; } } //for (j=0; j<d; j++) cout << "j = " << j << "\talx = " << alx[j] << "\t"; cout << endl; 
	for (j=0; j<cfk; j++) cout << "j = " << j << "\talsf = " << alsf[j] << "\tra_ch = " << rcf[j] << "\t"; cout << endl; 
	x = 0.0; yp = 0.0; for (j = 0; j < cfk; j++) { x = x + alsf[j]; yp = yp + rcf[j]; } 
	asp=koefPoglSred(tc, vybves, vybmar, dkoal*dkosp);
	if (fabs(yp)>e) { x = x / yp; if (fabs(asp)>e) x=(x+asp)/asp; } else x=0.0;
	q=6; for (k=0; k<q; k++) { pp=unau[k]; if (pp) delete[]pp; } if (unau) delete[]unau;
	if (rcf) delete[]rcf; if (pn) delete[]pn; if (alx) delete[]alx; 
	if (alsf) delete[]alsf; if (porm) delete[]porm;
	printf("sr_okp = %0.10lf\talpha_sr = %0.10lf\n", x, asp); //ослабление КП за счет пористой структуры вермикулита //k=getchar();
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
	if (((n==1) && (x<eps)) || (x<0.0)) { cout << "Beskonechnost!" << endl; k=getchar(); exit(1); }
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
E=E-s; }
return E; }