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
struct derevo {
	int otre; //1 - отражение или 0 - пропускание
	int ste; //номер стенки
	int vis; //видимость: 1 - виден, 0 - нет
	int lev; //номер уровня
	struct derevo *back; //указатель назад
	struct derevo *next; //указатель вперед
};
//-----
double *podschchiele(int, int, int, double *, double *, char *, int, double **);
double rasotprpovs(int *, int *, int, int, int, int, double *, double *, double **);
double otrprovs(int, int, double *, double *, int, double **, int *, int);
double **opredTempLPSten(double, int, double *, double *, double *, double *, int, int, double, double, double, int, char *, double);
double **RaschVneshIzluch(double *, double *, double ***, int, double, double, char *);
double *opredTempStenFragmMasTem(int, double *, double *, double *, double *, int, int, double, double, double, double, double, char *, int, int);
double opredKTPTKToch(double *, double *, double, int);
double *usrednen(double *, double *, int, int);
double *SeryeStenkiRasIzl(double, double, double, double *, double *, double *, double *, double *, int *, int, char *, int);
double *SeryeStenkiRasIzlDifPov(double, double, double, double *, double *, double *, double *, int *, int, char *);
double *RasIzlSerStenNacSok(double *, double *, double *, double ***, double ***, double, double, double, double, double, double, double, int, 
	double *, double *, int, double *, double *, char *, double *, int, double, double, double, double, double, double, double, int, int, int,
	double *, int, double, int, int, int, int);
double ***zadrkt(int, char *, int, double *, double *, double, double, double, double, int);
double *izstN(int, int, int, int, double *, double *, int, char *, double **);
double *opredTempTvKarFragm(int, double *, double *, int, double, double, double, double, double, char *, double, double *);
double *KoefPoglRosselNac(double *, int, double, double, int, int);
double opredTempStenFragm(double, int, double *, double *, double *, double *, int, int, double, double, double, double);
double **RaschSobLuchPlotTepPot(int, double *, double *, double *, double *, double *, char *, double ***);
double integroexponen(int, double);
double VIEW(int, int, double *, int);
void PrintMatr(double **, int);
double **preobMas(double *, int, int, int);
//-----
double *izstN(int izst, int kst, int l, int ocs, double *R, double *T, int N, char *snm, double **vf)
{
	int k=0, j=izst-kst; if (j<0) j=-j; 
	double *o=NULL;
	if (j<=N) o=podschchiele(izst, kst, ocs, R, T, snm, N, vf);
	else {
		o=new double[l]; if (!o) { cout << snm << endl; k=getchar(); exit(1); }
		for (k=0; k<l; k++) o[k]=0.0; }
	return o;
}
double ***zadrkt(int kost, char *snm, int N, double *Ref, double *Tra, double ht, double hv, double ha, double hb, int fl)
{ //расчет коэффициентов излучения, дошедшего от j до k-ой стенки
		int k=0, kst=0, q=3, m=4, b=0, f=m*kost*kost, s=0, p=0, x=0, y=0, z=0, v=0, j=0;
		double hf=1e0, *prs=new double[f], *pr=new double[m], Er=0.0, t=0.0, w=0.0;
		double **vf=NULL, *arg=NULL, *ptr=NULL, ***vm=new double**[m];
		if ((!prs) || (!pr)) { cout << snm << endl; k=getchar(); exit(1); }
		//-----
		if (fl>0) { q=38; p=3; x=0; y=x+1; z=y+1; 
			vf=new double*[kost]; arg=new double[p];
			if ((!vf) || (!arg)) { cout << snm << endl; k=getchar(); exit(1); }
				for (k=0; k<kost; k++) {
					ptr=new double[kost]; if (!ptr) { cout << snm << endl; k=getchar(); exit(1); } 
					vf[k]=ptr; for (j=0; j<kost; j++) ptr[j]=0.0; }
				for (k=0; k<kost; k++) {
					for (j=k+1; j<kost; j++) {
						v=k-j; if (v<0) v=-v; w=hv;
						if (v>1) { s=0; while (s<v) { w=w+(hv+ht); s++; }
						arg[x]=ha; arg[y]=hb; arg[z]=w;
						t=VIEW(q, p, arg, v); } else t=hf;
						vf[k][j]=t; vf[j][k]=t; } } }
		//-----
	for (k=0; k<f; k++) prs[k]=0.0;
		q=0; for (kst=1; kst<=kost; kst++) { //конечная стенка
		for (b=0; b<m; b++) pr[b]=0.0;
			for (k=1; k<=kost; k++) { //излучающая стенка
				pr=izstN(k, kst, m, kost, Ref, Tra, N, snm, vf);
				for (b=0; b<m; b++) { prs[q+b]=pr[b]; } q+=m; } }
		for (k=0; k<m; k++) vm[k]=preobMas(prs, m, kost, k); //PrintMatr(vf, kost); kst=m; for (b=0; b<kst; b++) { cout << "Matritsa " << b << endl; PrintMatr(vm[k], kost); }
		//-----
		if (pr) delete[]pr; if (arg) delete[]arg; if (prs) delete[]prs; //if ((vf) && (fl>0)) { PrintMatr(vf, kost); for (k=0; k<m; k++) PrintMatr(vm[k], kost); t=0.0; p=0; x=p+1; y=x+1; z=y+1; for (k=0; k<kost; k++) { for (j=0; j<kost; j++) { t=vm[p][k][j]+vm[x][k][j]-vm[y][k][j]-vm[z][k][j]; cout << "A ( " << k << " , " << j << " ) = " << t << "\t"; } cout << endl; } k=getchar(); }
	if (vf) { for (k=0; k<kost; k++) { ptr=vf[k]; if (ptr) delete[]ptr; } if (vf) delete[]vf; }
	return vm;
}
double **preobMas(double *vhma, int m, int n, int f)
{
	int k=0, j=0, s=f;
	double **matr=new double*[n], *p=NULL;
	for (k=0; k<n; k++) {
		p=new double[n];
		for (j=0; j<n; j++) { p[j]=vhma[s]; s+=m; }
		matr[k]=p; }
	return matr;
}
double *podschchiele(int no, int kste, int ocs, double *Ref, double *Tra, char *snm, int N, double **vf)
{
	int kol=2, *kl=new int[N], p=0, pe=p, ee=1, qs=p, kk=p, k=p, f=0, *ns=new int[N], nn=0;
	int j=p, jj=1, w=4, b=p, bb=p, gg=1, *st=new int[N], *ot=new int[N];
	double pp=0.0, *na=new double[w], nat=pp, hf=1e0, r=pp, eps=1e-6;
	struct derevo *prev=NULL, *prevv=NULL, *roo=NULL, *poir=NULL, *poil=NULL, *poi=NULL;
	struct derevo **pruk=new struct derevo*[ee], **ceuk=new struct derevo*[kol], **unpe=new struct derevo*[N];
	prev=new struct derevo;
	if ((!na) || (!kl) || (!prev) || (!pruk) || (!ceuk) || (!st) || (!ot) || (!unpe)) 
	{ cout << snm << endl; p = getchar(); exit(1); }
	kk=4; for (j=0; j<kk; j++) na[j]=0.0;
	kol=0; kk=kol; prev->lev=gg; gg++; prev->ste=no; prev->vis=1; prev->next=NULL; roo=prev; kol++; 
	prev->otre=-1; prev->back=NULL; bb=0; pruk[bb]=prev; unpe[b]=prev; b++; kl[k]=kol-kk; k++;
	//-----------
	bb=2; for (j=0; j<bb; j++) ceuk[j]=NULL; kk=kol; pe=0;
	poir=new struct derevo; poil=new struct derevo; 
	if ((!poir) || (!poil)) { cout << snm << endl; j=getchar(); exit(1); }
	poil->back=prev; poil->otre=0; poir->lev=gg; poil->lev=gg; 
	poir->back=prev; poir->otre=1; poir->ste=prev->ste+1; poil->ste=prev->ste-1;
	poir->vis=1; if ((poir->ste<1) || (poir->ste>ocs)) poir->vis=0; 
	if (poir->vis==1) { ceuk[pe] = poir; kol++; pe++; }
	else delete poir;
	poil->vis=1; if ((poil->ste<1) || (poil->ste>ocs)) poil->vis=0; 
	if (poil->vis == 1) { ceuk[pe] = poil; kol++; pe++; }
	else delete poil;
	j=0; poi=ceuk[j]; prev=unpe[b-1]; gg++;
	while ((poi) && (j<pe)) {
		if (poi->vis==1) {
			prev->next=poi;
			prev=poi; }
		j++; poi=ceuk[j]; }
	unpe[b]=prev; b++;
	kl[k]=kol-kk; k++;
	if (pruk) delete[]pruk; pruk=new struct derevo*[pe]; if (!pruk) { cout << snm << endl; j=getchar(); exit(1); }
	qs=0; for (j=0; j<pe; j++) { pruk[j]=NULL; poi=ceuk[j]; if ((poi) && (poi->vis==1)) { pruk[qs]=poi; qs++; } } 
	jj=qs; if (ceuk) delete[]ceuk;
	//------------//0 - пропускание (left), 1 - отражение (right)
	pp=3.0; 
	for (p=3; p<=N; p++) {
		w=0; r=eps; while (r<pow(2e0, pp-hf)) { w++; r=r+hf; } //w=floor(pow(2e0, pp-hf)); 
		bb=2*pe; pp=pp+hf;
		if (bb>0) { ceuk=new struct derevo*[bb]; 
		if (!ceuk) { cout << snm << endl; j=getchar(); exit(1); } 
		else for (j=0; j<bb; j++) ceuk[j]=NULL; }
		pe=0; qs=0; ee=0; kk=kol;
		while (ee<w) {
			if ((qs<jj) && (pruk) && (ceuk)) {
				prev=pruk[qs];
				if (prev) {
					prevv=prev->back;
					if (prevv) {
						poir=new struct derevo; poil=new struct derevo; 
						if ((!poir) || (!poil)) { cout << snm << endl; j=getchar(); exit(1); }
						poir->back=prev; poir->otre=1; poil->back=prev; poil->otre=0; poir->next=NULL; 
						poil->next=NULL; poir->lev=gg; poil->lev=gg;
						if (prev->otre==1) { if (prevv->ste>prev->ste) { poil->ste=prev->ste-1; poir->ste=prev->ste+1; } //отражение
						else { poil->ste=prev->ste+1; poir->ste=prev->ste-1; } }
						if (prev->otre==0) if (prevv->ste>prev->ste) { poil->ste=prev->ste-1; poir->ste=prev->ste+1; } //пропускание
						else { poil->ste=prev->ste+1; poir->ste=prev->ste-1; }
						if (prev->otre==-1) if (prevv->ste>prev->ste) { poil->ste=prev->ste-1; poir->ste=prev->ste+1; } //излучение j-ой стенки на соседние стенки
						else { poil->ste=prev->ste+1; poir->ste=prev->ste-1; }
						poir->vis=1; if ((poir->ste<1) || (poir->ste>ocs) || (prev->vis==0)) poir->vis=0; //проверка на видимость, что мы не вышли за пределы многослойной стенки, иначе - удаляем
						if (poir->vis==1) { ceuk[pe]=poir; kol++; pe++; } //kol - число видимых элементов
						else delete poir;
						poil->vis=1; if ((poil->ste<1) || (poil->ste>ocs) || (prev->vis==0)) poil->vis=0; //сначала вправо, потом влево
						if (poil->vis==1) { ceuk[pe]=poil; kol++; pe++; } //в ceuk сначала записываем все, потом в pruk - видимые
						else delete poil;
						ee+=2; qs++; } } } 
			ee++; }
		if (pruk) delete[]pruk; 
		pruk=new struct derevo*[pe]; if (!pruk) { cout << snm << endl; j=getchar(); exit(1); }
		qs=0; for (j=0; j<pe; j++) { poi=ceuk[j]; if ((poi) && (poi->vis==1)) { pruk[qs]=poi; qs++; } } jj=qs;
		kl[k]=kol-kk; k++;
		j=0; poi=pruk[j]; prev=unpe[b-1];
		while ((poi) && (j<pe)) {
			if (poi->vis==1) {
				prev->next=poi; prev=poi; } //спускаемся вниз по дереву
			j++; poi=pruk[j]; } 
		unpe[b]=prev; b++; gg++; 
		if (ceuk) delete[]ceuk; }
	//-----
	bb=1; poi=roo;
	for (k=1; k<=N; k++) {
		kk=bb-1; b=0;
		while ((b<kl[kk]) && (poi)) {
			if (poi->lev==bb) {
				for (w=0; w<N; w++) { ot[w]=0; st[w]=0; }
				prev=poi; nn=0;
				for (w=0; w<bb; w++) {
					if (poi) {
						st[kk-w]=poi->ste;
						ns[nn]=st[kk-w];
						ot[kk-w]=poi->otre;
						poi=poi->back; nn++; } }
				f=1; if (vf) { w=0; j=ns[w]; gg=j; for (w=1; w<nn; w++) { pe=ns[w]; if (pe>j) j=pe; if (pe<gg) gg=pe; }
				if (abs(j-gg)==1) f=-1; } //между соседними не считаем матрицы
				poi=prev; w=0; if (f>0) nat=rasotprpovs(st, ot, kk, kste, no, bb, Ref, Tra, vf); else nat=0.0;
				if (nat>0.0) { //cout << "nat = " << nat << "\t";
					if (st[kk]==kste) {
						if (st[kk-1]<kste) w=0; else w=1; 
						if (kste==1) w=1; 
						if (kste==ocs) w=0; 
						na[w]=na[w]+nat;
						w=2; j=1; if (st[j]<no) w=2; else w=3; 
						if (no==1) w=3; if (no==ocs) w=2; na[w]=na[w]+nat; //3 - излучение идет вправо от излучающей стенки, 2 - влево
					} } //0 - lr (излучение падает на рассматриваемую стенку слева, т.е. идет слева - направо), 1 - rl (излучение от других стенок идет справа, т.е. распространяется справа - налево)
				poi=poi->next; b++; } } 
		bb++; }
	poi=roo; while (poi)  { prev=poi->next; delete poi; poi=prev; }
	if (ot) delete[]ot; if (st) delete[]st; if (pruk) delete[]pruk; 
	if (kl) delete[]kl; if (ns) delete[]ns; if (unpe) delete[]unpe;
	return na;
}
double rasotprpovs(int *stt, int *ott, int gg, int kst, int izlst, int kost, double *R, double *T, double **vf)
{
	int e=1, ee=0, f=1, w=0, b=0, v=0; 
	double hf=1e0, po=hf, pot=hf, dk=hf, r=0.0;
	if ((stt[gg]==kst) && (kost>1)) {
		while (e<=gg) {
			w=stt[ee]; if (izlst==stt[b]) { 
				if (f>0) { e++; ee++; f=-1; continue; } }
			pot=otrprovs(w-1, ott[e], R, T, v, vf, stt, ee);
			po=po*pot; ee++; e++; } }
	else po=0.0;
	if ((vf) && (v==1)) { e=izlst-1; ee=kst-1; dk=vf[e][ee]; }
	return po*dk;
}
double otrprovs(int ste, int otr, double *R, double *T, int v, double **vf, int *ms, int ns)
{
	double ko=1e0, dk=ko;
	int k=0, j=0;
	if (otr==1) ko=R[ste];
	else if (!otr) ko=T[ste];
	if ((!v) && (vf) && (ns>0)) { k=ms[ns]-1; j=ms[ns-1]-1; dk=vf[j][k]; }
	return (ko*dk);
}
double **RaschVneshIzluch(double *Tra, double *Ref, double ***prot, int kost, double hlr0, double hrl0, char *snm)
{
	int k = 3, j=0, m=0, q = 2; 
	double *khrl=new double[kost], *khlr=new double[kost], *ri=new double[q], **mu=new double*[k], e=1e-20;
	double *hrl1 = new double[kost], *hlr1 = new double[kost], hf=1e0;
	double **mu0=NULL, **mu1=NULL, **mu2=NULL, **mu3=NULL;
	double w0=0.0, w1=w0, w2=w0, w3=w0, wk=w0, slr=w0, srl=w0, tlr=w0, trl=w0, tt=w0;
	if ((!khrl) || (!khlr) || (!ri) || (!hrl1) || (!hlr1) || (!mu)) { cout << snm << endl; k = getchar(); exit(1); }
	k=0; mu0=prot[k]; k++; mu1=prot[k]; k++; mu2=prot[k]; k++; mu3=prot[k];
	for (k = 0; k < kost; k++) { khrl[k] = 0.0; khlr[k] = 0.0; hrl1[k] = 0.0; hlr1[k] = 0.0; } for (k = 0; k < q; k++) ri[k] = 0.0;
	k=0; khlr[k] = hf; j = kost - 1; khrl[j] = hf;
	for (k = 1; k < kost; k++) {
		khlr[k] = khlr[k - 1] * Tra[k - 1];
		khrl[j - 1] = khrl[j] * Tra[j];
		j--; } //for (k=0; k<kost; k++) cout << "khrl " << k << " = " << khrl[k] << "\tkhlr " << k << " = " << khlr[k] << endl;
	q = 0; slr = 0.0; srl = 0.0; 
	for (j = 0; j < kost; j++) { //узнаем, какая доля внешнего излучения дошла до j-ой стенки - учитываем вклад от всех стенок
		hlr1[j] = 0.0; hrl1[j] = 0.0; //внешнее излучение падает с обеих сторон
		tlr = 0.0; trl = 0.0; m = q;
		for (k = 0; k<kost; k++) {
			w0=mu0[j][k]; 
			w1=mu1[j][k]; 
			w2=mu2[j][k]; 
			w3=mu3[j][k]; 
			wk = w2 + w3; //для внешнего излучения, распространяющегося слева направо
			tt = Ref[k] * w2 + Tra[k] * w3; //доля энергии, отраженная от k-ой стенки, и прошедшая через нее, пришла справа к j-ой стенке
			if (fabs(wk)>e) {
				tlr = tlr + tt*khlr[k] * w0 / wk; //доля энергии, отраженная от k-ой стенки, и прошедшая через нее, пришла слева к j-ой стенке
				trl = trl + tt*khlr[k] * w1 / wk;
			} //для учета при подсчете лучистой энергии, распространяющейся справа налево
	} q = m + 4 * kost - 1; m = q; //считаем, что стенки, до которых дошло излучение, как бы "сами" излучают, из-за того, что на них попало внешнее излучение
		for (k = (kost - 1); k >= 0; k--) {
			w0=mu0[j][k]; 
			w1=mu1[j][k]; 
			w2=mu2[j][k]; 
			w3=mu3[j][k]; 
			wk = w2 + w3; //для внешнего излучения, распространяющегося справа налево
			tt = Ref[k] * w3 + Tra[k] * w2; //деление на два потока
			if (fabs(wk) > e) {
				tlr=tlr+tt*khrl[k]*w0/wk; //для учета лучистой энергии, распространяющейся слева направо, падает слева от j-ой стенки
				trl=trl+tt*khrl[k]*w1/wk;
			} //то, что попадает справа налево на стенку, справа от j-ой стенки
		} 
		hlr1[j]=tlr*hlr0; //узнаем величину в абсолютных единицах этой лучистой энергии - слева от пластинки, падает слева направо
		hrl1[j]=trl*hrl0; //справа от пластинки, падает справа налево
	} 
	k=kost-1; slr = slr + Tra[k] * hlr1[k] + hrl0*Ref[k]; k=0; ri[k] = slr; //то, что вышло из последней стенки вправо
	k=0; srl = srl + Tra[k] * hrl1[k] + Ref[k] * hlr1[k]; k=1; ri[k] = srl; //то, что вышло из первой стенки влево
	k=0; hlr1[k] = hlr0; k=kost-1; hrl1[k] = hrl0; //for (k=0; k<ks; k++) cout << "hrl " << k << " = " << hrl1[k] << "\thlr " << k << " = " << hlr1[k] << endl;
	if (khlr) delete[]khlr; if (khrl) delete[]khrl; //cout << "Ras Vn Iz Sha k" << endl;
	k = 0; mu[k] = hlr1; k++; mu[k] = hrl1; k++; mu[k] = ri; 
	return mu;
}
double *RasIzlSerStenNacSok(double *Refl, double *Tran, double *Abso, double ***mpov, double ***mpos, double temc, double hlr0, double hrl0, double poristo, 
	double dpctp, double htch, double hvoz, int kost, double *ete, double *kttkve, int cemve, double *ktpvove, double *tevove, char *snm, 
	double *stchsrve, int dmkvove, double qobvee, double hko, double p1, double p2, double tx, double dkosups, double dkokpiko, int vybves, 
	int vybmar, int vp,	double *ktr, int dmktr, double epstem, int nomst, int indic, int vfqi, int kolPok)
{
	int j=2, k=0, f=6, *nost=new int[j], ki=0, m=3, u=1, x=4, y=5, v=0, q=0;
	double **mqi=new double*[kost-1], *qi=NULL, hf=1e0, d=0.0, w2=d, w3=d, w0=d, w1=d, wo=d;
	double **mu=NULL, mo=-hf, *Ts=NULL, *Tss=NULL, *hrl1=NULL, *hlr1=NULL, *reiz=NULL, qn=d;
	double ra=d, alks=d, *po=NULL, *sil=NULL, *sir=NULL, silv=0.0, sirv=0.0;
	double **mu0=NULL, **mu1=NULL, **mu2=NULL, **mu3=NULL;
	for (j=0; j<kost; j++) d=d+hf; 
	double *Rs=new double[f], *Ap=new double[f], *Tpctve=new double[kost]; //cout << "dp = " << dpctp << "\tht = " << htch << "\thv = " << hvoz << "\t";
	double *qir=new double[kost], *qil=new double[kost], e=1e-10, dt=0.0; //cout << "por = " << poristo << "\thlr0 = " << hlr0 << "\thrl0 = " << hrl0 << endl;
	double *tvs=new double[f], *epst=new double[f], *hvi=new double[f]; //cout << "d = " << d << "\tqob = " << qobvee << "\tkost = " << kost << "\t";
	double slr=0.0, srl=0.0, h=hvoz, w=p2*h, l=p1*h, t=0.0, r=t, ws=r, s=r;
	if ((!tvs) || (!epst) || (!hvi) || (!Rs) || (!Ap) || (!mqi) || (!qil) || (!qir) || (!Tpctve)) 
	{ cout << snm << endl; k = getchar(); exit(1); }
	for (j=0; j<f; j++) { tvs[j]=0.0; epst[j]=0.0; hvi[j]=0.0; Rs[j]=0.0; Ap[j]=0.0; }
	for (j=0; j<kost; j++) { qir[j]=0.0; qil[j]=0.0; Tpctve[j]=0.0; }
	//----- //j=0; cout << "u = " << u << "\tdko_1 = " << dkosups << "\tdko_2 = " << dkokpiko << "\tvybves = " << vybves << "\tvm = " << vybmar << endl;
	if (hlr0<e) { t=hvoz*(d-hf)+d*htch;
	dt=opredTempStenFragm(tx, 2*kost, ktpvove, tevove, ete, kttkve, cemve, dmkvove, htch, hvoz, qobvee, mo); 
	j=0; hlr0=ktr[j]*dt/t; } 
	Tpctve=opredTempTvKarFragm(kost, ete, kttkve, cemve, htch, hvoz, dpctp*qobvee, temc, mo, snm, tx, Tpctve); //for (k=0; k<kost; k++) cout << "t_v_s ( " << k << " ) = " << Tpctve[k] << "\t"; cout << endl;
	k=0; Ts=opredTempStenFragmMasTem(2*kost, ktpvove, tevove, ete, kttkve, cemve, dmkvove, htch, hvoz, (hf-dpctp)*qobvee, temc, mo, snm, kost, k); 
	k++; Tss=opredTempStenFragmMasTem(2*kost, ktpvove, tevove, ete, kttkve, cemve, dmkvove, htch, hvoz, (hf-dpctp)*qobvee, temc, mo, snm, kost, k); //for (k=0; k<kost; k++) cout << "t_l_s ( " << k << " ) = " << Ts[k] << "\tt_p_s = " << Tss[k] << "\t"; cout << endl;
	if ((epstem>e) && (indic>=0)) { if (!(indic)) Ts[nomst]=Ts[nomst]+epstem; else Tss[nomst]=Tss[nomst]+epstem; }
	mu=RaschVneshIzluch(Tran, Refl, mpov, kost, hlr0, hrl0, snm); 
	k=0; hlr1=mu[k]; k++; hrl1=mu[k]; k++; reiz = mu[k]; if (mu) delete[]mu; 
	k=0; slr=reiz[k]; k++; srl=reiz[k]; if (reiz) delete[]reiz; //for (k=0; k<kost; k++) cout << "hlr1 ( " << k << " ) = " << hlr1[k] << "\thrl1 = " << hrl1[k] << "\t"; cout << endl;
	mu=RaschSobLuchPlotTepPot(kost, Ts, Tss, Tran, Refl, Abso, snm, mpos);
	k=0; sir=mu[k]; k++; sil=mu[k]; k++; reiz=mu[k]; if (mu) delete[]mu;
	k=0; sirv=reiz[k]; k++; silv=reiz[k]; if (reiz) delete[]reiz; //for (k=0; k<kost; k++) cout << "silr ( " << k << " ) = " << sir[k] << "\tsirl = " << sil[k] << "\t"; cout << endl; k=getchar();
	k=0; ki=k; Ap[k]=w*l; nost[ki]=k; ki++; //k=0
	k++; u=k; Ap[k]=h*l; //u=1
	k++; Ap[k]=Ap[ki]; nost[ki]=k; //k=2
	k++; m=k; Ap[m]=Ap[u]; //m=3
	k++; x=k; Ap[x]=w*h; //x=4
	k++; y=k; Ap[y]=Ap[x]; //y=5 //k=4*kost*(kost-1); ki=4*kost*kost; for (j=k; j<ki; j++) { if (!(j%4)) cout << endl; cout << "w ( " << j << " ) = " << mpo[j] << "\t"; }
	for (j=0; j<(kost-1); j++) {
		for (k=0; k<f; k++) tvs[k]=Tpctve[j]; k=0; tvs[k]=Ts[j]; k=2; tvs[k]=Tss[j]; //0 - более высокая температура у стенки, 2 - менее высокая
		for (k=0; k<f; k++) { t=opredKTPTKToch(stchsrve, ete, tvs[k], cemve); 
		if (t<0.0) t=0.0; if (t>hf) t=hf; epst[k]=t; }
		for (k=0; k<f; k++) Rs[k]=hf-epst[k];
		for (k=0; k<f; k++) hvi[k]=0.0; k=0; hvi[k]=hrl1[j+1]+sil[j+1]; k=2; hvi[k]=hlr1[j]+sir[j]; 
	if (!vfqi) qi=SeryeStenkiRasIzl(w, h, l, tvs, epst, hvi, Rs, Ap, nost, f, snm, v);  //x=0; y=x+1; for (k=0; k<f; k++) cout << "\tte ( " << k+1 << " ) = " << tvs[k] << "\teps = " << epst[k] << "\tAp = " << Ap[k] << "\tq_1 = " << qi[x] << "\tq_3 = " << qi[y] << endl; cout << endl;
	else if (vfqi==1) qi=SeryeStenkiRasIzlDifPov(w, h, l, tvs, epst, hvi, Ap, nost, f, snm); //N серых диффузных поверхностей
	mqi[j]=qi;
	k=0; qil[j+1]=-qi[k]; //l - выходит со следующей и приходит на текущую слева (вправо)
	k++; qir[j]=-qi[k]; } //r - выходит с текущей стенки и приходит справа (влево) на следующую
	k=0; qil[k]=-(srl+silv); //приходит слева (направо) на первую стенку
	k=kost-1; qir[k]=-(slr+sirv); //приходит справа (налево) на последнюю стенку //for (k=0; k<kost; k++) cout << "ql ( " << k << " ) = " << qil[k] << "\tqr = " << qir[k] << "\t"; cout << endl;
	for (k=0; k<(kost-1); k++) { qi=mqi[k]; if (qi) delete[]qi; } if (mqi) delete[]mqi;
	if (tvs) delete[]tvs; if (epst) delete[]epst; if (hvi) delete[]hvi; //for (k=0; k<f; k++) cout << "Ap ( " << k << " ) = " << Ap[k] << "\t";
	if (nost) delete[]nost; if (Rs) delete[]Rs; if (Ap) delete[]Ap; if (Tpctve) delete[]Tpctve; //cout << "w = " << w << "\th = " << h << "\tl = " << l << endl; 
	//-----
	if ((!vp) || (vp==1)) { //уточнение коэффициента поглощения
	slr=0.0; srl=0.0; k=0; ki=kost-1; x=2; y=x+1;
	mu2=mpov[x]; mu3=mpov[y];
	for (j=0; j<kost; j++) { 
		r=qil[j]; //то, что идет влево
		s=qir[j]; //то, что идет вправо
		w2=mu2[ki][j];
		w3=mu3[ki][j];
		wo=w2+w3;
		if (fabs(wo)>e) slr=slr+(w2*r+s*w3)/wo;
		w2=mu2[k][j]; //влево
		w3=mu3[k][j]; //вправо
		wo=w2+w3;
		if (fabs(wo)>e) srl=srl+(w2*r+s*w3)/wo; }
	t=hvoz*(d-hf)+d*htch; qn=hlr0-hrl0; ra=slr-srl; //dt=opredTempStenFragm(tx, 2*kost, ktpvove, tevove, ete, kttkve, cemve, dmkvove, htch, hvoz, qobvee, mo); j=0; qn=ktr[j]*dt/t; qn=qn-hrl0; //j=0; cout << "te_KTR = " << tx << "\tktr = " << ktr[j] << "\tqin = " << qn << "\tsrl = " << srl << "\tslr = " << slr << "\tra = " << ra << "\tq_v_lr = " <<  hlr0 << "\tq_v_rl = " << hrl0 << "\ttol = " << t << "\tal = " << alks << endl;
	if ((ra>e) && (qn>e)) alks=-log(ra/qn)/t; else alks=0.0; }
	//-----
	if ((!vp) || (vp==1)) { k=1; po=new double[k]; if (!po) { cout << snm << endl; k=getchar(); exit(1); }
	k=0; if (!vp) po[k]=alks; else po[k]=ra; } //for (k=0; k<kost; k++) cout << "ql ( " << k << " ) = " << qil[k] << "\tqr = " << qir[k] << "\t"; cout << endl; k=getchar();
	if (vp==2) { po=qir; if (qil) delete[]qil; } 
	else if (vp==3) { po=qil; if (qir) delete[]qir; }
	//-----
	return po;
}
double **RaschSobLuchPlotTepPot(int kost, double *Ts, double *Tss, double *Tra, double *Ref, double *Ab, char *snm, double ***prot)
{
	int j=0, q=3, k=0, x=0, y=0, z=0, v=0, s=0, p=2;
	double tlr=0.0, trl=0.0, ta=0.0, tt=0.0, tb=0.0, w0=0.0, w1=0.0, w2=0.0, slr=0.0, srl=0.0;
	double w3=0.0, wk=0.0, sig=5.67e-8, **mu=new double*[q], ra=0.0;
	double *sislr=new double[kost], *sisrl=new double[kost], *ao=new double[p];
	double **mu0=NULL, **mu1=NULL, **mu2=NULL, **mu3=NULL, e=1e-20;
	k=0; mu0=prot[k]; k++; mu1=prot[k]; k++; mu2=prot[k]; k++; mu3=prot[k];
	if ((!sisrl) || (!sislr) || (!mu)) { cout << snm << endl; k=getchar(); exit(1); }
	q=0; x=q+1; y=x+1; z=y+1; v=z+1;
	for (j=0; j<kost; j++) { //узнаем, какая энергия собственного излучения дошла до j-ой стенки
		sislr[j]=0.0; sisrl[j]=0.0; tlr=0.0; trl=0.0;
		for (k=0; k<kost; k++) {
			w0=mu0[j][k]; 
			w1=mu1[j][k]; 
			w2=mu2[j][k]; 
			w3=mu3[j][k]; 
			wk=w2+w3;
			tt=pow(Ts[k], 4e0)*w2+pow(Tss[k], 4e0)*w3;
			if (fabs(wk)>e) {
				ta=sig*Ab[k]*w0*tt/wk; tlr=tlr+ta; // то, что попадает слева направо
				tb=sig*Ab[k]*w1*tt/wk; trl=trl+tb; } // то, что попадает справа налево
		}
		sislr[j]=tlr; //собственное излучение стенок слева направо и справа налево
		sisrl[j]=trl;
	}
	k=kost-1; slr=slr+Tra[k]*sislr[k]; slr=slr+Ab[k]*sig*pow(Tss[k], 4e0); //из последней стенки
	k=0; srl=srl+Tra[k]*sisrl[k]; srl=srl+Ab[k]*sig*pow(Ts[k], 4e0); //из первой стенки
	k=0; sislr[k]=0.0; k=kost-1; sisrl[k]=0.0; //стенок вне многослойной стенки нет, излучение ниоткуда не падает, есть только внешнее излучение 
	ra=slr-srl; //cout << "slr = " << slr << "\tsrl = " << srl << "\tra = " << ra << "\t"; //излучение за крайней правой стенкой и за крайней левой
	k=0; ao[k]=slr; k++; ao[k]=srl;
	k=0; mu[k]=sislr; k++; mu[k]=sisrl; k++; mu[k]=ao; 
	return mu;
}