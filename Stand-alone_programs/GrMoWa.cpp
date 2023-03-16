#include <iostream>
#include <math.h>
#include <conio.h>
#include <cstdio>
#include <stdlib.h>
using namespace std;
const int N=6;
const double pi = acos(-1.0), epsmin=1e-15;

double *GRAYDIFFSPEC(int, double *, double *, double *, double *,double **, int *, double *, double *);
double VIEW(int, int, double *);
double perppltf(double, double, double, double, double, double, double);
double parlpltf(double, double, double, double, double, double, double);
double *GAUSS(double **, double *, double *);
void SeryeStenkiRasIzl(double, double, double, double *, double *, double *, double *, double *);
void main()
{
	int j, f=N;
	double *tvs=new double[f], *epst=new double[f], *hvi=new double[f], *Rs=new double[f], *A=new double[f]; 
	if ((!tvs) || (!epst) || (!hvi) || (!Rs) || (!A)) { cout << "No memory!" << endl; j=getchar(); exit(1); } 
	else for (j=0; j<f; j++) { tvs[j]=0.0; epst[j]=0.0; hvi[j]=0.0; Rs[j]=0.0; A[j]=0.0; }
	double w=0.4, h=0.3, l=1e6;
	for (j=0; j<f; j++) tvs[j]=6e2; tvs[0]=1e3; tvs[2]=1e3;
epst[0]=0.3; Rs[0]=0.7; epst[1]=0.8; Rs[1]=0.2; epst[2]=0.3; Rs[2]=0.0;
epst[3]=0.8; Rs[3]=0.0; epst[4]=0.8; Rs[4]=0.0; epst[5]=0.8; Rs[5]=0.0;
		A[0]=w*l; A[2]=w*l; A[1]=h*l; A[3]=h*l; A[4]=h*w; A[5]=h*w; 	
	SeryeStenkiRasIzl(w, h, l, tvs, epst, hvi, Rs, A);
	delete []tvs; delete []epst; delete []hvi; delete []Rs; j=getchar();
}
void SeryeStenkiRasIzl(double ww, double hh, double ll, double *T, double *EPS, double *HOs, double *RHOs, double *A)
{ int i, j, f=3; double **Fs=new double*[N];
double *PIN=new double[N], *POUT=new double[N], *arg=new double[f], *p, t, r, b;
for (i=0; i<N; i++) { p=new double[N]; for (j=0; j<N; j++) p[j]=0.0; Fs[i]=p; }
	double sumq=0.0, sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	int *id=new int[N], iclsd; for (i=0; i<N; i++) id[i]=1; // Surfaces 1&3 (bottom and top)
    	for (i=0; i<N; i++) if ((id[i])==0) PIN[i]=q[i]; else PIN[i]=sigma*pow(T[i],4.0); //Заполнение массива PIN ПТП q и температурой T - перевод температур в ППИ, Fs(1,2)=F1-2, 
	iclsd=1; //угловые коэффициенты; конфигурация замкнута (iclsd=1), диагональные элементы не нужны
	arg[0]=fabs(hh); arg[1]=fabs(ll); arg[2]=fabs(ww); //h,l,w
	Fs[0][1]=VIEW(39,f,arg); //Fs(1,2)=F1-2, F12
    	arg[0]=fabs(ww); arg[1]=fabs(ll); arg[2]=fabs(hh);     // a,b,c
		r=VIEW(38,f,arg);
		t=parlpltf(ww,ww,2.0*ww,ll,0.0,ll,hh);
		b=RHOs[1];
    	Fs[0][2]=r+b*t; //Fs1-3=F1-3+rhos2*F1(2)-3, Fs13
    	Fs[0][3]=Fs[0][1]+RHOs[1]*perppltf(ww,2.0*ww,0.0,hh,ll,0.0,ll); //Fs1-4=F1-4+rhos2*F1(2)-4
    	Fs[0][4]=0.5*(1.0-(1.0-RHOs[1])*Fs[0][1]-Fs[0][2]-Fs[0][3]); //смещены два прямоугольника на ww по оси X
    	Fs[0][5]=Fs[0][4]; //Fs23
    	arg[0]=fabs(ww); arg[1]=fabs(ll); arg[2]=fabs(hh); // h,l,w
    	Fs[1][2]=VIEW(39,3,arg)+RHOs[0]*perppltf(hh,2.0*hh,0.0,ww,ll,0.0,ll);  // Fs2-3=F2-3+rhos1*F1(1)-3 //смещены два прямоугольника на hh по оси X, т.к. у них разные размеры
	   	arg[0]=fabs(ll); arg[1]=fabs(hh); arg[2]=fabs(ww); //a,b,c
    	Fs[1][3]=VIEW(38,3,arg)+RHOs[0]*parlpltf(hh,hh,2.0*hh,ll,0.0,ll,ww);  // Fs2-4=F2-4+rhos1*F2(1)-4, Fs24
    	Fs[1][4]=0.5*(1.0-(1.0-RHOs[0])*A[0]*Fs[0][1]/A[1]-Fs[1][2]-Fs[1][3]); // Fs25
    	Fs[1][5]=Fs[1][4];
    	arg[0]=fabs(ww); arg[1]=fabs(ll); arg[2]=2.0*fabs(hh);  // a,b,c
    	Fs[2][2]=RHOs[0]*(VIEW(38,3,arg)+RHOs[1]*parlpltf(ww,ww,2.0*ww,ll,0.0,ll,2.0*hh)); //Fs33=rhos1*F3(1)-3 + rhos1*rhos2*F3(12+21)-3 //излучение стенки саму на себя
    	Fs[2][3]=Fs[0][3]+RHOs[0]*(perppltf(0.0,ww,hh,2.0*hh,ll,0.0,ll)+RHOs[1]*perppltf(ww,2.0*ww,hh,2.0*hh,ll,0.0,ll)); // Fs3-4=Fs1-4+rhos1*F3(1)-3+rhos1*rhos2*F3(12+21)-4
    	Fs[2][4]=0.5*(1.0-(1.0-RHOs[0])*A[0]*Fs[0][2]/A[2]-(1.0-RHOs[1])*A[1]*Fs[1][2]/A[2]-Fs[2][2]-Fs[2][3]);
    	Fs[2][5]=Fs[2][4];
    	arg[0]=fabs(hh); arg[1]=fabs(ll); arg[2]=2.0*fabs(ww); // a,b,c
    	Fs[3][3]=RHOs[1]*(VIEW(38,3,arg)+RHOs[0]*parlpltf(hh,hh,2.0*hh,ll,0.0,ll,2.0*ww)); //Fs44=rhos2*F4(2)-4+rhos1*rhos2*F4(12+21)-4
    	Fs[3][4]=0.5*(1.0-((1.0-RHOs[0])*A[0]*Fs[0][3]+(1.0-RHOs[1])*A[1]*Fs[1][3]+A[2]*Fs[2][3])/A[3]-Fs[3][3]);
    	Fs[3][5]=Fs[3][4];
    	Fs[4][5]=0.5*(1.0-((1.0-RHOs[0])*A[0]*Fs[0][4]+(1.0-RHOs[1])*A[1]*Fs[1][4]+A[2]*Fs[2][4]+A[3]*Fs[3][4])/A[4]); //Fs56 //for (j=0; j<N; j++) cout << "A ( " << j << " ) = " << A[j] << endl; cout << "Fs" << endl; for (i=0; i<N; i++) { for (j=0; j<N; j++) printf("%0.4lf\t",Fs[i][j]); cout << endl; } 
	POUT=GRAYDIFFSPEC(iclsd,A,EPS,RHOs,HOs,Fs,id,PIN,POUT); //Solve system of equations //for (i=0; i<N; i++) cout << POUT[i] << endl;
    	sumq=0.0; //проверка полной ППИ
    	for (i=0; i<N; i++) //Вывод - перевод в температуры
		{ if (!(id[i])) T[i]=pow((POUT[i]/sigma), 0.25); else q[i]=POUT[i];
        	QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; 
			cout << "Surface " << i << "\tT [K] = " << T[i] << "\tq [W/m2] = " << q[i] << "\tQA [W] = " << QA[i]/ll << endl;}
for (i=0; i<N; i++) { p=Fs[i]; delete []p; } delete []Fs; delete []A; delete []PIN; 
delete []POUT; delete []arg; delete []q; delete []QA; delete []id; j=getchar(); }
double *GRAYDIFFSPEC(int iclsd, double *A, double *EPS, double *RHOs, double *HOs, double **Fs, int *ID, double *PIN, double *POUT) //Routine to fill view factor matrix and solve for missing surface temperatures/fluxes
{ int i, j; 
double **qm=new double*[N], **em=new double*[N], **pm=new double*[N], *B=new double[N], *p, *idf=new double[N], hf=1e0, s, ikr; // Compute missing view factors - Lower left triangle by reciprocity
for (i=0; i<N; i++) { p=new double[N]; qm[i]=p; } 
for (i=0; i<N; i++) { p=new double[N]; em[i]=p; } 
for (i=0; i<N; i++) { p=new double[N]; pm[i]=p; }
for (i=0; i<N; i++) { s=0.0; for (j=0; j<ID[i]; j++) s=s+hf; idf[i]=s; }
for (i=1; i<N; i++)
	 for (j=0; j<=(i-1); j++)
		  if (fabs(A[i])>0.0) Fs[i][j]=A[j]/A[i]*Fs[j][i]; //Для закрытой конфигурации, необходимо рассчитать диагональные элементы из правила суммирования
	if (iclsd==1) { for (i=0; i<N; i++) {
		Fs[i][i]=1e0;
		for (j=0; j<N; j++) {
			if (j!=i) {
				Fs[i][i]=Fs[i][i]-(1e0-RHOs[j])*Fs[i][j]; } }
		Fs[i][i]=Fs[i][i]/(1e0-RHOs[i]); } }
	 for (i=0; i<N; i++) { // Заполнение матриц коэффициентов q и e
		for (j=0; j<N; j++) {
			if (i==j) ikr=1e0; else ikr=0.0; // Kronecker delta_ij //cout << "i = " << i << "\tj = " << j << "\tk = " << ikr << endl;
			qm[i][j]=ikr/EPS[j]-((1.0-RHOs[j])/EPS[j]-1.0)*Fs[i][j];
			em[i][j]=ikr-(1.0-RHOs[j])*Fs[i][j]; } }
	 for (i=0; i<N; i++) { // Заполнение матрицы выходных коэффициентов POUT и RHS
			B[i]=-HOs[i];
			for (j=0; j<N; j++) {
					pm[i][j]=qm[i][j]*idf[j]-em[i][j]*(1e0-idf[j]);
					B[i]=B[i]+(em[i][j]*idf[j]-qm[i][j]*(1e0-idf[j]))*PIN[j]; } }
	 POUT=GAUSS(pm, B, POUT); //for (i=0; i<N; i++) cout << "pout = " << POUT[i] << endl; 
	 for (i=0; i<N; i++) { p=qm[i]; delete []p; p=em[i]; delete []p; p=pm[i]; delete []p; } 
	 delete []qm; delete []em; delete []pm; delete []B; return POUT; } //Инвертирование матрицы коэффициентов POUT и умножение на матрицу RHS, чтобы получить POUT
double VIEW(int NO, int NARG, double *ARG)
{ double VIEW=0.0, A, B, C, X, Y, RTX, RTY, RT, H, L, W, HH, WW, W2, H2, HW2, HW, H12, W12, C1, C2, C3, e=epsmin; int k;
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
{ double perppltf, A=0.0, F=0.0, *ARG=new double[3], e=epsmin;
int NARG=3;
      ARG[0]=Y2; ARG[1]=Z3; ARG[2]=X2;
      A=X2*Z3;
      if (fabs(A*Y2)>e) F=A*VIEW(39,NARG,ARG);
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
      perppltf=F/(2.0*(X2-X1)*Z1); //cout << "per = " << perppltf << endl;
	  delete []ARG; return perppltf; }
double parlpltf(double X1, double X2, double X3, double Y1, double Y2, double Y3, double C)
{ int NARG=3; double parlpltf, A=0.0, F=0.0, *ARG=new double[NARG], e=epsmin, t;
      ARG[0]=X3; ARG[1]=Y3; ARG[2]=C; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=A*VIEW(38,NARG,ARG); //cout << "F 1 = " << F << "\tA 1 = " << A << endl;
      ARG[1]=Y2; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG); //cout << "F 2 = " << F << "\tA 2 = " << A << endl;
      ARG[1]=Y3-Y1; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG); //cout << "F 3 = " << F << "\tA 3 = " << A << endl;
	  ARG[1]=Y2-Y1; A=ARG[0]*ARG[1]; t=VIEW(38,NARG,ARG); //printf("t = %0.15lf\n",t);
      if(fabs(A)>e) F=F+A*VIEW(38,NARG,ARG); // cout << "F 4 = " << F << "\tA 4 = " << A << endl;
      ARG[0]=X2; ARG[1]=Y3; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG); //cout << "F 5 = " << F << "\tA 5 = " << A << endl;
      ARG[1]=Y2; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F+A*VIEW(38,NARG,ARG); //cout << "F 6 = " << F << "\tA 6 = " << A << endl;
      ARG[1]=Y3-Y1; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F+A*VIEW(38,NARG,ARG); //cout << "F 7 = " << F << "\tA 7 = " << A << endl;
      ARG[1]=Y2-Y1; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG); //cout << "F 8 = " << F << "\tA 8 = " << A << endl;
      ARG[0]=X3-X1; ARG[1]=Y3; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG); //cout << "F 9 = " << F << "\tA 9 = " << A << endl;
      ARG[1]=Y2; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F+A*VIEW(38,NARG,ARG); //cout << "F 10 = " << F << "\tA 10 = " << A << endl;
      ARG[1]=Y3-Y1; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F+A*VIEW(38,NARG,ARG); //cout << "F 11 = " << F << "\tA 11 = " << A << endl;
      ARG[1]=Y2-Y1; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG); //cout << "F 12 = " << F << "\tA 12 = " << A << endl;
      ARG[0]=X2-X1; ARG[1]=Y3; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F+A*VIEW(38,NARG,ARG); //cout << "F 13 = " << "\tA 13 = " << A << endl;
      ARG[1]=Y2; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG); //cout << "F 14 = " << F << "\tA 14 = " << A << endl;
      ARG[1]=Y3-Y1; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F-A*VIEW(38,NARG,ARG); //cout << "F 15 = " << F << "\tA 15 = " << A << endl;
      ARG[1]=Y2-Y1; A=ARG[0]*ARG[1];
      if(fabs(A)>e) F=F+A*VIEW(38,NARG,ARG); //cout << "F 16 = " << F << "\tA 16 = " << A << endl;
      parlpltf=F/(4.0*(X1*Y1)); //cout << "par = " << parlpltf << endl;
      delete []ARG; return parlpltf; }
double *GAUSS(double **A, double *B, double *X)
{ int I, *L=new int[N], K, J, LK;
double *S=new double[N], SMAX, RMAX, R, XMULT, SUM;
      for (I=0; I<N; I++)
			{ L[I]=I; SMAX=0.0;
			for (J=0; J<N; J++)
				if (fabs(A[I][J])>SMAX) SMAX=fabs(A[I][J]);
			S[I]=SMAX; }
      for (K=0; K<N-1; K++) 
        { RMAX=0.0;
        for (I=K; I<N; I++) 
		{ R=fabs(A[L[I]][K])/S[L[I]];
		if (R<=RMAX) continue;
		else { J=I; RMAX=R; } }
		LK=L[J]; L[J]=L[K]; L[K]=LK;
		for (I=K+1; I<N; I++)
			{ XMULT = A[L[I]][K]/A[LK][K];
				for (J=K+1; J<N; J++)
					A[L[I]][J]=A[L[I]][J]-XMULT*A[LK][J];
				A[L[I]][K] = XMULT; } }
      for (K=0; K<N-1; K++)  
        for (I=K+1; I<N; I++)      
          B[L[I]] = B[L[I]] - A[L[I]][K]*B[L[K]];
      X[N-1] = B[L[N-1]]/A[L[N-1]][N-1];
      for (I=N-1; I>=0; I--) {
		  SUM = B[L[I]];       
		  for (J=I+1; J<N; J++)
          SUM = SUM - A[L[I]][J]*X[J];    
        X[I] = SUM/A[L[I]][I]; } 
	  delete []L; delete []S; return X; }