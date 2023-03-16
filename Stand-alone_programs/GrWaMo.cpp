#include <iostream>
#include <math.h>
#include <conio.h>
#include <cstdio>
#include <stdlib.h>

using namespace std;
const int N=6;
void GRAYDIFFSPEC(int,double [],double [],double [],double [],double [][N+1],int [],double [],double []);
double VIEW(int,int,double []);
double perppltf(double,double,double,double,double,double,double);
double parlpltf(double,double,double,double,double,double,double);
void GAUSS(double [][N+1],double [],double []);
void SeryeStenkiRasIzl(double, double, double, double *, double *, double *, double *);

void main()
{
	int j, f=6;
	double *tvs=new double[f], *epst=new double[f], *hvi=new double[f], *Rs=new double[f]; 
	if ((!tvs) || (!epst) || (!hvi) || (!Rs)) { cout << "No memory!" << endl; j=getchar(); exit(1); } 
	else for (j=0; j<f; j++) { tvs[j]=0.0; epst[j]=0.0; hvi[j]=0.0; Rs[j]=0.0; }
	double w=0.4, h=0.3, l=1e6;
	for (j=0; j<f; j++) tvs[j]=6e2; tvs[0]=1e3; tvs[2]=1e3;
epst[0]=0.3; Rs[0]=0.7; epst[1]=0.8; Rs[1]=0.2; epst[2]=0.3; Rs[2]=0.0;
epst[3]=0.8; Rs[3]=0.0; epst[4]=0.8; Rs[4]=0.0; epst[5]=0.8; Rs[5]=0.0;
	SeryeStenkiRasIzl(w, h, l, tvs, epst, hvi, Rs);
	delete []tvs; delete []epst; delete []hvi; delete []Rs; j=getchar();
}
void SeryeStenkiRasIzl(double w, double h, double l, double *tvs, double *epst, double *hvi, double *Rs)
{ double A[N+1],Fs[N+1][N+1],EPS[N+1],RHOs[N+1],HOs[N+1],PIN[N+1],POUT[N+1],arg[6],ww,hh,ll;
	double sumq,sigma=5.670E-8,T[N+1],q[N+1],QA[N+1];
	int id[N+1],iclsd,i;
	for (i=1; i<=N; i++) { EPS[i]=epst[i-1]; T[i]=tvs[i-1]; HOs[i]=hvi[i-1]; RHOs[i]=Rs[i-1]; id[i]=1; }
	
	// Dimensions
    	ww=w;
    	hh=h;
    	ll=l;
	
	// Surfaces 1&3 (bottom and top)
    	A[1]=ww*ll; 
	A[3]=ww*ll;          // ww x ll
    	//HOs[1]=0.; 
	//HOs[3]=0.;           // no external irradiation
    	//EPS[1]=0.3; 
	//RHOs[1]= 0.7; 
	//EPS[3]=0.3; 
	//RHOs[3]= 0.;
   	//id[1]=1; 
	//id[3]=1;             // T specified
    	//T[1]=1000.; 
	//T[3]=1000.;          // T in K
	
	// Surfaces 2&4 (left and right)
    	A[2]=hh*ll; 
	A[4]=hh*ll;         // hh x ll
    	//HOs[2]=0.; 
	//HOs[4]=0.;          // no external irradiation
    	//EPS[2]=0.8; 
	//RHOs[2]= 0.2; 
	//EPS[4]=0.8; 
	//RHOs[4]= 0.;
    	//id[2]=1; 
	//id[4]=1;            // T specified! Fs(1,2)=F1-2
    	//T[2]= 600.; 
	//T[4]= 600.;         // T in K
	
	// Surfaces 5&6 (front and back)
    	A[5]=hh*ww; 
	A[6]=hh*ww;         // hh x ww
    	//HOs[5]=0.; 
	//HOs[6]=0.;          // no external irradiation
    	//EPS[5]=0.8; 
	//RHOs[5]= 0.; 
	//EPS[6]=0.8; 
	//RHOs[6]= 0.;
    	//id[5]=1; 
	//id[6]=1;            // T specified
    	//T[5]= 600.;         // Fs(1,2)=F1-2
	//T[6]= 600.;         // T in K
	
	//fill PIN array with q ant T
		for (i=1; i<=N; i=i+1)
			{
			if (id[i]==0)
				PIN[i]=q[i];
			else 
				PIN[i]=sigma*pow(T[i],4);  //convert temperatures to emmissive powers
			}
	
	// view factors; since configuration is closed (iclsd=1), diagonal terms are not needed
	iclsd=1;
	
	//F12
	arg[1]=hh;
	arg[2]=ll;
	arg[3]=ww;		//h,l,w
	Fs[1][2]=VIEW(39,3,arg); //Fs(1,2)=F1-2
	
	// Fs13
    	arg[1]=ww; 
    	arg[2]=ll; 
    	arg[3]=hh;     // a,b,c
    	Fs[1][3]=VIEW(38,3,arg)+RHOs[2]*parlpltf(ww,ww,2.*ww,ll,0.,ll,hh);  // Fs1-3=F1-3+rhos2*F1(2)-3
    	Fs[1][4]=Fs[1][2]+RHOs[2]*perppltf(ww,2.*ww,0.,hh,ll,0.,ll);         // Fs1-4=F1-4+rhos2*F1(2)-4
    	Fs[1][5]=.5*(1.-(1.-RHOs[2])*Fs[1][2]-Fs[1][3]-Fs[1][4]);
    	Fs[1][6]=Fs[1][5];
	
	// Fs23
    	arg[1]=ww; 
	arg[2]=ll; 
	arg[3]=hh;    // h,l,w
    	Fs[2][3]=VIEW(39,3,arg)+RHOs[1]*perppltf(hh,2.*hh,0.,ww,ll,0.,ll);  // Fs2-3=F2-3+rhos1*F1(1)-3
	
	// Fs24
    	arg[1]=ll; 
	arg[2]=hh; 
	arg[3]=ww;     // a,b,c
    	Fs[2][4]=VIEW(38,3,arg)+RHOs[1]*parlpltf(hh,hh,2.*hh,ll,0.,ll,ww);  // Fs2-4=F2-4+rhos1*F2(1)-4
	
	// Fs25
    	Fs[2][5]=.5*(1.-(1.-RHOs[1])*A[1]*Fs[1][2]/A[2]-Fs[2][3]-Fs[2][4]);
    	Fs[2][6]=Fs[2][5];
	
	// Fs33=rhos1*F3(1)-3 + rhos1*rhos2*F3(12+21)-3    
    	arg[1]=ww; 
	arg[2]=ll; 
	arg[3]=2.*hh;  // a,b,c
    	Fs[3][3]=RHOs[1]*(VIEW(38,3,arg)+RHOs[2]*parlpltf(ww,ww,2.*ww,ll,0.,ll,2.*hh));
    	Fs[3][4]=Fs[1][4]+RHOs[1]*(perppltf(0.,ww,hh,2.*hh,ll,0.,ll)       // Fs3-4=Fs1-4+rhos1*F3(1)-3
                   +RHOs[2]*perppltf(ww,2.*ww,hh,2.*hh,ll,0.,ll));     //       +rhos1*rhos2*F3(12+21)-4  
    	Fs[3][5]=.5*(1.-(1.-RHOs[1])*A[1]*Fs[1][3]/A[3]-(1.-RHOs[2])*A[2]*Fs[2][3]/A[3]-Fs[3][3]-Fs[3][4]);
    	Fs[3][6]=Fs[3][5];
	
	// Fs44=rhos2*F4(2)-4 + rhos1*rhos2*F4(12+21)-4 
    	arg[1]=hh; 
	arg[2]=ll; 
	arg[3]=2.*ww;  // a,b,c
    	Fs[4][4]=RHOs[2]*(VIEW(38,3,arg)+RHOs[1]*parlpltf(hh,hh,2.*hh,ll,0.,ll,2.*ww));
    	Fs[4][5]=.5*(1.-((1.-RHOs[1])*A[1]*Fs[1][4]+(1.-RHOs[2])*A[2]*Fs[2][4]+A[3]*Fs[3][4])/A[4]-Fs[4][4]);
    	Fs[4][6]=Fs[4][5];
	
	// Fs56
    	Fs[5][6]=.5*(1.-((1.-RHOs[1])*A[1]*Fs[1][5]+(1.-RHOs[2])*A[2]*Fs[2][5]+A[3]*Fs[3][5]+A[4]*Fs[4][5])/A[5]);

	//Solve system of equations
	GRAYDIFFSPEC(iclsd,A,EPS,RHOs,HOs,Fs,id,PIN,POUT);
	//Output
	//Convert emissive powers to temperatures
    	sumq=0;                              //Check total flux=0?
    	for(i=1; i<=N; i=i+1)
		{
        	if(id[i]==0) 
            		T[i]=pow((POUT[i]/sigma),.25);       
        	else
            		q[i]=POUT[i];
        	QA[i]=q[i]*A[i];
        	sumq=sumq+QA[i];
		}    	
    	cout<<"surface     T        q [W/m2]      Q [W]"<<endl;
	for (i=1;i<=N;i=i+1)
		cout<<"\t"<<i<<"\t"<<T[i]<<"\t"<<q[i]<<"\t"<<QA[i]<<endl;
	cout<<"Sum of all fluxes = "<<sumq<<endl;
	  }

//******************************************************************************
// Routine to fill view factor matrix and solve for missing surface temperatures/fluxes
void GRAYDIFFSPEC(int iclsd,double A[],double EPS[],double RHOs[],double HOs[],double Fs[][N+1],int ID[],double PIN[],double POUT[])
{
// INPUT:
// N     = number of surfaces in enclosure
// iclsd = closed or open configuration identifier
//         iclsd=1: configuration is closed; diagonal Fii evaluated from summation rule
//         iclsd/=1: configuration has openings; Fii must be specified
// A(N)  = vector containing surface areas, [m2]
// EPS(N)= vector containing surface emittances
// RHOs(N)= vector containing surface specular reflectance components
// HOs(N) = vector containing external irradiation, in [W/m2]
// Fs(N,N)= vector containing view factors; on input only Fsij with j>i (iclsd=1) or
//         j>=i (iclsd/=1) are required; remainder are calculated
// ID(N) = vector containing surface identifier:
//         ID=0: surface heat flux is specified, in [W/m2]
//         ID=1: surface temperature is specified, in [K]
// PIN(N)= vector containing surface emissive powers (id=1) and fluxes (id=2)
// OUTPUT: 
// POUT(N)= vector containing unknown surface fluxes (for surfaces with id=1) and 
//           temperatures (for surfaces with id=2)
int i,j,ikr;
double qm[N+1][N+1],em[N+1][N+1],pm[N+1][N+1],B[N+1];
// Compute missing view factors
// Lower left triangle by reciprocity
for (i=2;i<=N;i=i+1)
	 for (j=1;j<=(i-1);j=j+1)
		  Fs[i][j]=A[j]/A[i]*Fs[j][i];
// If closed configuration, need to also calculate diagonal terms by summation rule
	if(iclsd==1)
	{
	for (i=1;i<=N;i=i+1)
		{
		Fs[i][i]=1.;
		for (j=1;j<=N;j=j+1)
			if(j!=i)
				Fs[i][i]=Fs[i][i]-(1.-RHOs[j])*Fs[i][j];
		Fs[i][i]=Fs[i][i]/(1.-RHOs[i]);
		}
	}
// Fill q- and e-coefficient matrices
	 for (i=1;i<=N;i=i+1)
		{
		for (j=1;j<=N;j=j+1)
			{
			ikr=(i/j)*(j/i);         // Kronecker delta_ij
			qm[i][j]=ikr/EPS[j]-((1.-RHOs[j])/EPS[j]-1.)*Fs[i][j];
			em[i][j]=ikr-(1.-RHOs[j])*Fs[i][j];
			}
	}
// Fill POUT-coefficient matrix and RHS
	 for (i=1;i<=N;i=i+1)
			{
			B[i]=-HOs[i];
			for (j=1;j<=N;j=j+1)
					{
					pm[i][j]=qm[i][j]*ID[j]-em[i][j]*(1-ID[j]);
					B[i]=B[i]+(em[i][j]*ID[j]-qm[i][j]*(1-ID[j]))*PIN[j];
					}
			}
// Invert POUT-coefficient matrix and multiply by RHS to obtain POUT
	 GAUSS(pm,B,POUT);
}


//******************************************************************************
// View Factor routine
//
double VIEW(int NO,int NARG,double ARG[])
{
//  *************************************************************************
//  *  THIS SUBROUTINE EVALUATES ALL VIEW FACTORS LISTED IN APPENDIX D      *
//  *  NO = CONFIGURATION NUMBER (1 THROUGH 51),                            *
//  *  NARG = NUMBER OF NON-DIMENSIONAL ARGUMENTS FOR VIEW FACTOR,          *
//  *  ARG = ARRAY CONTAINING THE NARG DIFFERENT ARGUMENTS (in alphabetical * 
//  *          order as given in App.D, angles given in degrees),           *
//  *  VIEW = VIEW FACTOR RETURNED BY THE SUBROUTINE                        *
//  *      (FOR NUMBERS 1-9, DF/DL IS GIVEN, WHERE DL IS NONDIMENSIONAL     *
//  *       INFINITESIMAL DIMENSION OF RECEIVING SURFACE).                  *
//  *************************************************************************
      
      
double VIEW,PHI,HR,RTA,RTB,XX,YY,Z,ZZ,LL,PHI1,PHI2,XY,RR1,RR2;
double F,PI,H,W,HW,ALPHA2,R,S,X,R1,R2,RR,SS,C,PHIMAX;
double A,B1,B2,BB1,BB2,B,Y,RTX,RTY,RT,L,HH,WW,W2,H2,HW2;
double H12,W12,C1,C2,C3,AA,BB,BL,L1,L2,D1,D2,OM,CC,D,SR;
 
if(NO<1 || NO>51) goto i;
PI=3.1415926E0;
      
switch (NO)
      {      
  case 1:
      {if(NARG!=1) goto ii;
      PHI=ARG[1];
      VIEW=0.5*cos(PHI*PI/180.);
      return VIEW;}
  case 2:
      {if(NARG!=2) goto ii;
      L=ARG[1];
      R=ARG[2];
      RR=R/L;
      VIEW=2.*RR/pow((1.+RR*RR),2);
      return VIEW;}
  case 3:
      {if(NARG!=2) goto ii;
      R=ARG[1];
      X=ARG[2];
      XX=X/R;
      VIEW=2.*XX/pow((1.+XX*XX),2);
      return VIEW;}
  case 4:
      {if(NARG!=3) goto ii;
      R1=ARG[1];
      R2=ARG[2];
      Z=ARG[3];
      RR=R2/R1;
      if(RR>1.)
	{cout<<"R2 MUST BE LESS THAN R1"<<endl;
	VIEW=0.;
	return VIEW;}
      ZZ=Z/R1;
      X=1.+RR*RR+ZZ*ZZ;
      VIEW=2.*ZZ*(X-2.*RR*RR)*RR/pow((X*X-4.*RR*RR),1.5);
      return VIEW;}
  case 5:
      {if(NARG!=3) goto ii;
      X=ARG[1];
      Y=ARG[2];
      PHI=ARG[3];
      YY=Y/X;
      PHI=PHI*PI/180.;
      VIEW=0.5*YY*pow(sin(PHI),2)/pow((1.+YY*YY-2.*YY*cos(PHI)),1.5);
      return VIEW;}
  case 6:
      {if(NARG!=3) goto ii;
      B=ARG[1];
      R=ARG[2];
      PHI=ARG[3];
      BB=B/R;
      PHI=PHI*PI/180.;
      VIEW=atan(BB)*cos(PHI)/PI;
      return VIEW;}
  case 7:
      {if(NARG!=3) goto ii;
      L=ARG[1];
      R1=ARG[2];
      R2=ARG[3];
      RR=R1/R2;
      BL=L/R1;
      AA=BL*BL+RR*RR+1.;
      VIEW=2.*RR*BL*BL*AA/pow((AA*AA-4.*RR*RR),1.5);
      return VIEW;}
  case 8:
      {if(NARG!=3) goto ii;
      R1=ARG[1];
      R2=ARG[2];
      X=ARG[3];
      RR=R1/R2;
      XX=X/R2;
      AA=XX*XX-RR*RR+1.;
      BB=XX*XX+RR*RR+1.;
      VIEW=2.*XX*AA/pow((BB*BB-4.*RR*RR),1.5);
      return VIEW;}
  case 9:
      {if(NARG!=2) goto ii;
      R=ARG[1];
      X=ARG[2];
      XX=0.5*X/R;
      if(XX>.5E6)
      	{VIEW=0.;
	return VIEW;}
      VIEW=1.-0.5*XX*(2.*XX*XX+3.)/pow((XX*XX+1.),1.5);
      return VIEW;}
  case 10:
      {if(NARG!=3) goto ii;
      A=ARG[1];
      B=ARG[2];
      C=ARG[3];
      AA=A/C;
      BB=B/C;
      RTA=sqrt(1.+AA*AA);
      RTB=sqrt(1.+BB*BB);
      VIEW=(AA/RTA*atan(BB/RTA)+BB/RTB*atan(AA/RTB))/(2.*PI);
      return VIEW;}
  case 11: 
      {if(NARG!=3) goto ii;
      A=ARG[1];
      B=ARG[2];
      C=ARG[3];
      X=A/B;
      Y=C/B;
      RT=sqrt(X*X+Y*Y);
      VIEW=(atan(1./Y)-Y/RT*atan(1./RT))/(2.*PI);
      return VIEW;}
  case 12:
      {if(NARG!=2) goto ii;
      H=ARG[1];
      R=ARG[2];
      HR=H/R;
      VIEW=1./(1.+HR*HR);
      return VIEW;}
  case 13:
      {if(NARG!=3) goto ii;
      A=ARG[1];
      H=ARG[2];
      R=ARG[3];
      HH=H/A;
      RR=R/A;
      ZZ=1.+HH*HH+RR*RR;
      VIEW=.5*(1.-(ZZ-2.*RR*RR)/sqrt(ZZ*ZZ-4.*RR*RR));
      return VIEW;}
  case 14:
      {if(NARG!=3) goto ii;
      H=ARG[1];
      L=ARG[2];
      R=ARG[3];
      HH=H/L;
      RR=R/L;
      ZZ=1.+HH*HH+RR*RR;
      VIEW=.5*HH*(ZZ/sqrt(ZZ*ZZ-4.*RR*RR)-1.);
      return VIEW;}
  case 15:
      {if(NARG!=3) goto ii;
      H=ARG[1];
      L=ARG[2];
      R=ARG[3];
      HH=H/R;
      LL=L/R;
      X=pow((1.+HH),2)+LL*LL;
      Y=pow((1.-HH),2)+LL*LL;
      VIEW=LL/HH*(1./LL*atan(LL/sqrt(HH*HH-1.))
      +(X-2.*HH)/sqrt(X*Y)*atan(sqrt(X/Y*(HH-1.)/(HH+1.))) 
      -atan(sqrt((HH-1.)/(HH+1.))))/PI;
      return VIEW;}
  case 16:
      {if(NARG!=2) goto ii;
      H=ARG[1];
      R=ARG[2];
      VIEW=pow((R/H),2);
      return VIEW;}
  case 17:
      {if(NARG!=2) goto ii;
      H=ARG[1];
      R=ARG[2];
      HH=H/R;
      VIEW=(atan(sqrt(1./(HH*HH-1.)))-sqrt(HH*HH-1.)/HH/HH)/PI;
      return VIEW;}
  case 18:
      {if(NARG!=3) goto ii;
      H=ARG[1];
      R=ARG[2];
      PHI=ARG[3];
      PHI=PHI*PI/180.;
      PHIMAX=acos(R/H);
      if(PHI>PHIMAX)
	{VIEW=0.;
	cout<<"PHI MUST BE LESS THAN ACOS(R/H),"<<180.*PHIMAX/PI<<endl;
	return VIEW;}
      VIEW=pow((R/H),2)*cos(PHI);
      return VIEW;}
  case 19:
      {if(NARG!=3) goto ii;
      H=ARG[1];
      L=ARG[2];
      R=ARG[3];
      HH=H/R;
      LL=L/R;
      if(HH<-1.) 
        VIEW=0.;
      else if(HH<+1.) 
        {X=LL*LL+HH*HH;
        Y=sqrt(X-1.);
        VIEW=(HH/pow(X,1.5)*acos(-HH/LL/Y)-Y*sqrt(1.-HH*HH)/X 
          -asin(Y/LL/LL)+PI/2.)/PI;}
      else
        {X=LL*LL+HH*HH;
        VIEW=pow(H/X,1.5);}
      return VIEW;}
  case 20:
      {if(NARG!=3) goto ii;
      H=ARG[1];
      R=ARG[2];
      Z=ARG[3];
      if(Z<0. || Z>H)
	{VIEW=0.;
	cout<<"INVALID Z"<<endl;
	return VIEW;}
      HH=.5*H/R;
      ZZ=.5*Z/R;
      VIEW=1.+HH-(ZZ*ZZ+0.5)/sqrt(ZZ*ZZ+1.)-(pow((HH-ZZ),2)+0.5)/sqrt(pow((HH-ZZ),2)+1.);
      return VIEW;}
  case 21:
      {if(NARG!=2) goto ii;
      R=ARG[1];
      Z=ARG[2];
      if(Z<0.)
      	{cout<<"INVALID Z"<<endl;
	VIEW=0.;
	return VIEW;}
      ZZ=Z/R;
      if(ZZ>1.E4)
      	{VIEW=pow(ZZ,-3.);
	return VIEW;}
      VIEW=0.5*((ZZ*ZZ+2.)/sqrt(ZZ*ZZ+4.)-ZZ);
      return VIEW;}
  case 22:
      {if(NARG!=3) goto ii;
      R1=ARG[1];
      R2=ARG[2];
      Z=ARG[3];
      RR=R2/R1;
      if(RR>1.)
	{cout<<"R2 MUST BE LESS THAN R1"<<endl;
	VIEW=0.;
	return VIEW;}
      ZZ=Z/R1;
      X=1.+RR*RR+ZZ*ZZ;
      VIEW=0.5*ZZ*(X/sqrt(X*X-4.*RR*RR)-1.);
      return VIEW;}
  case 23:
      {if(NARG!=3) goto ii;
      L=ARG[1];
      X=ARG[2];
      XX=X/L;
      PHI=ARG[3];
      PHI=PHI*PI/180.;
      VIEW=0.5*(1.+(cos(PHI)-XX)/sqrt(1.+XX*XX-2.*XX*cos(PHI)));
      return VIEW;}
  case 24:
      {if(NARG!=2) goto ii;
      PHI1=ARG[1];
      PHI2=ARG[2];
      if(PHI2<PHI1)
	{cout<<"WARNING, PHI2 MUST BE GREATER THAN PHI1"<<endl;
	PHI=PHI1;
	PHI1=PHI2;
	PHI2=PHI;}
      if(PHI1>90.)
	{PHI1=90.;
	cout<<"PHI1 MUST BE LESS THAN 90 OR GREATER THAN -90"<<endl;}
      if(PHI1<-90.)
	{PHI1=-90.;
	cout<<"PHI1 MUST BE LESS THAN 90 OR GREATER THAN -90"<<endl;}
      PHI1=PHI1*PI/180.;
      PHI2=PHI2*PI/180.;
      VIEW=0.5*(sin(PHI2)-sin(PHI1));
      return VIEW;}
  case 25:
      {if(NARG!=3) goto ii;
      A=ARG[1];
      B=ARG[2];
      R=ARG[3];
      AA=A/R;
      BB=B/R;
      VIEW=AA/(AA*AA+BB*BB);
      return VIEW;}
  case 26:
      {if(NARG!=3) goto ii;
      A=ARG[1];
      B=ARG[2];
      C=ARG[3];
      X=A/C;
      Y=B/C;
      XX=sqrt(1.+X*X);
      YY=sqrt(1.+Y*Y);
      VIEW=(YY*atan(X/YY)-atan(X)+X*Y/XX*atan(Y/XX))/(PI*Y);
      return VIEW;}
  case 27:
      {if(NARG!=3) goto ii;
      A=ARG[1];
      B=ARG[2];
      C=ARG[3];
      X=A/B;
      Y=C/B;
      if (Y>1.E6)
      	{VIEW=0.;
	return VIEW;}
      XY=X*X+Y*Y;
      YY=1.+Y*Y;
      VIEW=(atan(1/Y)+.5*Y*log(Y*Y*(XY+1.)/YY/XY)-Y/sqrt(XY)*atan(1./sqrt(XY)))/PI;
      return VIEW;}
  case 28:
      {if(NARG!=4) goto ii;
      H=ARG[1];
      R=ARG[2];
      S=ARG[3];
      X=ARG[4];
      HH=H/R;
      SS=S/R;
      XX=X/R;
      C=SS*SS+XX*XX;
      CC=sqrt(C);
      A=HH*HH+C-1.;
      B=HH*HH-C+1.;
      VIEW=SS/C*(1.-((1./PI)*(acos(B/A)-0.5*(sqrt(A*A+4.*HH*HH)*acos(B/A/CC)-B*asin(1./CC))/HH)) 
           -0.25*A/HH);
      return VIEW;}
  case 29:
      {if(NARG!=1) goto ii;
      PHI=ARG[1];
      VIEW=0.5*(1.+cos(PHI*PI/180.));
      return VIEW;}
  case 30:
      {if(NARG!=3) goto ii;
      A=ARG[1];
      R1=ARG[2];
      R2=ARG[3];
      RR1=R1/A;
      RR2=R2/A;
      if(RR2>1.)
      	{cout<<"R@ MUST BE LESS THAN A"<<endl;
	VIEW=0.;
	return VIEW;}
      VIEW=RR2*RR2/pow((1.+RR1*RR1),1.5);
      return VIEW;}
  case 31:
      {if(NARG!=2) goto ii;
// IDENTICAL TO #21
      R=ARG[1];
      X=ARG[2];
      XX=0.5*X/R;
      if(XX>65000.0)
	{VIEW=0.;
	return VIEW;}
      VIEW=(XX*XX+0.5)/sqrt(1.+XX*XX)-XX;
      return VIEW;}
  case 32:
      {if(NARG!=2) goto ii;
      H=ARG[1];
      W=ARG[2];
      HW=H/W;
      VIEW=sqrt(1.+HW*HW)-HW;
      return VIEW;}
  case 33:
      {if(NARG!=2) goto ii;
      H=ARG[1];
      W=ARG[2];
      HW=H/W;
      VIEW=0.5*(1.+HW-sqrt(1.+HW*HW));
      return VIEW;}
  case 34:
      {if(NARG!=1) goto ii;
      if(ARG[1]>180.)
      	{VIEW=0.;
	return VIEW;}
      ALPHA2=ARG[1]*PI/360.;
      VIEW=1.-sin(ALPHA2);
      return VIEW;}
  case 35:
      {if(NARG!=2) goto ii;
      R=ARG[1];
      S=ARG[2];
      X=1+0.5*S/R;
      VIEW=(asin(1./X)+sqrt(X*X-1.)-X)/PI;
      return VIEW;}
  case 36:
      {if(NARG!=3) goto ii;
      R1=ARG[1];
      R2=ARG[2];
      S=ARG[3];
      RR=R2/R1;
      SS=S/R1;
      C=1.+RR+SS;
      VIEW=(PI+sqrt(C*C-pow((RR+1.),2))-sqrt(C*C-pow((RR-1.),2)) 
             +(RR-1)*acos((RR-1.)/C)-(RR+1)*acos((RR+1.)/C))/(2.*PI);
      return VIEW;}
  case 37:
      {if(NARG!=3) goto ii;
      A=ARG[1];
      B1=ARG[2];
      B2=ARG[3];
      BB1=B1/A;
      BB2=B2/A;
      VIEW=(atan(BB1)-atan(BB2))/(2.*PI);
      return VIEW;}
  case 38:
      {if(NARG!=3) goto ii;
      A=ARG[1];
      B=ARG[2];
      C=ARG[3];
      X=A/C;
      Y=B/C;
      if(X<0.002 || Y<0.002)
      	{VIEW=0.0;
	return VIEW;}
      RTX=sqrt(1.+X*X);
      RTY=sqrt(1.+Y*Y);
      RT=sqrt(1.+X*X+Y*Y);
      VIEW=(log(RTX*RTY/RT)+X*RTY*atan(X/RTY)+Y*RTX*atan(Y/RTX) 
                         -X*atan(X)-Y*atan(Y))*2./(PI*X*Y);
      return VIEW;}
  case 39:
      {if(NARG!=3) goto ii;
      H=ARG[1];
      L=ARG[2];
      W=ARG[3];
      HH=H/L;
      WW=W/L;
      W2=WW*WW;
      H2=HH*HH;
      HW2=H2+W2;
      HW=sqrt(H2+W2);
      H12=H2+1.;
      W12=W2+1.;
      C1=W12*H12/(H12+W2);
      C2=W2*(H12+W2)/W12/HW2;
      C3=H2*(H12+W2)/H12/HW2;
      VIEW=(WW*atan(1./WW)+HH*atan(1./HH)-HW*atan(1./HW) 
         +.25*(log(C1)+W2*log(C2)+H2*log(C3)))/(PI*WW);
      return VIEW;}
  case 40:
      {if(NARG!=3) goto ii;
      A=ARG[1];
      R1=ARG[2];
      R2=ARG[3];
      RR1=R1/A;
      RR2=R2/A;
      X=1.+(1.+RR2*RR2)/(RR1*RR1);
      VIEW=0.5*(X-sqrt(X*X-4.*pow((R2/R1),2)));
      return VIEW;}
  case 41:
      {if(NARG!=3) goto ii;
      L=ARG[1];
      R1=ARG[2];
      R2=ARG[3];
      RR=R1/R2;
      if(RR>1.) 
       {cout<<"R2 MUST BE LESS THAN R1"<<endl;
       VIEW=0.;
       return VIEW;}
      BL=L/R2;
      if(BL>1000.) 
       {VIEW=0.;
       return VIEW;}
      AA=BL*BL+RR*RR-1.;
      BB=BL*BL-RR*RR+1.;
      VIEW=BB/(8.*RR*BL)+(acos(AA/BB) 
        -sqrt(pow(((AA+2.)/RR),2)-4.)/(2.*BL)*acos(AA*RR/BB) 
        -AA/(2.*RR*BL)*asin(RR))/(2.*PI);
      return VIEW;}
  case 42:
      {if(NARG!=2) goto ii;
      H=ARG[1];
      R=ARG[2];
      HH=0.5*H/R;
      VIEW=1.+HH-sqrt(1.+HH*HH);
      return VIEW;}
  case 43:
      {if(NARG!=2) goto ii;
      H=ARG[1];
      R=ARG[2];
      HH=0.5*H/R;
		VIEW=2.*HH*(sqrt(1.+HH*HH)-HH);
      return VIEW;}
  case 44:
      {if(NARG!=3) goto ii;
      H=ARG[1];
      R1=ARG[2];
      R2=ARG[3];
      RR=R2/R1;
     if(RR<1.)
	{cout<<"R2 MUST BE GREATER THAN R1"<<endl;
	VIEW=0.;
	return VIEW;}
      HH=H/R1;
      AA=HH*HH+4.*RR*RR;
      BB=4.*(RR*RR-1.);
      VIEW=1.-1./RR-(sqrt(AA)-HH)/(4.*RR)+(2./RR*atan(sqrt(BB)/HH) 
        -HH/(2.*RR)*(sqrt(AA)/HH*asin((HH*HH*(1.-2./RR/RR)+BB)/(AA-4.)) 
        -asin((RR*RR-2.)/RR/RR)))/PI;
      return VIEW;}
  case 45:
		{if(NARG!=3) goto ii;
      H=ARG[1];
      R1=ARG[2];
      R2=ARG[3];
      RR=R2/R1;
      if(RR<1.)
	{cout<<"R2 MUST BE GREATER THAN R1"<<endl;
        VIEW=0.;
	return VIEW;}
      HH=H/R1;
      if(HH>6250.)
	{VIEW=1.;
	return VIEW;}
      AA=HH*HH+RR*RR-1.;
      BB=HH*HH-RR*RR+1.;
      CC=HH*HH+RR*RR+1.;
      VIEW=1.-AA/(4.*HH)-(acos(BB/AA) 
        -sqrt(CC*CC-4.*RR*RR)/(2.*HH)*acos(BB/RR/AA) 
        -BB/(2.*HH)*asin(1./RR))/PI;
      return VIEW;}
  case 46:
      {if(NARG!=3) goto ii;
      H=ARG[1];
      R1=ARG[2];
      R2=ARG[3];
      RR=R1/R2;
      if(RR>1.)
	{cout<<"R2 MUST BE GREATER THAN R1"<<endl;
	VIEW=0.;
	return VIEW;}
      HH=H/R1;
      AA=1.-RR*RR-HH*HH;
      BB=1.-RR*RR+HH*HH;
      CC=1.+HH*HH+RR*RR;
      X=sqrt(1.-RR*RR);
      Y=RR*AA/BB;
      VIEW=(RR*(atan(X/HH)-atan(2.*X/HH))+0.25*HH*(asin(2.*RR*RR-1.)-asin(RR))
        +0.25*X*X/HH*(0.5*PI+asin(RR))-0.25*sqrt(CC*CC-4.*RR*RR)/HH*(0.5*PI+asin(Y)) 
	+0.25*sqrt(4.+HH*HH)*(0.5*PI+asin(1.-2.*RR*RR*HH*HH/(4.*X*X+HH*HH))))/PI;
      return VIEW;}
  case 47:
      {if(NARG!=3) goto ii;
      D=ARG[1];
      L1=ARG[2];
      L2=ARG[3];
      D1=D/L1;
      D2=D/L2;
      VIEW=0.25*atan(1./sqrt(D1*D1+D2*D2+D1*D1*D2*D2))/PI;
      return VIEW;}
  case 48:
      {if(NARG!=2) goto ii;
      A=ARG[1];
      R=ARG[2];
      RR=R/A;
      VIEW=0.5*(1.-1./sqrt(1.+RR*RR));
      return VIEW;}
  case 49:
      {if(NARG!=2) goto ii;
      A=ARG[1];
      R=ARG[2];
      RR=R/A;
      VIEW=1./sqrt(1.+RR*RR);
      return VIEW;}
  case 50:
      {if(NARG!=4) goto ii;
      R1=ARG[1];
      R2=ARG[2];
      S=ARG[3];
      OM=ARG[4];
      SS=S/R1;
      RR=R2/R1;
      A=asin(1./(SS+1.))*180./PI;
      if(OM<A) cout<<"Formula only valid for OM > "<<A<<"deg"<<endl;
      OM=OM*PI/180.;
      SR=1.+SS+RR/tan(OM);
      VIEW=0.5*(1.-SR/sqrt(SR*SR+RR*RR));
      return VIEW;}
 case 51:
     {if(NARG!=2) goto ii;
      D=ARG[1];
      S=ARG[2];
      SS=D/S;
      if(SS>1)
	{cout<<"D MUST BE LESS THAN S"<<endl;
	VIEW=0.;
	return VIEW;}
      VIEW=SS*acos(SS)+1.-sqrt(1.-SS*SS);
      return VIEW;}
      }
i:    cout<<"***ILLEGAL VALUE FOR NO (NO ="<<NO<<"), MUST BE 1<=NO>=51 ***"<<endl;
VIEW=0.;
return VIEW;
ii:   cout<<"*** WRONG VALUE OF NARG (NARG ="<<NARG<<") FOR NO ="<<NO<<"***"<<endl;
VIEW=0.;
return VIEW;
}



//***********************************************************************************
double perppltf(double X1,double X2,double Y1,double Y2,double Z1,double Z2,double Z3)
{
//  *************************************************************************
// *  THIS SUBROUTINE EVALUATES THE VIEW FACTOR BETWEEN TWO PERPENDICULAR, *
//  *  RECTANGULAR pltf OF SIZE (X2-X1)xZ1 AND (Y2-Y1)x(Z3-Z2), DISPLACED *
//  *  FROM ANOTHER BY Z2 IN THE Z-DIRECTION, BY X1 IN THE X-DIRECTON, AND  *
//  *  BY Y1 IN THE Y-DIRECTION, IN TERMS OF VIEW FACTORS OF DIRECTLY       *
//  *  ADJACENT, IDENTICAL RECTANGLES (CONFIG. 39), AS SHOWN IN FIG.4-15a   *
//  *************************************************************************
    
double perppltf,A,F,ARG[6];
int NARG;
      
      NARG=3;
      F=0.;
      ARG[1]=Y2;
      ARG[2]=Z3;
      ARG[3]=X2;
      A=X2*Z3;
      if(fabs(A*Y2)>1E-6) F=A*VIEW(39,NARG,ARG);
      ARG[1]=Y1;
      A=X2*Z3;
      if(fabs(A*Y1)>1E-6) F=F-A*VIEW(39,NARG,ARG);
      ARG[1]=Y2;
      ARG[3]=X1;
      A=X1*Z3;
      if(fabs(A*Y2)>1E-6) F=F-A*VIEW(39,NARG,ARG);
      ARG[1]=Y1;
      A=X1*Z3;
      if(fabs(A*Y1)>1E-6) F=F+A*VIEW(39,NARG,ARG);
      ARG[1]=Y2;
      ARG[2]=Z2;
      A=X1*Z2;
      if(fabs(A*Y2)>1E-6) F=F+A*VIEW(39,NARG,ARG);
      ARG[1]=Y1;
      A=X1*Z2;
      if(fabs(A*Y1)>1E-6) F=F-A*VIEW(39,NARG,ARG);
      ARG[1]=Y2;
      ARG[3]=X2;
      A=X2*Z2;
      if(fabs(A*Y2)>1E-6) F=F-A*VIEW(39,NARG,ARG);
      ARG[1]=Y1;
      A=X2*Z2;
      if(fabs(A*Y1)>1E-6) F=F+A*VIEW(39,NARG,ARG);
      ARG[1]=Y2;
      ARG[2]=(Z3-Z1);
      A=X2*(Z3-Z1);
      if(fabs(A*Y2)>1E-6) F=F-A*VIEW(39,NARG,ARG);
      ARG[1]=Y1;
      A=X2*(Z3-Z1);
      if(fabs(A*Y1)>1E-6) F=F+A*VIEW(39,NARG,ARG);
      ARG[1]=Y2;
      ARG[3]=X1;
      A=X1*(Z3-Z1);
      if(fabs(A*Y2)>1E-6) F=F+A*VIEW(39,NARG,ARG);
      ARG[1]=Y1;
      A=X1*(Z3-Z1);
      if(fabs(A*Y1)>1E-6) F=F-A*VIEW(39,NARG,ARG);
      ARG[1]=Y2;
      ARG[2]=(Z2-Z1);
      ARG[3]=X2;
      A=X2*(Z2-Z1);
      if(fabs(A*Y2)>1E-6) F=F+A*VIEW(39,NARG,ARG);
      ARG[1]=Y1;
      A=X2*(Z2-Z1);
      if(fabs(A*Y1)>1E-6) F=F-A*VIEW(39,NARG,ARG);
      ARG[1]=Y2;
      ARG[3]=X1;
      A=X1*(Z2-Z1);
      if(fabs(A*Y2)>1E-6) F=F-A*VIEW(39,NARG,ARG);
      ARG[1]=Y1;
      A=X1*(Z2-Z1);
      if(fabs(A*Y1)>1E-6) F=F+A*VIEW(39,NARG,ARG);
      perppltf=F/(2.*(X2-X1)*Z1);
      return perppltf;
}
//***********************************************************************************
double parlpltf(double X1,double X2,double X3,double Y1,double Y2,double Y3,double C)
{
//  *************************************************************************
//  *  THIS SUBROUTINE EVALUATES THE VIEW FACTOR BETWEEN TWO PARALLEL,      *
//  *  RECTANGULAR pltf OF SIZE X1xY1 AND (X3-X2)x(Y3-Y2), DISPLACED      *
//  *  FROM ANOTHER BY C IN THE Z-DIRECTION, BY X2 IN THE X-DIRECTON, AND   *
//  *  BY Y2 IN THE Y-DIRECTION, IN TERMS OF VIEW FACTORS OF DIRECTLY       *
//  *  OPPOSITE, IDENTICAL RECTANGLES (CONFIG. 38), AS SHOWN IN FIG.4-15b   *
//  *************************************************************************
double parlpltf,A,F,ARG[6];
int NARG;

      NARG=3;
      ARG[1]=X3;
      ARG[2]=Y3;
      ARG[3]=C;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=A*VIEW(38,NARG,ARG);
      ARG[2]=Y2;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F-A*VIEW(38,NARG,ARG);
      ARG[2]=Y3-Y1;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F-A*VIEW(38,NARG,ARG);
      ARG[2]=Y2-Y1;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F+A*VIEW(38,NARG,ARG);
      ARG[1]=X2;
      ARG[2]=Y3;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F-A*VIEW(38,NARG,ARG);
      ARG[2]=Y2;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F+A*VIEW(38,NARG,ARG);
      ARG[2]=Y3-Y1;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F+A*VIEW(38,NARG,ARG);
      ARG[2]=Y2-Y1;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F-A*VIEW(38,NARG,ARG);
      ARG[1]=X3-X1;
      ARG[2]=Y3;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F-A*VIEW(38,NARG,ARG);
      ARG[2]=Y2;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F+A*VIEW(38,NARG,ARG);
      ARG[2]=Y3-Y1;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F+A*VIEW(38,NARG,ARG);
      ARG[2]=Y2-Y1;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F-A*VIEW(38,NARG,ARG);
      ARG[1]=X2-X1;
      ARG[2]=Y3;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F+A*VIEW(38,NARG,ARG);
      ARG[2]=Y2;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F-A*VIEW(38,NARG,ARG);
      ARG[2]=Y3-Y1;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F-A*VIEW(38,NARG,ARG);
      ARG[2]=Y2-Y1;
      A=ARG[1]*ARG[2];
      if(fabs(A)>1E-6) F=F+A*VIEW(38,NARG,ARG);
      parlpltf=F/(4.*(X1*Y1));
      return parlpltf;
}

//***********************************************************************************
void GAUSS(double A[][N+1],double B[],double X[])
{
int I,L[N+1],K,J,LK;
double S[N+1];
double SMAX,RMAX,R,XMULT,SUM;
      for (I=1;I<=N;I=I+1)
			{L[I]=I;
			SMAX=0.0;
			for (J=1;J<=N;J=J+1)
				{if (fabs(A[I][J])>SMAX)
				SMAX=fabs(A[I][J]);
				}
			S[I]=SMAX;
			}
      for (K=1;K<=N-1;K=K+1) 
        {RMAX=0.0;
        for (I=K;I<=N;I=I+1) 
		{R = fabs(A[L[I]][K])/S[L[I]];
		if (R<=RMAX)
		{}
		else
		{J=I;
		RMAX=R;
		}}
		LK=L[J];
		L[J]=L[K];
		L[K]=LK;
		for (I=K+1;I<=N;I=I+1)
			{XMULT = A[L[I]][K]/A[LK][K];
				for (J=K+1;J<=N;J=J+1)
					A[L[I]][J] = A[L[I]][J] - XMULT*A[LK][J];
				A[L[I]][K] = XMULT;
			}
          	}    
 
      for (K=1;K<=N-1;K=K+1) 
        { 
        for (I=K+1;I<=N;I=I+1)      
          {B[L[I]] = B[L[I]] - A[L[I]][K]*B[L[K]];
	  //cout<<B[L[I]]<<" "<<L[I]<<endl;
	  } 
        }
      X[N] = B[L[N]]/A[L[N]][N];
      for (I=N-1;I>=1;I=I-1)      
        {SUM = B[L[I]];       
		  for (J=I+1;J<=N;J=J+1)
          SUM = SUM - A[L[I]][J]*X[J];    
        X[I] = SUM/A[L[I]][I];
        }         
}