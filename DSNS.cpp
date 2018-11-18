#include<cmath>
#include<iostream>
#include<fstream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;
typedef struct FCOMPLEX {double r,i;} fcomplex;


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% A program for finding the current in the diffusive wire sandwiched between two 
//% supercondctors for different phase difference of these two superconductors
//% (josephson junction)
//% We consider that the normal region consists of n-1 nodes.
//% Give guessed values of Green's Functions
//% Energy gaps of two superconductors: Delta1=Delta2=1
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main(){

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Data Decleration
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double DvEth;

// Part 1------------------------------------------------------------------------------------
int dd,n,zz,vv;
// errIn: error which the iteration will terminate if the difference of current calculated
// from n nodes and n+1 nodes is smaller than it. 
double errIn=pow(10.0,-3.0);
// t[2]: One dimensional array with size two which contains Temperatures of Superconductors 1&2
// versus their critical temperature. t[1]=T_1/T_c_1, t[2]=T_2/T_c_2
//double t[2];
//--------------------------------------------------------------------------------------------

// Part 2-------------------------------------------------------------------------------------
int i,s;
// Gn; Conductance of Normal wire
// Gn=sigma*A/L; sigma: conductivity of wire; A: cross cestion area; L: Length of the wire
double Gn;
// Conductance: conductance of connector n
// RelConductance: conductance of connector n versus Normal wire conductance
double RelConductance[200],Conductance[200];
//--------------------------------------------------------------------------------------------

// Part 3-------------------------------------------------------------------------------------
// k: kth Matsubara frequency which will increase till the convergence of current in 
// this step is satisfied.
int tt,k,j,q,rr;
double phi,pi=3.14159265358979323846264338327950288419716939937510;
// errIkn: error which the iteration will terminate if the difference of current calculated
// from Matsubara's kth frequency and k+1th frequency is smaller than it.
// AIn: Array of spectral current when we have n nodes. Array relates to phase difference which we have 
// specified.
double errIkn[200],AIn[55][60];
// G: Green's function of each node which will get nearer to real one in each iteration.
		
//---------------------------------------------------------------------------------------------

// Part 4--------------------------------------------------------------------------------------

// f:number of iterations needed for convergence
int f,h,u,v,w,y,ee,ii,jj,ll,qq,ww,yy,qqq,www;
int sigmaz[2][2]={1,0,0,-1};
double  cte,lkg,lkgncte;
// Gdec: Green's function related to leakage current, Gdc2: Gdec to the power of 2
double	Gdecn[2][2],Gdcn2[2][2];
// w(k)=(2.*k+1).*pi
// Gap of energy in superconductor: Delta(T)=1.76.*Tc.*tanh(1.74.*((Tc./T)-1).^1./2)
// TH=tanh(1.74.*((Tc./T)-1).^1./2)
double TH[2],WkvDelta[2],TvEth[2],cte2[2],alpha[50000],tolI[55][200][60],I[55][200][60];
//Imxk1&2: Matrix current for kth frequency of Matsubara when we have n nodes.
// Conductance isn't considered here yet.
//complex <float> ZZ[5];
fcomplex Imxkn1[2][2][200];
fcomplex Imxkn2[2][2][200];
fcomplex Imxkn[2][2][200];
// SgImk= sigmaz*Imxk
fcomplex SgImkn[2][2][200];
// iSgImkn= i*SgImk
fcomplex iSgImkn[2][2][200];
// TriSgImkn: Trace of Matrix iSgImk
fcomplex TriSgImkn[200];
// Ikn:This is spectral current for Matsubara's kth frequency when we have n nodes.
fcomplex Ikn[200];
// hIn: half of current for n node. half, becouse positive frequencies are considered.
fcomplex hIn[200];
// In: current related to n nodes, it's hIn times 2. Becouse current related to negative
// frequencies are equal to current related to positive frequencies.
fcomplex In[200];
// tolIn.r=Ikn[n-1].r/In[n-1].r; 
fcomplex tolIn[200];
//---------------------------------------------------------------------------------------------
			
// Part 5--------------------------------------------------------------------------------------
int d,kk,kkk;
double one[2][2]={1.0,0.0,0.0,1.0};
double max[210];
fcomplex Gkn[2][2][200];
fcomplex Gnewkn[2][2][200];
fcomplex Mkn[2][2][200];
fcomplex ANCOMMsmspkn1[2][2][200];
fcomplex ANCOMMsmspkn2[2][2][200];
fcomplex ANCOMMsmspkn[2][2][200];
fcomplex ANCOMMsmdeckn[2][2][200];
fcomplex ANCOMMdecspkn[2][2][200];
fcomplex Dkn[2][2][200];
fcomplex tolGkn[2][2][200];
fcomplex errGkn[2][2][200];
//----------------------------------------------------------------------------------------------


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
//printf("input the ratio of critical temperature to Thouless energy  ");
DvEth=20;

//Open output files
//|||||||||||||||||||||||||||||||||||
ofstream out_stream;out_stream.open("2DvEth20p2.dat");
ofstream out_stream9;out_stream9.open("tolI49DvEth20p2.dat");
ofstream out_stream10;out_stream10.open("2maximumDvEth20p2.dat");
//||||||||||||||||||||||||||||||||||||
out_stream<<"DvEth, the ratio of gap to Thouless energy is:"<<endl;
out_stream<<DvEth<<endl;
out_stream<<"Connector 0: connects S1 to Node'1' in Normal wire"<<endl;
out_stream<<"Connector 1: connects Node'1' in Normal to Node'2' in Normal wire"<<endl;
out_stream<<"Connector 2: connects Node'2'in Normal to S2"<<endl;

//********************************************************************************************
// PART 1
//A loop which goes over different values of t.
for (dd=1; dd<=5; dd++){
    
	n=60;
	TvEth[0]=0.1*dd;
	TvEth[1]=0.1*dd;
	Gn=1.0;
	zz=1;

	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	// PART 2
	// This part will find the number of nodes which converges our results \\
	
	// This loop will continue till the n which satisfies the convergence of current is found.
	while (zz==1){
		n=n+1;
        Conductance[0]=n*Gn;
		Conductance[n-1]=n*Gn;
		// Here, number of nodes in Normal wire is n-1
		for (i=1; i<=n-2; i++){
			Conductance[i]=((n-1)+1)*Gn;
		}
		for (i=0; i<=n-1; i++){
			RelConductance[i]=Conductance[i]/Gn;
		}
      
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		// PART 3
		// This part will find the current for each specified phase difference, phi.\\
		
		
		for (tt=0; tt<=50; tt++){
			phi=(tt*(pi))/50;
			k=0;
			// initial guess for Green's functions of nodes in Normal wire
			for (q=1; q<=n-1; q++){
				Gkn[0][0][q].r=1.0;Gkn[0][0][q].i=0.0;Gkn[0][1][q].r=0.0;Gkn[0][1][q].i=0.0;Gkn[1][0][q].r=0.0;
				Gkn[1][0][q].i=0.0;Gkn[1][1][q].r=-1.0;Gkn[1][1][q].i=0.0;
			}
	  		for (j=0; j<=n-1; j++){
				tolIn[j].r=0.0;tolIn[j].i=0.0;errIkn[j]=pow(10.0,-6.0);hIn[j].r=0.0;hIn[j].i=0.0;
          		}
		
             		rr=1;
			//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
			//PART 4
			// This part will iterate till the convergence of current related to Matsubara's frequency
			// for specified phi is satisfied
			
			while  (rr==1){
				f=0;
				k=k+1;
	   			
				alpha[k]=(2*(k)-1)*pi;
				for (h=0; h<=1; h++){
					//TH[h]=(tanh(1.74*(sqrt((1/t[h])-1))));
					cte2[h]=((pow(alpha[k],2))*(pow(TvEth[h],2))+pow(DvEth,2));
				}
				cte=(1/(sqrt(cte2[0])));
	    			fcomplex a;
				a.r=0.0;a.i=-1.0;
				fcomplex b;
				if (phi==pi){
					b.i=0.0;b.r=-1.0;
				}
				else if (phi==pi/2){
					b.i=1.0;b.r=0.0;
				}else{
				b.r=cos(phi);b.i=sin(phi);
				}
				fcomplex c;
				c.r=a.r*b.r-a.i*b.i;
				c.i=a.i*b.r+a.r*b.i;
				fcomplex cc;
				cc.r=a.r*b.r-a.i*b.i;
				cc.i=-a.i*b.r-a.r*b.i;

			
				// Green's functions of superconductors:
				Gkn[0][0][n].r=cte*(alpha[k])*(TvEth[0]); Gkn[0][0][n].i=0.0;
				Gkn[0][1][n].r=cte*(DvEth)*(c.r); Gkn[0][1][n].i=cte*(DvEth)*(c.i);
				Gkn[1][0][n].r=cte*(DvEth)*cc.r; Gkn[1][0][n].i=cte*(DvEth)*cc.i;
				Gkn[1][1][n].r=-cte*(alpha[k])*(TvEth[0]); Gkn[1][1][n].i=0.0;
				Gkn[0][0][0].r=cte*(alpha[k])*(TvEth[0]); Gkn[0][0][0].i=0.0;
				Gkn[0][1][0].r=0.0;Gkn[0][1][0].i=cte*(DvEth)*a.i;
				out_stream<<Gkn[0][1][0].r<<endl;
				Gkn[1][0][0].r=0.0;Gkn[1][0][0].i=-cte*(DvEth)*a.i;
				Gkn[1][1][0].r=-cte*(alpha[k])*(TvEth[0]);Gkn[1][1][0].i=0.0;
                
                for (qqq=0; qqq<=1; qqq++){
					for (qq=0; qq<=1; qq++){
						for (q=1; q<=n-1; q++){
							errGkn[qqq][qq][q].r=pow(10.0,-10.0);
							errGkn[qqq][qq][q].i=pow(10.0,-10.0);
						}
					}
               	}
								 
	
				// leakage
				lkgncte=-(2.0/(n))*pi*((2.0)*k-1.0);
				lkg=lkgncte*TvEth[0];
				double Gdecn[2][2]={-lkg,0.0,0.0,lkg};//lkgcte*t[0]*(-TcvEth);Gdecn[0][1]=0.0;Gdecn[1][0]=0.0;Gdecn[1][1]=lkgcte*t[0]*TcvEth;
                                

				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				// PART 5
				// This part will calculate Green's function for each node with kth Matsubara's frequency by iteration
			
				d=1;
				while (d==1){
					for (kk=0; kk<=1; kk++){
						for (j=0; j<=1; j++){
							for (s=1; s<=n-1; s++){
								Mkn[kk][j][s].r=RelConductance[s-1]*(Gkn[kk][j][s-1].r)+RelConductance[s]*(Gkn[kk][j][s+1].r)+Gdecn[kk][j];
								Mkn[kk][j][s].i=RelConductance[s-1]*(Gkn[kk][j][s-1].i)+RelConductance[s]*(Gkn[kk][j][s+1].i);
						
								ANCOMMsmspkn1[kk][j][s].r=0.0;ANCOMMsmspkn1[kk][j][s].i=0.0;
								ANCOMMsmspkn2[kk][j][s].r=0.0;ANCOMMsmspkn2[kk][j][s].i=0.0;
								ANCOMMsmdeckn[kk][j][s].r=0.0;ANCOMMsmdeckn[kk][j][s].i=0.0;
								ANCOMMdecspkn[kk][j][s].r=0.0;ANCOMMdecspkn[kk][j][s].i=0.0;

								for (i=0; i<=1; i++){
									ANCOMMsmspkn1[kk][j][s].r=(Gkn[kk][i][s+1].r)*(Gkn[i][j][s-1].r)-(Gkn[kk][i][s+1].i)*(Gkn[i][j][s-1].i)+ANCOMMsmspkn1[kk][j][s].r;
									ANCOMMsmspkn1[kk][j][s].i=(Gkn[kk][i][s+1].i)*(Gkn[i][j][s-1].r)+(Gkn[kk][i][s+1].r)*(Gkn[i][j][s-1].i)+ANCOMMsmspkn1[kk][j][s].i;
									ANCOMMsmspkn2[kk][j][s].r=(Gkn[kk][i][s-1].r)*(Gkn[i][j][s+1].r)-(Gkn[kk][i][s-1].i)*(Gkn[i][j][s+1].i)+ANCOMMsmspkn2[kk][j][s].r;
									ANCOMMsmspkn2[kk][j][s].i=(Gkn[kk][i][s-1].i)*(Gkn[i][j][s+1].r)+(Gkn[kk][i][s-1].r)*(Gkn[i][j][s+1].i)+ANCOMMsmspkn2[kk][j][s].i;
									ANCOMMsmspkn[kk][j][s].r=ANCOMMsmspkn2[kk][j][s].r+ANCOMMsmspkn1[kk][j][s].r;
									ANCOMMsmspkn[kk][j][s].i=ANCOMMsmspkn2[kk][j][s].i+ANCOMMsmspkn1[kk][j][s].i;

									ANCOMMsmdeckn[kk][j][s].r=(Gkn[kk][i][s-1].r)*Gdecn[i][j]+Gdecn[kk][i]*(Gkn[i][j][s-1].r)+ANCOMMsmdeckn[kk][j][s].r;
									ANCOMMsmdeckn[kk][j][s].i=(Gkn[kk][i][s-1].i)*Gdecn[i][j]+Gdecn[kk][i]*(Gkn[i][j][s-1].i)+ANCOMMsmdeckn[kk][j][s].i;
			
									//c.r=a.r*b.r-a.i*b.i;
									//c.i=a.i*b.r+a.r*b.i;
							
									ANCOMMdecspkn[kk][j][s].r=(Gkn[kk][i][s+1].r)*Gdecn[i][j]+Gdecn[kk][i]*(Gkn[i][j][s+1].r)+ANCOMMdecspkn[kk][j][s].r;
									ANCOMMdecspkn[kk][j][s].i=(Gkn[kk][i][s+1].i)*Gdecn[i][j]+Gdecn[kk][i]*(Gkn[i][j][s+1].i)+ANCOMMdecspkn[kk][j][s].i;
								}
							}
						}
					}
	
     
					for (i=0; i<=1;i++){
						for (j=0; j<=1; j++){
							for (s=1; s<=n-1; s++){
								Gdcn2[i][j]=0.0;
								for (kkk=0; kkk<=1; kkk++){
									Gdcn2[i][j]=Gdecn[i][kkk]*Gdecn[kkk][j]+Gdcn2[i][j];
									Dkn[i][j][s].r=(pow(RelConductance[s-1],2.0))*one[i][j]+(pow(RelConductance[s],2.0))*one[i][j]+Gdcn2[i][j]+
										RelConductance[s]*RelConductance[s-1]*(ANCOMMsmspkn[i][j][s].r)+(RelConductance[s-1])*(ANCOMMsmdeckn[i]
										[j][s].r)+(RelConductance[s])*(ANCOMMdecspkn[i][j][s].r);
                   					Dkn[i][j][s].i=RelConductance[s]*RelConductance[s-1]*(ANCOMMsmspkn[i][j][s].i)+(RelConductance[s-1])*
										(ANCOMMsmdeckn[i][j][s].i)+(RelConductance[s])*(ANCOMMdecspkn[i][j][s].i);
								}
						
								Gnewkn[i][j][s].r=(Mkn[i][j][s].r)/(sqrt(Dkn[0][0][s].r));
								Gnewkn[i][j][s].i=(Mkn[i][j][s].i)/(sqrt(Dkn[0][0][s].r));
                        
								tolGkn[i][j][s].r=Gnewkn[i][j][s].r-Gkn[i][j][s].r;
								tolGkn[i][j][s].i=Gnewkn[i][j][s].i-Gkn[i][j][s].i;
								tolGkn[i][j][s].r=abs(tolGkn[i][j][s].r);
								tolGkn[i][j][s].i=abs(tolGkn[i][j][s].i);
								}
						}
					}
	
					
					f=f+1;
			
					// Here, improved Green's function is replaced with the old one
					for (i=0; i<=1; i++){
						for (j=0; j<=1; j++){
							for (s=1; s<=n-1; s++){
								Gkn[i][j][s].r=Gnewkn[i][j][s].r;
								Gkn[i][j][s].i=Gnewkn[i][j][s].i;
							}
						}

					}

					d=0;
					
					for (i=0; i<=1; i++){
						for (j=0; j<=1; j++){ 
							for (s=1; s<=n-1; s++) {
								if (tolGkn[i][j][s].r>errGkn[i][j][s].r || tolGkn[i][j][s].i>errGkn[i][j][s].i){
									d=1;
									break;
								}
								
							}
						}
						
					}
				}

				//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

				for (i=0; i<=1; i++){
					for (j=0; j<=1; j++){
						for (v=0; v<=n-1; v++){
							Imxkn1[i][j][v].r=0.0;
							Imxkn1[i][j][v].i=0.0;
							Imxkn2[i][j][v].r=0.0;
							Imxkn2[i][j][v].i=0.0;
		  
							for (ll=0; ll<=1; ll++){
		  
   								Imxkn1[i][j][v].r=(Gkn[i][ll][v+1].r)*(Gkn[ll][j][v].r)-(Gkn[i][ll][v+1].i)*(Gkn[ll][j][v].i)+Imxkn1[i][j][v].r;
								Imxkn1[i][j][v].i=(Gkn[i][ll][v+1].i)*(Gkn[ll][j][v].r)+(Gkn[i][ll][v+1].r)*(Gkn[ll][j][v].i)+Imxkn1[i][j][v].i;
								Imxkn2[i][j][v].r=(Gkn[i][ll][v].r)*(Gkn[ll][j][v+1].r)-(Gkn[i][ll][v].i)*(Gkn[ll][j][v+1].i)+Imxkn2[i][j][v].r;
								Imxkn2[i][j][v].i=(Gkn[i][ll][v].i)*(Gkn[ll][j][v+1].r)+(Gkn[i][ll][v].r)*(Gkn[ll][j][v+1].i)+Imxkn2[i][j][v].i;
								Imxkn[i][j][v].r=Imxkn2[i][j][v].r-Imxkn1[i][j][v].r;
								Imxkn[i][j][v].i=Imxkn2[i][j][v].i-Imxkn1[i][j][v].i;
							}
						}
					}
				}

 
				for (w=0; w<=1; w++){
					for (ww=0; ww<=1; ww++){
						for (v=0; v<=n-1; v++){
							SgImkn[w][ww][v].r=0;
							SgImkn[w][ww][v].i=0;
							for (www=0; www<=1; www++){
								SgImkn[w][ww][v].r=sigmaz[w][www]*(Imxkn[www][ww][v].r)+ (SgImkn[w][ww][v].r);
								SgImkn[w][ww][v].i=sigmaz[w][www]*(Imxkn[www][ww][v].i)+ (SgImkn[w][ww][v].i);
							}
						}
					}
				}


				for (w=0; w<=1; w++){
					for (ww=0; ww<=1; ww++){
						for (v=0; v<=n-1; v++){
							iSgImkn[w][ww][v].i=SgImkn[w][ww][v].r;
							iSgImkn[w][ww][v].r=-SgImkn[w][ww][v].i;
						}
					}
				}

				for (v=0; v<=n-1; v++){
					TriSgImkn[v].r=0.0;
					TriSgImkn[v].i=0.0;
				}

				for (y=0; y<=1; y++){
					for (yy=0; yy<=1; yy++){
						for (v=0; v<=n-1; v++){
							if (y==yy){
								TriSgImkn[v].r=iSgImkn[y][yy][v].r+TriSgImkn[v].r;
								TriSgImkn[v].i=iSgImkn[y][yy][v].i+TriSgImkn[v].i;
							}
						}
					}
				}
		
				for (v=0; v<=n-1; v++){

					Ikn[v].r=(pi/4)*(TvEth[0])*RelConductance[v]*(TriSgImkn[v].r);
					Ikn[v].i=(pi/4)*(TvEth[0])*RelConductance[v]*(TriSgImkn[v].i);
					hIn[v].r=Ikn[v].r+hIn[v].r;
 					hIn[v].i=Ikn[v].i+hIn[v].i;
					In[v].r=(2.0)*hIn[v].r;
					In[v].i=(2.0)*hIn[v].i;
			
					if ((In[v].r)!=0.0){
						tolIn[v].r=(Ikn[v].r)/(In[v].r);
					}
				}
  				
				rr=0;	   
				for (ee=0; ee<=n-1; ee++){
					if (In[ee].r==0.0){
						rr=0;
					}
					else if (tolIn[ee].r>errIkn[ee]){
						rr=1;
						break;
					}
				}
	 				
	 			
			}

			//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

			AIn[dd][tt]=In[0].r;
		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		for (vv=0; vv<=50; vv++){
			I[dd][n][vv]=AIn[dd][vv];
			//out_stream<<"test2"<<endl;
			out_stream<<"I["<<dd<<"]["<<n<<"]"<<"["<<vv<<"]"<<endl;
			out_stream<<I[dd][n][vv]<<endl;
			//out_stream<<"\n"<<endl;
		}
	
     	 
	        for (ii=0; ii<=50; ii++){
			if (I[dd][n][ii]!=0.0){
				tolI[dd][n][ii]=(I[dd][n][ii]-I[dd][n-1][ii])/I[dd][n][ii];
			 	tolI[dd][n][ii]=abs(tolI[dd][n][ii]);
            }
			else{
				tolI[dd][n][ii]=0.0;
			}		
			 	
		}
		 
	 	out_stream9<<tolI[dd][n][49]<<endl;
		
		zz=0;
		for (s=1; s<=49; s++){
			if (tolI[dd][n][s]>errIn){
				zz=1;
			}
			break;
		}
		
	}

    	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    	max[dd]=AIn[dd][0];
	for (vv=0; vv<=50; vv++){
		 out_stream<<"converged current with respect to n"<<endl;
		 out_stream<<"I("<<vv*(pi/50.0)<<")="<<endl;
		 out_stream<<"AIn["<<dd<<"]["<<vv<<"]"<<endl;
		 out_stream<<AIn[dd][vv]<<endl;
		 if (AIn[dd][vv]>max[dd]){
		 	max[dd]=AIn[dd][vv];
		 }
	}
	out_stream10<<max[dd]<<endl;
	out_stream<<"number of nodes needed for convergence for TvEth="<<TvEth[1]<<"is:"<<endl;
	out_stream<<n<<endl;
	out_stream<<"\n"<<endl;
}

//********************************************************************************************
  
out_stream.close();
out_stream9.close();
out_stream10.close(); 

}
