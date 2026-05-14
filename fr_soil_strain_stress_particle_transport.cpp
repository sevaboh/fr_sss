// H equation + elasticity equation + transport for sand/silt/clay particles in soil
// vgm through modelled sand/silt/clay using Python Rosetta module
// porousity equal to theta_s from vgm
// strain-stress parameters of soil are fixed
// density of soil particles depends on sand/silt/clay through empirical relationship

#include "gamma.h"
#include "basic_solver.h"
#include "particle_transport.h"
//////////////////////////////////////
///// global variables ///////////////
//////////////////////////////////////
const double HtoPmult=9.8;

// strain-stress parameters
double rho_w=1000; // density of water
double lambda=19788461.54;
double mu=13192307.69;
double B=0.4,B2=0.4;
double sigma=0.6; // starting porosity
const double g=9.81;
// particle transport parameters for sand/silt/clay
double pt_l_a[3]={0.001,0.001*2,0.001/2};
double pt_l_s[3]={0.1,0.1,0.2};
double pt_k_det[3]={0.0001/86400,0.0002/86400,0.0001/86400};
double pt_c0[3]={0.25,0.35,0.4};
double ls_eps3=1e-8;
// fractional orders
double alpha=1;
double beta=1;
// for modelling irrigation
double irr_depth=0.3; // control depth, m
double irr_p0=-2; // lower limit of average water head,m
double irr_p1=-0.5; // upper limit of average water head,m
double irr_flow=1e-6; // surface irrigation flow, m/s
double irr_tau=60; // tau when irrigating
int is_irr=0;
double irr_time0=120*86400,irr_time1=270*86400; // range of time in which to apply irrigation, s
////////// auxiliary
// space-fractional derivative non-local part for (ii,x)
double SFD(int ii,double *U,double gamma,double *&sp_precalc,int N)
{
	double nlp=0.0;
	// precalc coefficients
#pragma omp critical
	if (gamma!=1)
	if (sp_precalc==NULL)
	{
		if (sp_precalc) delete [] sp_precalc;
		sp_precalc = new double[N + 2];
		for (int i = 0;i < N + 2;i++)
			sp_precalc[i] = pow((double)i + 1, 1.0-gamma) - pow((double)i , 1-gamma);
	}
	if (gamma!=1)
	    for (int i = 0;i <ii;i++)
		    nlp+=sp_precalc[ii-i]*U[i];
	return nlp;
}

class basic_fr_solver: public basic_solver {
public:
	double Dg;
	double alpha;
	double *sp_precalc=NULL;
	double balance;
	basic_fr_solver(double _L,double _N,double e):basic_solver(_L,_N,e)
	{
	    balance=1e30;
	}
};
//////////////////////////////////////
// Solver for Richards equation
//////////////////////////////////////
class H_solver: public basic_fr_solver {
public:
	double *b_V; // water movement velocity
	double *b_Vd; // d velocity / dz
	double *pr_D2=NULL;
	double *pr_D3=NULL;
	double *pr_zd=NULL;
	// auxiliary
	double B_;
	basic_fr_solver *Zs;
	particle_transport *sand,*silt,*clay;
	double H0_; // initial H parameter
	FILE *v_out;
	// inputs
	// EV per day
	double *EV_T; // times
	double *EV_F; // values
	int nEV_T; // size of EV_T and EV_F
	double pr_ET; // current evapotranspiration
	// precipitation
	double *perc_T, *perc_A; // times and values
	int nperc; // size
	double pr_perc; // current precipitation
	// grownwater level
	double *flT, *flV; // times and values
	int nfl; // size of flT,flV
	double pr_fl; // precalculated growndwater level

	// values precalculated on the start of time step
	double *pr_dwdh=NULL; // d theta / dh
	double *pr_w=NULL; // wetness
	double *pr_K=NULL; // hydraulic conductivity

	// modelling irrigation
	int irr_on=0;
	double pr_ah;
	/// calc wetness on the base of van Genuchten model ////////
	double wetness(int i,double P)
	{
		P *= HtoPmult;
		if (P <= 0.0)
			return pr_vgm_s1s[i];
		return pr_vgm_s0s[i]+((pr_vgm_s1s[i] - pr_vgm_s0s[i]) / pow(1 + pow(pr_vgm_as[i]*P*10.19716, pr_vgm_ns[i]), (1 - 1 / pr_vgm_ns[i])));
	}
	// precalulculate wetness in all cells in pr_w
	void pr_wetness()
	{
	    // alloc
#pragma omp critical
	    if (pr_w==NULL)
		pr_w=new double[N+2];
#pragma omp barrier
#pragma omp for
	    for (int i=0;i<=N;i++)
		pr_w[i]=wetness(i,-b_U[i]);
#pragma omp barrier
	}
	// precalculate d theta / dh in all cells in pr_dwdh
	double dw_dh(double P,int i)
	{
		double Ch;
		if (P <= 0.0)
			Ch = 0.0;
		else
		{
			P*=10.19716;
			P*=HtoPmult;
			double aPn=pow(pr_vgm_as[i]*P, pr_vgm_ns[i]);
			Ch=-(1.0/(P/(HtoPmult*10.197196)))*(((1.0/pr_vgm_ns[i])-1)*pr_vgm_ns[i]*aPn*pow(aPn+1.0,(1.0/pr_vgm_ns[i])-2.0)*(pr_vgm_s1s[i]-pr_vgm_s0s[i]));
		}
		return Ch;
	}
	void pr_dw_dh()
	{
	    // alloc
#pragma omp critical
	    if (pr_dwdh==NULL)
		pr_dwdh=new double[N+2];
#pragma omp barrier
#pragma omp for
	    for (int i=0;i<=N;i++)
	    {
		double P;
		double Ch = 0.0;
		int idx=0;
		P = -b_U[idx=i];
		// calculate
		Ch = dw_dh(P,i);
		if (!isfinite(Ch)) Ch=0;
		pr_dwdh[idx]=Ch;
	    }
#pragma omp barrier
	}
	double avg_h()
	{
	    double ah=0;
	    int nh=0;
	    for (int i=0;i<=N;i++)
		if (i*dL<irr_depth)
		{
		    ah+=b_U[i];
		    nh++;
		}
	    if (nh) ah/=nh;
	    return ah;
	}
	double Kr(double w,int i)
	{
		// Mualem's model
		double kr = pow(((w - pr_vgm_s0s[i]) / (pr_vgm_s1s[i] - pr_vgm_s0s[i])), -pr_vgm_power[i])
		    *pow(1.0 - pow(1.0 - pow(((w - pr_vgm_s0s[i]) / (pr_vgm_s1s[i] - pr_vgm_s0s[i])), 
			    (1.0 / (1.0 - 1.0 / pr_vgm_ns[i]))), (1.0 - 1.0 / pr_vgm_ns[i])), 2.0);
		return kr;
	}
	// precalculate hydraulic conductivity in all cells in pr_K
	void pr_KK()
	{
	    // alloc
#pragma omp critical
	    if (pr_K==NULL)
		pr_K=new double[N+2];
#pragma omp barrier
#pragma omp for
	    for (int i=0;i<=N;i++)
		pr_K[i]=Kr(pr_w[i],i)*pr_vgm_k[i];
#pragma omp barrier
	}
	// sets current ET, precipitation, root system depth, groundwater level
	void precalc_values()
	{
		double v=0.0,k;
		int i=0;
		// precipitation
		pr_perc=linear_from_pair(perc_T,perc_A,nperc);
		// Evapotranspiration
		pr_ET=linear_from_pair(EV_T,EV_F,nEV_T);
		// grownwater level
		pr_fl=linear_from_pair(flT,flV,nfl);
		// head in control layer
		pr_ah=avg_h();

		if (pr_fl>0) pr_perc=pr_ET=0; // no flow if soil is submerged under water
	}
	// upper flow
	double flow0(double *pickard_prevU)
	{
	    double flow=pr_ET-pr_perc;
	    // irrigation
	    if (is_irr)
		if (irr_on==1)
		    flow-=irr_flow;
	    return flow;
	}
	// upper boundary condition (DbU=Uc)
	double Uc(double *pickard_prevU)
	{
	    double flow=flow0(pickard_prevU);
	    // condition for flux
	    double kk;
	    kk=(pr_K[1]+Km1(0,pickard_prevU))/2;
	    double ret=dL*((flow/kk)+1);
	    return ret;
	}
	// coefficients of three-diagonal linear equations system
	double Km1(int ii,double *pickard_prevU) // K in point ii-1
	{
		double km1;
		if (ii!=0) // G(2-a)*z^(a-1)*D^a f | z->0 -> f'(0)
		    km1=pr_K[ii-1];
		else // in -1 through linear change of h
		    km1=pr_vgm_k[0]*Kr(wetness(0,-(2*pickard_prevU[0]-pickard_prevU[1])),0);
		return km1;
	}
	double B(int i,double *pickard_prevU) // for H_i+1
	{
		if (i==0) (1.0/(dL*dL))*(-(0.5*pr_K[i+1]-0.5*Km1(i,pickard_prevU))-pr_K[i]);
		return Dg*pow(dL,1-alpha)*pow((i+1)*dL,alpha-1)*pow((i+1)*dL,alpha-1)*(-(0.5*pr_K[i+1]-0.5*Km1(i,pickard_prevU))-pr_K[i]);
	}
	double A_(int i,double *pickard_prevU) // for H_i-1
	{
		if (i==0) return -pr_K[i]/(dL*dL);
		return -Dg*pow(dL,1-alpha)*pow((i+1)*dL,alpha-1)*pow(i*dL,alpha-1)*pr_K[i];
	}
	double R(int i,double *pickard_prevU) // for H_i
	{
		if (i==0)
			return -(1.0/(dL*dL))*(-0.5*pr_K[i+1]+0.5*Km1(i,pickard_prevU))+pr_K[i]*2.0/(dL*dL)+
				((steady_state==0)?((pr_dwdh[i]+pr_w[i]*pr_vgm_cw[i])/tau):0);
		return -Dg*pow(dL,1-alpha)*pow((i+1)*dL,alpha-1)*pow((i+1)*dL,alpha-1)*(-0.5*pr_K[i+1]+0.5*Km1(i,pickard_prevU))+
			Dg*pow(dL,1-alpha)*pow((i+1)*dL,alpha-1)*pr_K[i]*(pow((i+1)*dL,alpha-1)+pow(i*dL,alpha-1))+
			((steady_state==0)?((pr_dwdh[i]+pr_w[i]*pr_vgm_cw[i])/tau):0);
	}
	// right part 
	double Rp(int ii,double *b_Uold,double *pickard_prevU)
	{
		double ret=0.0;
		ret = ((steady_state==0)?(((pr_dwdh[ii]+pr_w[ii]*pr_vgm_cw[ii])/tau)*b_Uold[ii]):0);
		if (ii==0)
			ret-=(pr_K[ii+1]-Km1(ii,pickard_prevU))/(2*dL);
		else
			ret-=Dg*pow(dL,1-alpha)*pow((ii+1)*dL,alpha-1)*dL*dL*pow((ii+1)*dL,alpha-1)* (pr_K[ii+1]-Km1(ii,pickard_prevU))/(2*dL);
		if (ii!=0)
			ret+=Dg*pow(dL,1-alpha)*pow((ii+1)*dL,alpha-1)*dL*pr_K[ii]*(pow(ii*dL,alpha-1)-pow((ii+1)*dL,alpha-1));
		return ret;
	}
	// linear solver
	// matrix multiplication of vector UU
	void mmult(double *UU)
	{
		mmult_main(UU);
#pragma omp single
		{
			// bottom boundary condition - h level
			MM[N]=UU[N];
			// upper boundary condition - flux
			MM[0]=UU[0]+((pr_A[0]+pr_B[0])/pr_R[0])*UU[1]; // second order condition through "ghost point"
		}
#pragma omp barrier
	}
	// precalculations for current time step
	void __precalc(double *pickard_prevU)
	{
#pragma omp single
		{ // precalculate non-linear coefficients using the solution of previous iteration
		    memcpy(sb_U,b_U,sizeof(double)*(N+2));
		    memcpy(b_U,pickard_prevU,sizeof(double)*(N+2));
		}
#pragma omp barrier
		pr_wetness();
		pr_dw_dh();
		pr_KK();
		// precalc matrix coefficients
#pragma omp single
	    {
		precalc_values();
		if (pr_A==NULL) pr_A=new double[N+2];
		if (pr_B==NULL) pr_B=new double[N+2];
		if (pr_R==NULL) pr_R=new double[N+2];
		if (pr_D2==NULL) pr_D2=new double[N+2];
		if (pr_D3==NULL) pr_D3=new double[N+2];
		if (pr_zd==NULL) pr_zd=new double[N+2];
		if (steady_state)
		    calculate_velocity();
		memcpy(b_U,sb_U,sizeof(double)*(N+2));
		// irrigation
		if (is_irr)
		{
		    if (irr_on==1)
		    {
		     if ((pr_ah>irr_p1)||(T>irr_time1))
			irr_on=0;
		    }
		    else
			if ((pr_ah<irr_p0)&&((T>irr_time0)&&(T<irr_time1)))
			    irr_on=1;
		}
		if (irr_on)
		if (tau>irr_tau)
		{
			T-=tau;
			tau=irr_tau;
			T+=tau;
		}
		// recalc vgm
		particle_transport *ssc[]={sand,silt,clay};
		python_get_mean_vgm(ssc,N,(T-tau)/86400);
		for (int i=0;i<N;i++)
		{
		    pr_vgm_s0s[i]=vgm_mean[5*i+0];
		    pr_vgm_s1s[i]=vgm_mean[5*i+1];
		    pr_vgm_as[i]=vgm_mean[5*i+2];
		    pr_vgm_ns[i]=vgm_mean[5*i+3];
		    pr_vgm_k[i]=vgm_mean[5*i+4];
		}
	    }
		// for D^a H
#pragma omp barrier
		if (alpha!=1)
#pragma omp for nowait
		for (int i = 0;i <= N ;i++)
		{
		    if ((i!=0)&&(i!=N))
			pr_D2[i]= ((pickard_prevU[i+1] - pickard_prevU[i])/dL)-1;
		    pr_D2[0]= flow0(pickard_prevU)/(0.5*(pr_K[0]+pr_K[1]));
		    pr_D2[N]= ((pickard_prevU[N] - pickard_prevU[N-1])/dL)-1;
		}
		// for D^b theta d/dt u_z
#pragma omp barrier
#pragma omp for nowait
		for (int i = 0;i < N ;i++)
		{
			if (i!=0)
				pr_zd[i]= 0.5*dL*Dg*((pr_w[i]*((Zs->b_U[i+1]-Zs->b_U[i-1])-(Zs->sb_U[i+1]-Zs->sb_U[i-1]))/Zs->tau)+
					((Zs->b_U[i]-Zs->sb_U[i])/Zs->tau)*(pr_w[i+1]-pr_w[i-1])); // d/dx (theta* dz/dt)
			else
				pr_zd[i]= dL*Dg*((pr_w[i]*((Zs->b_U[i+1]-Zs->b_U[i])-(Zs->sb_U[i+1]-Zs->sb_U[i]))/Zs->tau)+
					((Zs->b_U[i]-Zs->sb_U[i])/Zs->tau)*(pr_w[i+1]-pr_w[i])); // d/dx (theta* dz/dt)
		}
#pragma omp barrier
		// precalc linear system values
		__precalc_main(pickard_prevU);
		// D^a H-z
#pragma omp barrier
		if (alpha!=1)
#pragma omp for nowait
		for (int i = 0;i <=N ;i++)
		{
			if (i==0) {pr_D3[i]=0;continue;}
			double v1=pow((i+1)*dL,alpha-1)*SFD(i,pr_D2,alpha,sp_precalc,N);
			double v2=pow(i*dL,alpha-1)*SFD(i-1,pr_D2,alpha,sp_precalc,N);
			if (i==1) v2=0;
			// d/dz (K D^a H)=K_i d/dz(D^aH) + D^a[i]*dK/dz
			// d/dz(D^a H)=(D^a[i]-D^a[i-1])/dL
			// dK/dz=(K[i+1]-K[i-1])/2dL
			double v=Dg*dL*dL*((pr_K[i]*(v1-v2)/dL)+((pr_K[i+1]-pr_K[i-1])/(2*dL))*v1); 
			double vfull=v+Dg*dL*dL*((pr_K[i]*(pow((i+1)*dL,alpha-1)*pr_D2[i]-pow(i*dL,alpha-1)*pr_D2[i-1])/dL)+
					    ((pr_K[i+1]-pr_K[i-1])/(2*dL))*pow((i+1)*dL,alpha-1)*pr_D2[i]); 
			v*=pow(dL,1-alpha)*pow((i+1)*dL,alpha-1);
			pr_D3[i]=vfull;
			if ((i>0)&&(i<N))
				RP[i]+=((pr_R[i]!=0)?v/pr_R[i]:v);
		}
		// D^a D^a H-z
#pragma omp barrier
		if (alpha!=1)
#pragma omp for nowait
		for (int i = 1;i <N ;i++)
		{
			double v=Dg*dL*dL*pow((i+1)*dL,alpha-1)*SFD(i,pr_D3,alpha,sp_precalc,N);
			RP[i]+=((pr_R[i]!=0)?v/pr_R[i]:v);
		}
		// influence of deformations - D^b theta d/dt u_z 
#pragma omp barrier
		if (steady_state==0)
#pragma omp for nowait
		for (int i = 1;i < N ;i++)
		{
			double v=B_*pow((i+1)*dL,alpha-1)*(SFD(i,pr_zd,alpha,sp_precalc,N)+pr_zd[i]);
			RP[i]+=((pr_R[i]!=0)?v/pr_R[i]:v);
		}
#pragma omp barrier
#pragma omp single
		{
			// bottom boundary condition - groundwater level through given h
			RP[N]=L+pr_fl;
			// upper boundary condition - flux
			RP[0]=(Rp(0,b_U,pickard_prevU)+pr_A[0]*2.0*Uc(pickard_prevU))/pr_R[0];  // second order condition through "ghost point"
		}
#pragma omp barrier
	}
	void calculate_velocity()
	{
		for (int i=1;i<N;i++)
		    b_V[i]=pr_K[i]*(((b_U[i+1]-b_U[i-1])/(2.0*dL))-1.0);
		b_V[0]=pr_K[0]*(((b_U[1]-b_U[0])/dL)-1);
		b_V[N]=pr_K[N]*(((b_U[N]-b_U[N-1])/dL)-1);
		for (int i=1;i<N;i++)
		    b_Vd[i]=(b_V[i+1]-b_V[i-1])/(2.0*dL);
		b_Vd[0]=(b_V[1]-b_V[0])/dL;
		b_Vd[N]=(b_V[N]-b_V[N-1])/dL;
	}
	// constructor
	H_solver(int _N,double e,double _L,double _b,double _bb,basic_fr_solver *_Zs) : basic_fr_solver(_L,_N,e),B_(_bb),Zs(_Zs)
	{
		b_V=new double [N+2];
		b_Vd=new double [N+2];
		alpha=_b;
		Dg=pow(dL,-1-alpha);
		precalc_values();
		// initial conditions
		for (int i = 0;i < N + 1;i++)
			b_U[i] = L+pr_fl-(N-i)*dL;
		out=fopen("out_H.txt","wt");
		v_out=fopen("out_V.txt","wt");
		printf("H_solver: N %d L %g alpha %g B %g\n",N,L,alpha,B_);
	}
	// desctructor
	~H_solver()
	{
		delete b_V;
		delete b_Vd;
		if (pr_K) delete [] pr_K;
		if (pr_w) delete [] pr_w;
		if (pr_D2) delete [] pr_D2;
		if (pr_D3) delete [] pr_D3;
		if (pr_zd) delete [] pr_zd;
		if (pr_dwdh) delete [] pr_dwdh;
		if (vgm_ns) delete [] vgm_ns;
		if (vgm_s0s) delete [] vgm_s0s;
		if (vgm_s1s) delete [] vgm_s1s;
		if (vgm_as) delete [] vgm_as;
		if (vgm_h0) delete [] vgm_h0;
		if (vgm_h1) delete [] vgm_h1;
		if (vgm_k) delete [] vgm_k;
		if (vgm_power) delete [] vgm_power;
		if (vgm_cw) delete [] vgm_cw;
		if (pr_vgm_ns) delete [] pr_vgm_ns;
		if (pr_vgm_s0s) delete [] pr_vgm_s0s;
		if (pr_vgm_s1s) delete [] pr_vgm_s1s;
		if (pr_vgm_as) delete [] pr_vgm_as;
		if (pr_vgm_k) delete [] pr_vgm_k;
		if (pr_vgm_power) delete [] pr_vgm_power;
		if (pr_vgm_cw) delete [] pr_vgm_cw;
		if (EV_T) delete [] EV_T;
		if (EV_F) delete [] EV_F;
		if (perc_T) delete [] perc_T;
		if (perc_A) delete [] perc_A;
		if (flT) delete [] flT;
		if (flV) delete [] flV;
		fclose(out);
		fclose(v_out);
	}
	// output
	void output_solution()
	{
		if (pr_K) calculate_velocity();
		fprintf(out,"t(days) %g tau(seconds) %g ET %g Prec %g avgH %g(%d) Fl %g - H: ",(T-tau)/86400.0,tau,pr_ET,pr_perc,pr_ah,irr_on,pr_fl);
		if (T==tau)
		{
		    for (int i=0;i<N+1;i+=zoutstep)
			fprintf(out,"%g ",i*dL);
		    fprintf(out,"\n");
		    fprintf(out,"t(days) %g tau(seconds) %g ET %g Prec %g avgH %g(%d) Fl %g - H: ",(T-tau)/86400.0,tau,pr_ET,pr_perc,pr_ah,irr_on,pr_fl);
		}
		for (int i=0;i<N+1;i+=zoutstep)
			fprintf(out,"%g ",b_U[i]);
		// calculate water table value (saturated zone starting from the bottom)
		double wt;
		for (wt=N;wt>0;wt--)
			if (b_U[(int)wt]<0)
				break;
		wt*=dL;
		fprintf(out,"water_table %g ",wt);
		// velocities
		fprintf(v_out,"t(days) %g V: ",(T-tau)/86400.0);
		if (T==tau)
		{
		    for (int i=0;i<N+1;i+=zoutstep)
			fprintf(v_out,"%g ",i*dL);
		    fprintf(v_out,"\n");
		    fprintf(v_out,"t(days) %g V: ",(T-tau)/86400.0);
		}
		for (int i=0;i<N+1;i+=zoutstep)
			fprintf(v_out,"%g ",b_V[i]);
		fprintf(v_out,"\n");
		// balance of mass
		if ((pr_w)&&(pr_D2))
		{
		    double total_theta=0;
		    double mult=1.0/Gamma(1+alpha);
		    double mult2=1.0/Gamma(2-alpha);
		    for (int i=0;i<N;i++)
		    {
			double b=pow((N-i)*dL,alpha)-pow((N-(i+1))*dL,alpha);
			total_theta+=0.5*(pr_w[i]+pr_w[i+1])*b*mult2*pow((i+1)*dL,1-alpha);
		    }
		    total_theta*=mult;
		    double flux1=-flow0(b_U);
		    double flux2=pr_w[0]*(Zs->b_U[0]-Zs->sb_U[0])/Zs->tau;
		    if (balance==1e30)
		    {
			 balance=total_theta;
			 flux2=0;
		    }
		    for (int i = 0;i <= N ;i++)
		    {
			if ((i!=0)&&(i!=N))
			    pr_D2[i]= ((b_U[i+1] - b_U[i])/dL)-1;
			pr_D2[0]= flow0(b_U)/(0.5*(pr_K[0]+pr_K[1]));
			pr_D2[N]= ((b_U[N] - b_U[N-1])/dL)-1;
		    }
		    double flux3=pr_K[N]*pow((N+1)*dL,alpha-1)*Dg*dL*dL*(SFD(N,pr_D2,alpha,sp_precalc,N)+pr_D2[N]);
		    fprintf(out,"total_theta %g balance %g fluxE %g fluxDef %g fluxBot %g",total_theta,balance,flux1,flux2,flux3);
		    balance+=(flux1+flux2+flux3)*(T-tau-last_o);
		}
		fprintf(out,"\n");
	}
};
///////////////////////////////////////////////////////
// Solver for strain-stress equation ////////
///////////////////////////////////////////////////////
//dz/dt = v
// ((1-n)r_s+n*r_w) dv/dt = (l+2m)*D_z D_z z - B D_z( theta/n p ) + ((1-n)r_s+n*r_w) g
// dn/dt=D_z((1-n)v)
// [0,...,N-1] - z
// [N,...,2*N-1] - v
// steady state - (l+2m)*d2 z /dz2 = rho g
class Z_solver: public basic_fr_solver {
public:
	// inputs
	double lambda,mu,B_;
	H_solver *Hs; // solver for H equation

	FILE *out2;
	FILE *out3;

	int N1,O1; // last index for z, first for v

	double *pr_D2=NULL;
	double *pr_D3=NULL;
	double *pr_D4=NULL;

	// density of soil particles (from https://doi.org/10.1590/0103-8478cr20160762 , Fig 1)
	double rho_s(int i)
	{
	    if (i>N1-2) i=N1-2;
	    double v0=Hs->sand->theta[i]*Hs->sand->b_U[i]+Hs->sand->b_U[O1+i]+Hs->sand->b_U[O1+N1+1+i];
	    double v1=Hs->silt->theta[i]*Hs->silt->b_U[i]+Hs->silt->b_U[O1+i]+Hs->silt->b_U[O1+N1+1+i];
	    double v2=Hs->clay->theta[i]*Hs->clay->b_U[i]+Hs->clay->b_U[O1+i]+Hs->clay->b_U[O1+N1+1+i];
	    double rel_silt=1000*v1/(v0+v1+v2);
	    return 1000*(1.567-0.0008*rel_silt);
	}
	// coefficients of three-diagonal linear equations system
	double A_(int i,double *pickard_prevU)
	{
		return 0;
	}
	double B(int i,double *pickard_prevU)
	{
		return 0;
	}
	double R(int i,double *pickard_prevU)
	{
		if ((steady_state==1)&&(i<N1))
			return 0.0;
		if ((steady_state==1)&&(i>=O1))
			return 1.0;
		if (i<O1) // z
			return 1.0/tau;
		return ((1-Hs->pr_vgm_s1s[i-O1])*rho_s(i-O1)+Hs->pr_vgm_s1s[i-O1]*rho_w)/tau; //v
	}
	// right part 
	double Rp(int ii,double *b_Uold,double *pickard_prevU)
	{
		if (steady_state==1) 
				return ((ii<N1)?((1-sigma)*rho_s(ii)+sigma*rho_w)*g:0.0);
		if (ii<O1) // z
			return b_Uold[ii]/tau;
		// v
		double ret=((1-Hs->pr_vgm_s1s[ii-O1])*rho_s(ii-O1)+Hs->pr_vgm_s1s[ii-O1]*rho_w)*((b_Uold[ii]/tau)+g);
		return ret;
	}
	// linear solver
	// matrix multiplication on vector UU
	void mmult(double *UU)
	{
		if ((B_==0)&&(Hs->B_==0)) return; // don't solve if Biot coef=0
		mmult_main(UU);
		// z: -v
#pragma omp barrier
		if (steady_state==0)
#pragma omp for nowait
		for (int i=1;i<N1;i++)
			MM[i]-=((pr_R[i]!=0.0)?UU[i+O1]/pr_R[i]:UU[i+O1]);
		// v: -d2 z
#pragma omp barrier
#pragma omp for nowait
		for (int i=1;i<N1;i++)
		{
			double a=pow(i*dL,alpha-1);
			double b=pow((i+1)*dL,alpha-1);
			double v=(lambda+2*mu)*pow(dL,1-alpha)*pow((i+1)*dL,alpha-1)*Dg*(a*UU[i-1] - (a+b)*UU[i] + b*UU[i+1]);
			if (steady_state==0)
				MM[i+O1]-=((pr_R[i+O1]!=0.0)?v/pr_R[i+O1]:v);
			else
				MM[i]-=((pr_R[i]!=0.0)?v/pr_R[i]:v);
		}
#pragma omp barrier
#pragma omp single
		{
		    // u_z
			// bottom boundary condition
			MM[N1]=UU[N1];
			// upper boundary condition
			MM[0]=UU[0]-UU[1];
		    // v_z
			// bottom boundary condition
			MM[N]=UU[N];
			// upper boundary condition
			MM[O1]=UU[O1]-UU[O1+1];
			if (steady_state==1)
				MM[O1]=UU[O1];
		}
#pragma omp barrier
#pragma omp for nowait
		for (int i=1;i<N1;i++)
			if (steady_state==0)
				MM[i+O1]/=(lambda+2*mu)*Dg;
			else
				MM[i]/=(lambda+2*mu)*Dg;
#pragma omp barrier
	}
	// precalculations for current time step
	void __precalc(double *pickard_prevU)
	{
		if ((B_==0)&&(Hs->B_==0)) return; // don't solve if Biot coef=0
#pragma omp single
		{ // precalculate non-linear coefficients using the solution of previous iteration
		    memcpy(sb_U,b_U,sizeof(double)*(N+2));
		    memcpy(b_U,pickard_prevU,sizeof(double)*(N+2));
		    if (pr_A==NULL) pr_A=new double[N+2];
		    if (pr_B==NULL) pr_B=new double[N+2];
		    if (pr_R==NULL) pr_R=new double[N+2];
		    if (pr_D2==NULL) pr_D2=new double[N+2];
		    if (pr_D3==NULL) pr_D3=new double[N+2];
		    if (pr_D4==NULL) pr_D4=new double[N+2];
		    memcpy(b_U,sb_U,sizeof(double)*(N+2));
		}
		// for D_z z
#pragma omp barrier
#pragma omp for nowait
		for (int i = 0;i <= N1 ;i++)
		{
			if (i!=N1)
				pr_D2[i]= Dg*(pickard_prevU[i+1] - pickard_prevU[i]);
			else
				pr_D2[i]=Dg*(pickard_prevU[i]-pickard_prevU[i-1]);
		}
		// for D_z D_z z
#pragma omp barrier
#pragma omp for nowait
		for (int ii = 0;ii <= N1 ;ii++)
		{
		    if (ii==0) 
			pr_D3[ii]=0;
		    else
		    {
			if (ii==1)
			    pr_D3[ii]=pow((ii+1)*dL,alpha-1)*(SFD(ii,pr_D2,alpha,sp_precalc,N1)+pr_D2[ii])-
				   pow((ii)*dL,alpha-1)*pr_D2[ii-1];
			else
			    pr_D3[ii]=pow((ii+1)*dL,alpha-1)*(SFD(ii,pr_D2,alpha,sp_precalc,N1)+pr_D2[ii])-
				   pow((ii)*dL,alpha-1)*(SFD(ii-1,pr_D2,alpha,sp_precalc,N1)+pr_D2[ii-1]);
		    }
		}
		// for D_z(theta/n p)
#pragma omp barrier
		if (Hs->pr_K)
#pragma omp for nowait
		for (int i = 0;i < N1 ;i++)
		{
		    if (i!=0)
			pr_D4[i]=((Hs->pr_w[i+1]/Hs->pr_vgm_s1s[i+1])*Hs->b_U[i+1]-(Hs->pr_w[i-1]/Hs->pr_vgm_s1s[i-1])*Hs->b_U[i-1])/(2.0*dL);
		    else
			pr_D4[i]=((Hs->pr_w[i+1]/Hs->pr_vgm_s1s[i+1])*Hs->b_U[i+1]-(Hs->pr_w[i]/Hs->pr_vgm_s1s[i])*Hs->b_U[i])/dL;
		}
#pragma omp barrier
	    // precalc linear system values
		__precalc_main(pickard_prevU);
		// D^a D^a z
#pragma omp barrier
#pragma omp for nowait
		for (int i = 1;i <N1 ;i++)
		{
		    double v=SFD(i,pr_D3,alpha,sp_precalc,N1)+pr_D3[i];
		    v-=(pow((i+1)*dL,alpha-1)*pr_D2[i]-pow(i*dL,alpha-1)*pr_D2[i-1]);
		    v*=(lambda+2*mu)*pow(dL,1-alpha)*pow((i+1)*dL,alpha-1);
		    if (steady_state==0)
			RP[i+O1]+=((pr_R[i+O1]!=0)?v/pr_R[i+O1]:v);
		    else
			RP[i]+=((pr_R[i]!=0)?v/pr_R[i]:v);
		}
		// D_z(theta/n p)
#pragma omp barrier
		if (Hs->pr_K)
#pragma omp for nowait
		for (int i = 1;i < N1 ;i++)
		{
			double v=B_*rho_w*g*Dg*dL*dL*pow((i+1)*dL,alpha-1)*(SFD(i,pr_D4,alpha,sp_precalc,N1)+pr_D4[i]);
			if (steady_state==0)
				RP[i+O1]-=((pr_R[i+O1]!=0)?v/pr_R[i+O1]:v);
			else
				RP[i]-=((pr_R[i]!=0)?v/pr_R[i]:v);
		}
#pragma omp barrier
#pragma omp single
		{
			// z
			// bottom boundary condition
			RP[N1]=0.0;
			// upper boundary condition
			RP[0]=0;
			if (steady_state==1) RP[0]=0;
			// v
			// bottom boundary condition
			RP[N]=0.0;
			// upper boundary condition
			RP[O1]=0.0;
			if (steady_state==1)
				RP[O1]=0;
		}
#pragma omp barrier
#pragma omp for nowait
		for (int i=1;i<N1;i++)
			if (steady_state==0)
				RP[i+O1]/=(lambda+2*mu)*Dg;
			else
				RP[i]/=(lambda+2*mu)*Dg;
	}
	// constructor
	Z_solver(int _N,double e,double _L,H_solver *_Vs,double _l,double _m,double _a,double _b) : basic_fr_solver(_L,_N,e), Hs(_Vs),lambda(_l),mu(_m),B_(_b)
	{
		out=fopen("out_Z.txt","wt");
		out2=fopen("out_Vz.txt","wt");
		out3=fopen("out_n.txt","wt");
		alpha=_a;
		N1=(N-1)/2;
		O1=N1+1;
		dL = L/N1;
		Dg=pow(dL,-1-alpha);
		for (int i = 0;i < N + 1;i++)
			b_U[i] = 0;
		printf("Z_solver: N %d(%d,%d) L %g alpha %g B %g lambda %g mu %g\n",N,N1,O1,L,alpha,B_,lambda,mu);
	}
	// destructor
	~Z_solver()
	{
		if (pr_D2) delete [] pr_D2;
		if (pr_D3) delete [] pr_D3;
		if (pr_D4) delete [] pr_D4;
		fclose(out);
		fclose(out2);
		fclose(out3);
	}
	// output
	void output_solution()
	{
		fprintf(out,"t(days) %g tau(seconds) %g Z: ",(T-tau)/86400.0,tau);
		fprintf(out2,"t(days) %g tau(seconds) %g V: ",(T-tau)/86400.0,tau);
		fprintf(out3,"t(days) %g tau(seconds) %g n: ",(T-tau)/86400.0,tau);
		if (T==tau)
		{
		    for (int i=0;i<=N1;i+=zoutstep)
			{
				fprintf(out,"%g ",i*dL);
				fprintf(out2,"%g ",i*dL);
				fprintf(out3,"%g ",i*dL);
			}
		    fprintf(out,"\n");
		    fprintf(out,"t(days) %g tau(seconds) %g Z: ",(T-tau)/86400.0,tau);
		    fprintf(out2,"\n");
		    fprintf(out2,"t(days) %g tau(seconds) %g V: ",(T-tau)/86400.0,tau);
		    fprintf(out3,"\n");
		    fprintf(out3,"t(days) %g tau(seconds) %g n: ",(T-tau)/86400.0,tau);
		}
		for (int i=0;i<=N1;i+=zoutstep)
		{
		    fprintf(out,"%g ",b_U[i]);
		    fprintf(out2,"%g ",b_U[i+O1]);
		    fprintf(out3,"%g ",Hs->pr_vgm_s1s[i]);
		}
		// balance of momentum
		if ((Hs->pr_w!=NULL)&&(pr_D2!=NULL))
		{
		    double total_momentum=0,total_rhog=0;
		    double mult=1.0/Gamma(1+alpha);
		    double mult2=1.0/Gamma(2-alpha);
		    for (int i=0;i<N1;i++)
		    {
			double b=pow((N1-i)*dL,alpha)-pow((N1-(i+1))*dL,alpha);
			total_momentum+=0.5*(((1-Hs->pr_vgm_s1s[i])*rho_s(i)+Hs->pr_vgm_s1s[i]*rho_w)*b_U[O1+i]+((1-Hs->pr_vgm_s1s[i+1])*rho_s(i+1)+Hs->pr_vgm_s1s[i+1]*rho_w)*b_U[O1+i+1])*
					 b*mult2*pow((i+1)*dL,1-alpha);
			total_rhog+=0.5*(((1-Hs->pr_vgm_s1s[i])*rho_s(i)+Hs->pr_vgm_s1s[i]*rho_w)*g+((1-Hs->pr_vgm_s1s[i+1])*rho_s(i+1)+Hs->pr_vgm_s1s[i+1]*rho_w)*g)*b*mult2*pow((i+1)*dL,1-alpha);
		    }
		    total_momentum*=mult;
		    total_rhog*=mult;
		    for (int i = 0;i <= N1 ;i++)
		    {
			if (i!=N1)
			    pr_D2[i]= (b_U[i+1] - b_U[i])/dL;
			else
			    pr_D2[N1]= (b_U[N1] - b_U[N1-1])/dL;
		    }
		    double flux1=(lambda+2.0*mu)*pr_D2[0];
		    double flux1a=(lambda+2.0*mu)*pow((N1+1)*dL,alpha-1)*Dg*dL*dL*(SFD(N1,pr_D2,alpha,sp_precalc,N1)+pr_D2[N1]);
		    pr_D4[0]=((Hs->pr_w[1]/Hs->pr_vgm_s1s[1])*Hs->b_U[1]-(Hs->pr_w[0]/Hs->pr_vgm_s1s[0])*Hs->b_U[0])/dL;
		    pr_D4[N1]=((Hs->pr_w[N1]/Hs->pr_vgm_s1s[N1])*Hs->b_U[N1]-(Hs->pr_w[N1-1]/Hs->pr_vgm_s1s[N1-1])*Hs->b_U[N1-1])/dL;
		    double flux2=-B_*rho_w*g*pr_D4[0];
		    double flux2a=-B_*rho_w*g*pr_D4[N1];
		    if (balance==1e30)
			 balance=total_momentum;
		    fprintf(out,"total_momentum %g balance %g flux_sigma %g %g fluxP %g %g total_rhog %g ",total_momentum,balance,flux1,flux1a,flux2,flux2a,total_rhog);
		    balance+=(flux1+flux1a+flux2+flux2a+total_rhog)*(T-tau-last_o);
		}
		fprintf(out,"\n");
		fprintf(out2,"\n");
		fprintf(out3,"\n");
	}
};
////////////////////////////////////////
// reads data from <time,value> pair file into rlT,rlV (nrl - their size) from filename
////////////////////////////////////////
void read_pair_file(double *&rlT,double *&rlV,int &nrl,char *filename)
{
    FILE *fi=fopen(filename,"rt");
    if (fi)
    {
	char str[1024];
	nrl=0;
	while (fgets(str,1024,fi)) nrl++;
	rlT=new double[nrl];
	rlV=new double[nrl];
	nrl=0;
	fseek(fi,0,SEEK_SET);			
	while (fscanf(fi,"%lg %lg\n",rlT+nrl,rlV+nrl)==2)
	    if (rlT[nrl]>=0)
		nrl++;
	printf("%d pairs read from %s\n",nrl,filename);
    }
    else
	printf("error reading %s\n",filename);
}
////////////////////////////////////////
////////////////// main/////////////////
////////////////////////////////////////
int main(int argc,char **argv)
{
	double Tm = 365*24*3600.0; // ending time, s
	Om = 30*24*3600.0; // interval for solution output
	double L=3; // domain depth, m
	int N=50;
	// reading basic parameters
	for (int i=1;i<argc;i+=2)
	{
		// basic parameters
		if (strcmp(argv[i],"Tm")==0) // in days
			Tm=atof(argv[i+1]);
		if (strcmp(argv[i],"Om")==0) // in days
			Om=atof(argv[i+1]);
		if (strcmp(argv[i],"NB")==0)
			N=atoi(argv[i+1]);
		if (strcmp(argv[i],"LB")==0)
			L=atof(argv[i+1]);
		if (strcmp(argv[i],"debug_level")==0)
			debug_level=atoi(argv[i+1]);
		if (strcmp(argv[i],"Zoutstep")==0)
			zoutstep=atoi(argv[i+1]);
		// (v,z) equation parameters
		if (strcmp(argv[i],"lambda")==0)
			lambda=atof(argv[i+1]);
		if (strcmp(argv[i],"mu")==0)
			mu=atof(argv[i+1]);
		if (strcmp(argv[i],"porosity")==0)
			sigma=atof(argv[i+1]);
		if (strcmp(argv[i],"biot")==0)
			B=atof(argv[i+1]);
		if (strcmp(argv[i],"b2")==0)
			B2=atof(argv[i+1]);
		// particle transport parameters
		if (strcmp(argv[i],"pt_l_a0")==0)
			pt_l_a[0]=atof(argv[i+1]);
		if (strcmp(argv[i],"pt_l_a1")==0)
			pt_l_a[1]=atof(argv[i+1]);
		if (strcmp(argv[i],"pt_l_a2")==0)
			pt_l_a[2]=atof(argv[i+1]);
		if (strcmp(argv[i],"pt_l_s0")==0)
			pt_l_s[0]=atof(argv[i+1]);
		if (strcmp(argv[i],"pt_l_s1")==0)
			pt_l_s[1]=atof(argv[i+1]);
		if (strcmp(argv[i],"pt_l_s2")==0)
			pt_l_s[2]=atof(argv[i+1]);
		if (strcmp(argv[i],"pt_k_det0")==0)
			pt_k_det[0]=atof(argv[i+1]);
		if (strcmp(argv[i],"pt_k_det1")==0)
			pt_k_det[1]=atof(argv[i+1]);
		if (strcmp(argv[i],"pt_k_det2")==0)
			pt_k_det[2]=atof(argv[i+1]);
		if (strcmp(argv[i],"pt_c00")==0)
			pt_c0[0]=atof(argv[i+1]);
		if (strcmp(argv[i],"pt_c01")==0)
			pt_c0[1]=atof(argv[i+1]);
		if (strcmp(argv[i],"pt_c02")==0)
			pt_c0[2]=atof(argv[i+1]);
		// fractional derivative parameters
		if (strcmp(argv[i],"alpha")==0)
			alpha=atof(argv[i+1]);
		if (strcmp(argv[i],"beta")==0)
			beta=atof(argv[i+1]);
		// linear solver parameters
		if (strstr(argv[i],"ls_eps")!=NULL)
		    ls_eps=atof(argv[i+1]);
		if (strstr(argv[i],"ls_e2")!=NULL)
		    ls_eps2=atof(argv[i+1]);
		if (strstr(argv[i],"ls_e3")!=NULL)
		    ls_eps3=atof(argv[i+1]);
		if (strstr(argv[i],"pick_eps")!=NULL)
		    pick_eps=atof(argv[i+1]);
		if (strstr(argv[i],"pick_move")!=NULL)
		    pick_move=atof(argv[i+1]);
		if (strstr(argv[i],"fr_eps")!=NULL)
		    fr_eps=atof(argv[i+1]);
		if (strstr(argv[i],"ls_max_iter")!=NULL)
		    init_ls_max_iter=ls_max_iter=atoi(argv[i+1]);
		if (strstr(argv[i],"ls_percent")!=NULL)
		    ls_percent=atof(argv[i+1]);
		if (strstr(argv[i],"ls_min_tau")!=NULL) //  in seconds
		    ls_min_tau=atof(argv[i+1]);
		if (strstr(argv[i],"ls_mult")!=NULL)
		    ls_mult=atof(argv[i+1]);
		if (strstr(argv[i],"max_tau")!=NULL) // in seconds
		    max_tau=atof(argv[i+1]);
		if (strstr(argv[i],"st_tau")!=NULL) // in seconds
		    st_tau=atof(argv[i+1]);
		// irrigation modelling parameters
		if (strstr(argv[i],"irr_depth")!=NULL) // in seconds
		    irr_depth=atof(argv[i+1]);
		if (strstr(argv[i],"irr_p0")!=NULL) // in seconds
		    irr_p0=atof(argv[i+1]);
		if (strstr(argv[i],"irr_p1")!=NULL) // in seconds
		    irr_p1=atof(argv[i+1]);
		if (strstr(argv[i],"irr_flow")!=NULL) // in seconds
		    irr_flow=atof(argv[i+1]);
		if (strstr(argv[i],"irr_time0")!=NULL) // in seconds
		    irr_time0=atof(argv[i+1]);
		if (strstr(argv[i],"irr_time1")!=NULL) // in seconds
		    irr_time1=atof(argv[i+1]);
		if (strstr(argv[i],"irr_tau")!=NULL) // in seconds
		    irr_tau=atof(argv[i+1]);
		if (strstr(argv[i],"is_irr")!=NULL) // in seconds
		    is_irr=atoi(argv[i+1]);
	}
	printf("Grid size - %d, Tend %g Tsave %g \n",N,Tm,Om);
	// creating solvers
	H_solver *solver=new H_solver(N,ls_eps,L,beta,B2,NULL);
	Z_solver *Zsolver=new Z_solver(2*N+1,ls_eps2,L,solver,lambda,mu,alpha,B);
	solver->Zs=Zsolver;
	// reading inputs and passing them to the solver
	int check=0;
	for (int i=1;i<argc;i+=2)
	{
		// reading the file with van Genuchtem-Mualem model's parameters
		if (strstr(argv[i],"vgm")!=NULL)
		{
		    FILE *fi=fopen(argv[i+1],"rt");
		    if (fi)
		    {
			char str[1024];
			vgm_nlayers=-1;
			while (fgets(str,1024,fi)) vgm_nlayers++;
			vgm_ns=new double[vgm_nlayers];
			vgm_as=new double[vgm_nlayers];
			vgm_s0s=new double[vgm_nlayers];
			vgm_s1s=new double[vgm_nlayers];
			vgm_h0=new double[vgm_nlayers];
			vgm_h1=new double[vgm_nlayers];
			vgm_k=new double[vgm_nlayers];
			vgm_power=new double[vgm_nlayers];
			vgm_cw=new double[vgm_nlayers];
			soil_lambda=new double[vgm_nlayers];
			soil_Ct=new double[vgm_nlayers];
			for (int i=0;i<4;i++) soil_D[i]=new double[vgm_nlayers];
			vgm_nlayers=0;
			fseek(fi,0,SEEK_SET);
			fgets(str,1024,fi);
			while(fgets(str,1024,fi))
			{
				sscanf(str,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",vgm_s0s+vgm_nlayers,vgm_s1s+vgm_nlayers,
				vgm_ns+vgm_nlayers,vgm_as+vgm_nlayers,vgm_h0+vgm_nlayers,vgm_h1+vgm_nlayers,vgm_k+vgm_nlayers,vgm_power+vgm_nlayers,
				vgm_cw+vgm_nlayers,soil_lambda+vgm_nlayers,soil_Ct+vgm_nlayers,
				soil_D[0]+vgm_nlayers,soil_D[1]+vgm_nlayers,soil_D[2]+vgm_nlayers,soil_D[3]+vgm_nlayers);
				vgm_nlayers++;
			}
			printf("VGM data read (%d layers)\n",vgm_nlayers);
			check++;
		    }
		}
		// ET
		if (strstr(argv[i],"et_file")!=NULL)
		{
			read_pair_file(solver->EV_T,solver->EV_F,solver->nEV_T,argv[i+1]);
			check++;
		}
		// precipitation
		if (strstr(argv[i],"prec_file")!=NULL)
		{
			read_pair_file(solver->perc_T,solver->perc_A,solver->nperc,argv[i+1]);
			check++;
		}
		// groundwater level
		if (strstr(argv[i],"gw_file")!=NULL)
		{
			read_pair_file(solver->flT,solver->flV,solver->nfl,argv[i+1]);
			check++;
		}
	}
	if (check!=4)
	{
		printf("not all input file were given\n");
		return 1;
	}
	// initialize H solver and create particle transport solvers
	solver->soil_calc_coefs();
	solver->pr_wetness();
	python_initialize(N);
	particle_transport *pt1=new particle_transport(3*N,L/N,solver->tau,pt_l_a[0],pt_l_s[0],pt_k_det[0],solver->b_V,solver->pr_w,pt_c0[0],0,ls_eps3);
	particle_transport *pt2=new particle_transport(3*N,L/N,solver->tau,pt_l_a[1],pt_l_s[1],pt_k_det[1],solver->b_V,solver->pr_w,pt_c0[1],1,ls_eps3);
	particle_transport *pt3=new particle_transport(3*N,L/N,solver->tau,pt_l_a[2],pt_l_s[2],pt_k_det[2],solver->b_V,solver->pr_w,pt_c0[2],2,ls_eps3);
	solver->sand=pt1;
	solver->silt=pt2;
	solver->clay=pt3;
	basic_solver *others[4]={Zsolver,pt1,pt2,pt3};
	// solve steady state problem for H
	    double smi=ls_max_iter;
	    solver->steady_state=1;
	    Zsolver->steady_state=1;
	    ls_max_iter*=20;
	    solver->solve_up_to_T(solver->tau,solver->tau,others,1,1);
	    solver->steady_state=0;
	    Zsolver->steady_state=0;
	    solver->T=solver->tau;
	    Zsolver->T=Zsolver->tau;
	    ls_max_iter=smi;
	// run simulation
	solver->output_solution();
	Zsolver->output_solution();
	pt1->output_solution();
	pt2->output_solution();
	pt3->output_solution();

	solver->solve_up_to_T(Tm,solver->tau,others,4,1);

	solver->output_solution();
	Zsolver->output_solution();
	pt1->output_solution();
	pt2->output_solution();
	pt3->output_solution();

	python_finalize();
	return 0;
}
