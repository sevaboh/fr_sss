#include "gamma.h"

//////////////////////////////////////
///// global variables ///////////////
//////////////////////////////////////
const double HtoPmult=9.8;

int debug_level=0;
int zoutstep=1;
int testing_mode=0; // testing with constructed solution (1)
// linear solver (TFQMR)
int ls_max_iter=30; // maximal number of iterations
int init_ls_max_iter=30;
double ls_eps=1e-12; // accuracy threshold (H)
double ls_eps2=1e-14; // accuracy threshold (Z)
double pick_eps=1e-15; // accuracy threshold for pickard iteration
double pick_move=1;
double fr_eps=1e-20; // accuracy for numerical integration
double ls_min_tau=0.001; // minimal time step length, s
double max_tau=1*24*3600.0; // maximal time step length, s
double ls_mult=1.25; // time step multiplier/divisor
double st_tau=0.01;
double ls_percent=0.66; // threshold in per cent of max.iterations used to make a decision to change time step length
// soil parameters
// VGM parameters
double *vgm_ns,*vgm_s0s,*vgm_s1s,*vgm_as,*vgm_h0,*vgm_h1,*vgm_k,*vgm_power,*vgm_cw; // VGM coefs per layer of soil
// strain-stress parameters
double rho_s=1940; //density of soil particles kg/m^3
double rho_w=1000; // density of water
double lambda=19788461.54;
double mu=13192307.69;
double B=0.4,B2=0.4;
double sigma=0.6; // starting porosity
const double g=9.81;
// fractional orders
double alpha=1;
double beta=1;

int vgm_nlayers; // number of layers

double Om; // output interval
double last_o=0; // last output time

double *soil_lambda,*soil_Ct,*soil_D[4]; // to read from input file, not used later
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

///////////////////////////////////////
////////// basic solver class /////////
///////////////////////////////////////

class basic_solver;
class basic_solver {
public:
	int N;
	// storage
	double *b_U; // solution (pressure head/concentration/temperature)
	double *sb_U; // to save solution
	double *MM,*RP; // multiplication results and right part for TFQMR
	// matrix parts
	double *pr_A=NULL,*pr_B=NULL,*pr_R=NULL;
	// soil parameters
	// VGM parameters
	double *pr_vgm_ns=NULL,*pr_vgm_s0s,*pr_vgm_s1s,*pr_vgm_as,*pr_vgm_k,*pr_vgm_power,*pr_vgm_cw; // VGM coefs per cell
	double *pr_soil_lambda;
	double *pr_soil_Ct;
	double *pr_soil_D[4];

	// steps and auxiliary
	double tau; // current time step length 
	double T; // current time
	double L,dL; // domain depth and space variable step length
	int steady_state=0;
	double ls_eps;
	int iter_counter;
	// output
	FILE *out;
	// auxiliary
	// precalculate soil parameters for each cell in pr_*
	void soil_calc_coefs()
	{
	    if (pr_vgm_ns==NULL)
	    {
		pr_vgm_ns=new double[N+2];
		pr_vgm_s0s=new double[N+2];
		pr_vgm_s1s=new double[N+2];
		pr_vgm_as=new double[N+2];
		pr_vgm_k=new double[N+2];
		pr_vgm_power=new double[N+2];
		pr_vgm_cw=new double[N+2];
		pr_soil_lambda=new double[N+2];
		pr_soil_Ct=new double[N+2];
		for (int i=0;i<4;i++)
			pr_soil_D[i]=new double[N+2];
		memset(pr_vgm_s0s,0,(N+2)*sizeof(double));
		// calculate
		for (int i=0;i<=N;i++)
		{
		double vgm_s0,vgm_s1,vgm_a,vgm_n,k,avr_power,vgm_ss,soil_l,soil_ct,sD[4];
	        avr_power=3.5;
		if (vgm_nlayers!=0)
		{
		    // i inside a layer
		    for (int j=0;j<vgm_nlayers;j++)
		    if (((i*dL)>=vgm_h0[j])&&((i*dL)<vgm_h1[j]))
		    {
			vgm_s0=vgm_s0s[j];
			vgm_s1=vgm_s1s[j];
			vgm_a=vgm_as[j];
			vgm_n=vgm_ns[j];
			vgm_ss=vgm_cw[j];
			soil_l=soil_lambda[j];
			soil_ct=soil_Ct[j];
			for (int kk=0;kk<4;kk++) sD[kk]=soil_D[kk][j];
			k=vgm_k[j];
			avr_power=vgm_power[j];
			goto save;
		    }
		    double h1,h0,minh0=1e300,maxh1=0;
		    int i1,i0;
		    int mi1,mi0;
		    int first=1;
		    for (int j=0;j<vgm_nlayers;j++)
		    {
			if (vgm_h0[j]<minh0)
			{
			    minh0=vgm_h0[j];
			    mi0=j;
			}
			if (vgm_h1[j]>maxh1)
			{
			    maxh1=vgm_h1[j];
			    mi1=j;
			}
			for (int k=0;k<vgm_nlayers;k++)
			if (j!=k)
			if (((i*dL)>=vgm_h1[j])&&((i*dL)<vgm_h0[k]))
			{
			    if (first)
			    {
				h1=vgm_h1[j];
				i1=j;
				h0=vgm_h0[k];
				i0=k;
				first=0;
			    }
			    if (vgm_h1[j]>h1)
			    {
				h1=vgm_h1[j];
				i1=j;
			    }
			    if (vgm_h0[k]<h0)
			    {
				h0=vgm_h0[k];
				i0=k;
			    }
			}
		    }
		    // i between two layers
		    if (first==0)
		    {
			double _k=((i*dL)-h1)/(h0-h1);
			vgm_s0=_k*vgm_s0s[i0]+(1.0-_k)*vgm_s0s[i1];
			vgm_s1=_k*vgm_s1s[i0]+(1.0-_k)*vgm_s1s[i1];
			vgm_a=_k*vgm_as[i0]+(1.0-_k)*vgm_as[i1];
			vgm_n=_k*vgm_ns[i0]+(1.0-_k)*vgm_ns[i1];
			vgm_ss=_k*vgm_cw[i0]+(1.0-_k)*vgm_cw[i1];
			soil_l=_k*soil_lambda[i0]+(1.0-_k)*soil_lambda[i1];
			soil_ct=_k*soil_Ct[i0]+(1.0-_k)*soil_Ct[i1];
			for (int kk=0;kk<4;kk++) sD[kk]=_k*soil_D[kk][i0]+(1.0-_k)*soil_D[kk][i1];
			k=_k*vgm_k[i0]+(1.0-_k)*vgm_k[i1];
		    }
		    else
		    {
			// i below minimal
			if ((i*dL)<=minh0)
			{
			    vgm_s0=vgm_s0s[mi0];
			    vgm_s1=vgm_s1s[mi0];
			    vgm_a=vgm_as[mi0];
			    vgm_n=vgm_ns[mi0];
			    vgm_ss=vgm_cw[mi0];
			    soil_l=soil_lambda[mi0];
			    soil_ct=soil_Ct[mi0];
			    for (int kk=0;kk<4;kk++) sD[kk]=soil_D[kk][mi0];
			    k=vgm_k[mi0];
			    avr_power=vgm_power[mi0];
			}
			else
			{
			    // i above maximal
			    if ((i*dL)>=maxh1)
			    {
				vgm_s0=vgm_s0s[mi1];
				vgm_s1=vgm_s1s[mi1];
				vgm_a=vgm_as[mi1];
				vgm_n=vgm_ns[mi1];
				vgm_ss=vgm_cw[mi1];
				soil_l=soil_lambda[mi1];
				soil_ct=soil_Ct[mi1];
				for (int kk=0;kk<4;kk++) sD[kk]=soil_D[kk][mi1];
				k=vgm_k[mi1];
				avr_power=vgm_power[mi1];
			    }
			}
		    }
	    	}
		// save
save:
		pr_vgm_s0s[i]=vgm_s0;
		pr_vgm_s1s[i]=vgm_s1;
		pr_vgm_as[i]=vgm_a;
		pr_vgm_ns[i]=vgm_n;
		pr_vgm_cw[i]=vgm_ss;
		pr_vgm_k[i]=k;
		pr_vgm_power[i]=avr_power;
		pr_soil_lambda[i]=soil_l;
		pr_soil_Ct[i]=soil_ct;
		for (int kk=0;kk<4;kk++) pr_soil_D[kk][i]=sD[kk];
		}
	    }
	}
	// linearly interpolate a value from a list of [T,V] pairs for currebt T 
	double linear_from_pair(double *TaT,double *TaF,int nTa,double __T=-1)
	{
		double ret=0.0;
		if (steady_state==1) return TaF[0];
		if (nTa != 0)
		{
			// calc linear conbination
			int i;
			double p = 0.0;
			double _T=T;
			if (__T!=-1) _T=__T;
			while (_T>TaT[nTa-1]) _T-=TaT[nTa-1];		
			for (i = 0;i<nTa;i++)
				if (TaT[i]>_T)
					break;
			if ((i == 0) || (i == nTa))
			{
				if (i == nTa) i--;
				p= TaF[i];
			}
			else
			{
				double k = (_T - TaT[i - 1]) / (TaT[i] - TaT[i - 1]);
				p= k*TaF[i] + (1 - k)*TaF[i - 1];
			}
			ret=p;
		}
		return ret;
	}
	// linear solver
	// linear system coefficients - A*U(i-1)+R*U(i)+B*U(i+1)=Rp
	virtual double A_(int i,double *pickard_prevU)=0;
	virtual double B(int i,double *pickard_prevU)=0;
	virtual double R(int i,double *pickard_prevU)=0;
	virtual double Rp(int i,double *U,double *pickard_prevU)=0;
	// matrix multiplication of vector UU
	virtual void mmult(double *UU)=0;
	void mmult_main(double *UU)
	{
#pragma omp for nowait
		for (int i = 1;i < N ;i++)
		{
			double rr=pr_R[i];
			MM[i]=UU[i-1]*pr_A[i]+UU[i+1]*pr_B[i]+UU[i]*rr;
			// normalization
			if (rr!=0.0)
				MM[i]/=rr;
		}
	}
	// precalculations for current time step
	virtual void __precalc(double *pickard_prevU)=0;
	void __precalc_main(double *pickard_prevU)
	{
#pragma omp for nowait
		for (int i = 0;i <= N ;i++)
		{
			pr_A[i]=A_(i,pickard_prevU);
			pr_B[i]=B(i,pickard_prevU);
			pr_R[i]=R(i,pickard_prevU);
			// right part
			RP[i]=Rp(i,b_U,pickard_prevU);
			// normalization
			if (pr_R[i]!=0.0)
				RP[i]/=pr_R[i];
		}
	}
	// linear solver (TFQMR) (returns 1 if tau was increased at the end)
	int calc_step(double *pickard_prevU=NULL,int fix_tau=0)
	{
		if (pickard_prevU==NULL) pickard_prevU=b_U;
		// first step
		if (T==tau)
		{
		    for (int i =0;i <= N+1 ;i++)
		    {
			MM[i]=0.0;
			RP[i]=0.0;
		    }
		}
		//////////
		double *w=new double[(N+2)];
		double *y[2];
		y[0]=new double[(N+2)];
		y[1]=new double[(N+2)];
		double *rr=new double[(N+2)];
		double *v=new double[(N+2)];
		double *d=new double[(N+2)];
		double *x=new double[(N+2)];
		double theta=0,nu=0,tt=0,ro=0,ro1=0,c=0,b=0;
		double sgm,aa,rv;
		int n,m;
start:
		if (debug_level>=2) printf("--------- %g ----------\n",tau);
		theta=nu=tt=ro=ro1=c=b=0.0;
		memcpy(x,b_U,(N+2)*sizeof(double));
#pragma omp parallel
	    {
		////////// precalculations
		__precalc(pickard_prevU);
		mmult(x);
#pragma omp for reduction(+:tt)
		for (int i=0;i<(N+2);i++)
		{
			// w_1=y_1=r_0=b-Ax_0
			// rr0: ro=rr_0*r_0!=0 - rr0=r0
			w[i]=y[0][i]=RP[i]-MM[i];
			// d=0
			d[i]=0;
			// tau=||r0||
			tt+=w[i]*w[i];
		}
#pragma omp single
		tt=sqrt(tt);
		if (tt==0.0)
			goto eloop1;
#pragma omp barrier
		// random rr0, ro0=rr_0*r_0
#pragma omp for reduction(+:ro)
		for (int i=0;i<(N+2);i++)
		{
			rr[i]=tt*((rand() % 10000) / (10000.0 - 1.0));
			ro+=rr[i]*w[i];
		}
		// v=Ay_1
		mmult(y[0]);
#pragma omp for
		for (int i=0;i<(N+2);i++)
			v[i]=MM[i];
#pragma omp single
		n=1;
#pragma omp barrier
loop1:
		{
			int br=0;
			// sigma_n-1 - rr_0*v_n-1
#pragma omp single
			sgm=0;
#pragma omp barrier
#pragma omp for reduction(+:sgm)
			for (int i=0;i<(N+2);i++)
				sgm+=rr[i]*v[i];
			// a_n-1=po_n-1/sigma_n-1
#pragma omp single
			aa=ro/sgm;
#pragma omp barrier
			// y_2n=y_2n-1 - a_n-1 * v_n-1
#pragma omp for
			for (int i=0;i<(N+2);i++)
				y[1][i]=y[0][i]-aa*v[i];
#pragma omp single
			m=2*n-1;
#pragma omp barrier
loop2:
			{
				double ot=theta,onu=nu;
				// w_m+1=w_m-a_n-1Ay_m
				mmult(y[m-(2*n-1)]);
#pragma omp for
				for (int i=0;i<(N+2);i++)
					w[i]=w[i]-aa*MM[i];
				// theta_m=||w_m+1||/tau_m-1; c_m=1/sqrt(1+theta_m^2)
#pragma omp single
				theta=0;
#pragma omp barrier
#pragma omp for reduction(+:theta)
				for (int i=0;i<(N+2);i++)
					theta+=w[i]*w[i];
#pragma omp single
			    {
				theta=sqrt(theta)/tt;
				c=1.0/sqrt(1.0+theta*theta);
				// tau_m=tau_m-1 * theta_m *c_m; nu_m=c_m^2 *a_n-1
				tt=tt*theta*c;
				nu=c*c*aa;
				rv=0.0;
			    }
#pragma omp barrier
				// d_m = y_m+(theta_m-1^2 nu_m-1 / a_n-1)*d_m-1
				// x_m=x_m-1 + nu_m *d_m
#pragma omp for
				for (int i=0;i<(N+2);i++)
				{
					d[i]=y[m-(2*n-1)][i]+d[i]*(ot*ot*onu/aa);
					x[i]=x[i]+nu*d[i];
				}
				mmult(x);
#pragma omp for reduction(+:rv)
				for (int i=0;i<(N+2);i++)
					rv+=(RP[i]-MM[i])*(RP[i]-MM[i]);
#pragma omp single
				rv=sqrt(rv)/((N+2));
#pragma omp barrier
				if (rv<ls_eps)
				{
				    br=1;
				    goto eloop2;
				}
				if (!isfinite(rv))
					goto eloop1;
			}
#pragma omp single
			m++;
#pragma omp barrier
			if (m<=2*n)
			    goto loop2;
eloop2:
			if (br==1)
				goto eloop1;
			// ro_n=rr0*w_2n+1, beta_n=ro_n/ro_n-1
#pragma omp single
			ro1=0;
#pragma omp barrier
#pragma omp for reduction(+:ro1)
			for (int i=0;i<(N+2);i++)
				ro1+=rr[i]*w[i];
#pragma omp single
		    {
			b=ro1/ro;
			ro=ro1;
		    }
#pragma omp barrier
			// y_2n+1 = w_2n+1+beta_n*y_2n
#pragma omp for
			for (int i=0;i<(N+2);i++)
				y[0][i]=w[i]+b*y[1][i];
			// v_n=Ay_2n+1+b*(Ay_2n+b*v_n-1)
			mmult(y[1]);
#pragma omp for
			for (int i=0;i<(N+2);i++)
				v[i]=b*(MM[i]+b*v[i]);
			mmult(y[0]);
#pragma omp for
			for (int i=0;i<(N+2);i++)
				v[i]=MM[i]+v[i];
		}
#pragma omp single
		n++;
#pragma omp barrier
		if (n<ls_max_iter)
		    goto loop1;
eloop1:;
	    }
		iter_counter+=n;
		// change tau and recalc if needed
		if (fix_tau==0)
		if (n==ls_max_iter)
		{
			if (tau<ls_min_tau)
		     {
			 x[0]=NAN;
			 printf("minimal tau value reached %g %g\n",tau, ls_min_tau);
			 exit(2);
		     }
		     else
		     {
		    	T-=tau;
		    	tau/=ls_mult;
		    	T+=tau;
			// debug output
			if (debug_level>=3) printf("%p r n %d r %g T %g tau %g sol %g %g %g %g %g %g\n",this,n,rv,T,tau,x[0],x[1],x[2],b_U[0],b_U[1],b_U[2]);
		        goto start;
		     }
		}
		// save solution and free memory
		memcpy(b_U,x,(N+2)*sizeof(double));
		delete [] w;
		delete [] y[0];
		delete [] y[1];
		delete [] rr;
		delete [] v;
		delete [] d;
		delete [] x;
		// increase tau if needed and increase T
		int ret=0;
		if (fix_tau==0)
		if (n<ls_max_iter*ls_percent)
			if (tau<max_tau)
			{
				 tau*=ls_mult;
				 ret=1;
			}
		T+=tau;
		// debug output
		if (debug_level>=3) printf("%p e n %d rv %g T %g tau %g - ret %d fix %d mt %g nthr %g\n",this,n,rv,T,tau,ret,fix_tau,max_tau,ls_max_iter*ls_percent);
		return ret;
	}
	// perform Pickard iteration 
	void pickard_calc_step(basic_solver **other=NULL,int nother=0)
	{
	    double *prev_U=new double[N+2];
	    double *init_U=new double[N+2];
	    double **prev_other_U=NULL,**init_other_U=NULL;
	    double diff=0,prev_diff=1e30;
	    int iter=0;
	    int mul=0;
	    int *other_muls=NULL;
	    iter_counter=0;
	    memcpy(init_U,b_U,sizeof(double)*(N+2));
restart:
	    memcpy(b_U,init_U,sizeof(double)*(N+2));
	    prev_diff=1e30;
	    if (init_other_U)
	    for (int i=0;i<nother;i++)
		memcpy(other[i]->b_U,init_other_U[i],sizeof(double)*(other[i]->N+2));
	    iter=0;
		if (debug_level>=3) printf("%p pickard init 0 tau %g T %g other %d \n",this,tau,T-tau,nother);
	    mul=calc_step();
	    if (mul)
	    {	
		T-=tau;
		tau/=ls_mult;
		T+=tau;
	    }
		if (debug_level>=3) printf("%p pickard init 1 tau %g T %g other %d %d\n",this,tau,T-tau,nother,mul);
	    if (nother) // run other solvers
	    {
		if (prev_other_U==NULL) { prev_other_U=new double *[nother]; memset(prev_other_U,0,nother*sizeof(void *));}
		if (init_other_U==NULL) { init_other_U=new double *[nother]; memset(init_other_U,0,nother*sizeof(void *));}
		if (other_muls==NULL) other_muls=new int[nother];
		for (int i=0;i<nother;i++)
		{
			if (prev_other_U[i]==NULL) prev_other_U[i]=new double [other[i]->N+2];
			if (init_other_U[i]==NULL) init_other_U[i]=new double [other[i]->N+2];
			memcpy(init_other_U[i],other[i]->b_U,sizeof(double)*(other[i]->N+2));
			other_muls[i]=other[i]->calc_step();
			if (other_muls[i])
			{	
				other[i]->T-=other[i]->tau;
				other[i]->tau/=ls_mult;
				other[i]->T+=other[i]->tau;
		        }
		}
		// find minimal tau
		double min_tau=tau;
		for (int i=0;i<nother;i++)
			if (other[i]->tau<min_tau)
				min_tau=other[i]->tau;
		// solve all once more with fixed minimal tau
		memcpy(b_U,init_U,sizeof(double)*(N+2));
		for (int i=0;i<nother;i++)
			memcpy(other[i]->b_U,init_other_U[i],sizeof(double)*(other[i]->N+2));
		T-=2.0*tau;
		tau=min_tau;
		T+=tau;
		calc_step(NULL,1);
		for (int i=0;i<nother;i++)
		{
			other[i]->T-=2.0*other[i]->tau;
			other[i]->tau=min_tau;
			other[i]->T+=other[i]->tau;
			other[i]->calc_step(NULL,1);
		}
	    }
	    // check for not finite
	    int is_nf=0;
	    for (int i=0;i<N;i++)
		if (!isfinite(b_U[i]))
		{
		    is_nf=1;
		    break;
		}
	    if (is_nf==0)
	    for (int j=0;j<nother;j++)
	    {
		for (int i=0;i<other[j]->N;i++)
		if (!isfinite(other[j]->b_U[i]))
		{
		    is_nf=1;
		    break;
		} 
		if (is_nf)
		    break;
	    }
	    if (is_nf) // decrease tau and restart
		goto decrease;
	    // debug
	    if (debug_level>=2)
	    {
		 printf("%p pickard init tau %g T %g other %d ",this,tau,T-tau,nother);
		for (int i=0;i<nother;i++) printf(" %d tau %g T %g ",i,other[i]->tau,other[i]->T-other[i]->tau);
		printf("\n");
	    }
	    do
	    {
		// save solution on previous iteration and restore initial values
		memcpy(prev_U,b_U,sizeof(double)*(N+2));
		memcpy(b_U,init_U,sizeof(double)*(N+2));
		// solve on next iteration
		T-=tau;
		calc_step(prev_U,1);
		for (int i=0;i<N+1;i++) b_U[i]=prev_U[i]+pick_move*(b_U[i]-prev_U[i]); // move to new b_U
		for (int i=0;i<nother;i++)
		{
			memcpy(prev_other_U[i],other[i]->b_U,sizeof(double)*(other[i]->N+2));
			memcpy(other[i]->b_U,init_other_U[i],sizeof(double)*(other[i]->N+2));
		    other[i]->T-=other[i]->tau;
		    other[i]->calc_step(prev_other_U[i],1);
		    for (int j=0;j<other[i]->N+1;j++) 
			other[i]->b_U[j]=prev_other_U[i][j]+pick_move*(other[i]->b_U[j]-prev_other_U[i][j]); // move to new b_U
		}
		// calculate difference
		diff=0;
		int sn=N;
		for (int i=0;i<N+1;i++)
		    diff+=(b_U[i]-prev_U[i])*(b_U[i]-prev_U[i]);
		for (int j=0;j<nother;j++)
		{
		    for (int i=0;i<other[j]->N+1;i++)
			diff+=(other[j]->b_U[i]-prev_other_U[j][i])*(other[j]->b_U[i]-prev_other_U[j][i]);
		    sn+=other[j]->N;
		}
		diff=diff/sn;
		// debug
		if (debug_level>=2) printf("%p pickard iter %d diff %g T %g S %g %g %g - %g %g %g\n",this,iter,diff,T-tau,b_U[0],b_U[1],b_U[2],b_U[N-2],b_U[N-1],b_U[N]);
		if (diff>prev_diff) goto decrease;
		prev_diff=diff;
	    if (!isfinite(diff))
			break;
		if ((iter++)> ls_max_iter) 
			break;
	    }
	    while (diff>pick_eps);
	    if (!isfinite(diff))
	    {
			ls_max_iter/=2.0;
			goto restart;
	    }
	    if (iter>=ls_max_iter) // decrease step
	    {
decrease:
		 if (tau<ls_min_tau)
		 {
			 b_U[0]=NAN;
			 printf("minimal tau value reached\n");
			 exit(2);
		 }
		 else
		 {
		    T-=tau;
		    tau/=ls_mult;
		    T+=tau;
		    if (debug_level>=2) printf("%p pickard restart niter %d diff %g tau %g T %g\n",this,iter,diff,tau,T-tau);
		    goto restart;
		}
	    }
	    // increase tau if needed
	    if (mul)
	    {
		T-=tau;
		tau*=ls_mult;
		T+=tau;
	    }
	    if (nother)
		for (int i=0;i<nother;i++)
		    if (other_muls[i])
		    {
			other[i]->T-=other[i]->tau;
			other[i]->tau*=ls_mult;
			other[i]->T+=other[i]->tau;
		    }
	    // debug
	    if (debug_level>=1) printf("%p pickard niter %d diff %g tau %g T %g iters %d\n",this,iter,diff,tau,T-tau,iter_counter);
	    // clean up
	    delete [] prev_U;
	    delete [] init_U;
	    if (nother)
	    {
		for (int i=0;i<nother;i++)
		{
		    delete [] prev_other_U[i];
		    delete [] init_other_U[i];
		}
		delete [] prev_other_U;
		delete [] init_other_U;
		delete [] other_muls;
	    }
		// reset max iter
		ls_max_iter=init_ls_max_iter;
	}
	// solve up to the time point eT
	void solve_up_to_T(double eT,double _tau,basic_solver **other=NULL,int nother=0,int output=0)
	{
		double prev_tau=-1;
		while ((T-tau)<(eT-ls_min_tau))
		{
		    if (T+tau>eT)
		    {
			prev_tau=tau;
			T-=tau;
			tau=eT-T;
			T+=tau;
		    }
		    pickard_calc_step(other,nother);
		    if (output)
		    if ((T-tau-last_o)>=Om)
		    {
			output_solution();
			for (int i=0;i<nother;i++)
			    other[i]->output_solution();
			last_o=T-tau;
			printf("T %g\n",T/86400);
		    }
			if (steady_state==1)
				break;
		}
		if (prev_tau!=-1)
			if (prev_tau>tau)
			{ 
				T-=tau;
				tau=prev_tau; 
				T+=tau;
				for (int i=0;i<nother;i++)
				{
					other[i]->T-=other[i]->tau;
					other[i]->tau=prev_tau;
					other[i]->T+=other[i]->tau;
				}
			}
	}
	// constructor
	basic_solver(double _L,int _N,double e)
	{
		N=_N;
		ls_eps=e;
		b_U=new double[(N+2)];
		sb_U=new double[(N+2)];
		MM=new double[(N+2)];
		RP=new double[(N+2)];

		L = _L;
		dL = L/N;
		tau = st_tau;
		T=tau;
	}
	// desctructor
	~basic_solver()
	{
		delete [] b_U;
		delete [] sb_U;
		delete [] MM;
		delete [] RP;
		if (pr_A) delete [] pr_A;
		if (pr_B) delete [] pr_B;
		if (pr_R) delete [] pr_R;
	}
	// output
	virtual void output_solution()=0;
};
class basic_fr_solver: public basic_solver {
public:
	double Dg;
	double alpha;
	double *sp_precalc=NULL;
	double balance;
	virtual double testing_d2fdzdt(double z,double t)=0;
	virtual double testing_d2afdzdt(double z,double t)=0;
	virtual double testing_dfdt(double z,double t)=0;
	virtual double bf(double z,double t)=0;
	std::vector<double> ppv;
	// integrate bf lowering step to achieve less than eps difference
	double integrate_bf(double z,double T)
	{
	    double v=0;
	    double prev_v=1e30;
	    int N=1;
	    int maxN=10000;
#pragma omp critical
	    if (ppv.size()==0)
	    for (int i=0;i<maxN;i++)
		ppv.push_back(pow((double)(i+1),1-alpha)-pow(double(i),1-alpha));
	    while (N<maxN)
	    {
		double dz=z/N;
		double mult=pow(dz,1-alpha);
		v=0;
		for (int i=N-1;i>=0;i--)
		{
		    double m=ppv[N-(i+1)]*mult;
		    if ((alpha==1)&&(i==(N-1))) m=1;
		    v+=bf((i+0.5)*dz,T)*m;
		}
		if (fabs(v-prev_v)<fr_eps)
		    break;
		prev_v=v;
		N*=2;
	    }
	    return v;
	}
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
	    {
		pr_K[i]=Kr(pr_w[i],i)*pr_vgm_k[i];
		if (testing_mode==1) //k=k0*(1-h*h)
		{
		    if (b_U[i]>=0) 
			pr_K[i]=pr_vgm_k[i];
		    else
			pr_K[i]=pr_vgm_k[i]*(1.0-b_U[i]*b_U[i]);
		}
	    }
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
	// testing
	const double tfk=0.01;
	double testing_f(double z,double t) // solution for testing mode
	{
		return tfk+t*t*(z*z-z)/(30.0*86400*30.0*86400);
	}
	double testing_dfdt(double z,double t)
	{
		return  2*t*(z*z-z)/(30.0*86400*30.0*86400);
	}
	double testing_dfdz(double z,double t)
	{
		return  (2*z-1)*t*t/(30.0*86400*30.0*86400);
	}
	double testing_dfadz(double z,double t) // G(2-a)/x^(1-a) D^a_z
	{
		double i1=(Gamma(3)/(2-alpha))*z-1;
		return (t*t/(30.0*86400*30.0*86400))*i1; 
	}
	double testing_d2fdzdt(double z,double t)
	{
		return  2*(2*z-1)*t/(30.0*86400*30.0*86400);
	}
	double testing_d2afdzdt(double z,double t)
	{
		double i1=(Gamma(3)/(2-alpha))*z-1;
		return (2*t/(30.0*86400*30.0*86400))*i1; 
	}
	double testing_d2fdz2(double z,double t)
	{
		return  2*t*t/(30.0*86400*30.0*86400);
	}
	double testing_d2fadz2(double z,double t) // G(2-a)/x^(1-a) D^a_Z G(2-a)/x^(1-a) K D^a_Z H-z   K=1-H*H
	{
	    if (testing_f(z,t)>=0)
	    {
		double i1=(Gamma(3)/(2-alpha));
		return pr_vgm_k[int(z/dL)]*(t*t/(30.0*86400*30.0*86400))*i1; 
	    }
	    else
	    {
		double a=t*t/(30.0*86400*30.0*86400);
		double a1=((1-tfk*tfk)*2*a/(2-alpha))-2*tfk*a*(1+a);
		double a2=(2*tfk*2*a*a/(2-alpha))+(2*tfk*a+a*a)*(1+a);
		double a3=((2*tfk*a+a*a)*2*a/(2-alpha))+2*a*a*(1+a);
		double a4=(4*a*a*a/(2-alpha))+a*a*(1+a);
		double a5=2*a*a*a/(2-alpha);
		double ret= pr_vgm_k[int(z/dL)]*(a1+(2/(2-alpha))*a2*z-(Gamma(4)*Gamma(2-alpha)/Gamma(4-alpha))*a3*z*z+
			(Gamma(5)*Gamma(2-alpha)/Gamma(5-alpha))*a4*z*z*z-(Gamma(6)*Gamma(2-alpha)/Gamma(6-alpha))*a5*z*z*z*z);
		return ret;
	    }
	}
	double bf(double x,double t)
	{
	    return (wetness(0,-testing_f(x,t))*Zs->testing_d2fdzdt(x,t)+dw_dh(-testing_f(x,t),0)*testing_dfdz(x,t)*Zs->testing_dfdt(x,t));
	}
	double testing_F(double z,double t)
	{
		double F;
		F=(dw_dh(-testing_f(z,t),0)+wetness(0,-testing_f(z,t))*pr_vgm_cw[0])*testing_dfdt(z,t);
		F-=testing_d2fadz2(z,t);
		if ((alpha==1)||(z==0))
			F-=B_*bf(z,t);
		else
			F-=B_*pow(z,alpha-1)*integrate_bf(z,t);
		return F;
	}
	// upper flow
	double flow0(double *pickard_prevU)
	{
	    double flow=pr_ET-pr_perc;
	    // irrigation
	    if (is_irr)
		if (irr_on==1)
		    flow-=irr_flow;
	    if (testing_mode==1)
			flow=(pr_K[1]+pr_K[0])*0.5*(testing_dfadz(0,T)-1);
	    return flow;
	}
	// upper boundary condition (DbU=Uc)
	double Uc(double *pickard_prevU)
	{
	    double flow=flow0(pickard_prevU);
	    // condition for flux
	    double kk;
	    if (testing_mode==0)
		kk=(pr_K[1]+Km1(0,pickard_prevU))/2;
	    else
		kk=(pr_K[1]+pr_K[0])/2;
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
		{
		    km1=pr_vgm_k[0]*Kr(wetness(0,-(2*pickard_prevU[0]-pickard_prevU[1])),0);
		    if (testing_mode==1)
			km1=pr_vgm_k[0]*(1-(testing_f(ii*dL,T)-2*Uc(pickard_prevU))*(testing_f(ii*dL,T)-2*Uc(pickard_prevU)));
		}
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
		if (testing_mode==1) ret+=testing_F(ii*dL,T);
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
			if (testing_mode==1)
			    RP[N]=testing_f(N*dL,T);
		}
#pragma omp barrier
//		for (int i=0;i<2;i++)
//			printf("%g %d %g %g %g %g\n",T,i,pr_A[i],pr_B[i],pr_R[i],RP[i]);
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
		{
			b_U[i] = L+pr_fl-(N-i)*dL;
			if (testing_mode==1) b_U[i]=testing_f(i*dL,0);
		}
		out=fopen("out_H.txt","wt");
		v_out=fopen("out_V.txt","wt");
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
		if ((pr_w)&&(testing_mode==0))
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
		// testing
		if (testing_mode==1)
		{
		    double diff=0,sq1=0,sq2=0;
		    for (int i=0;i<N+1;i++)
		    {
			diff+=(b_U[i]-testing_f(i*dL,T-tau))*(b_U[i]-testing_f(i*dL,T-tau));
			sq1+=b_U[i]*b_U[i];
			sq2+=testing_f(i*dL,T-tau)*testing_f(i*dL,T-tau);
		    }
		    diff=sqrt(diff)/N;
		    sq1=sqrt(sq1);
		    sq2=sqrt(sq2);
		    printf ("testing diff %g rel %g at 0 %g %g at N %g %g\n",diff,fabs(sq1-sq2)/sq2,b_U[0],testing_f(0,T-tau),b_U[N],testing_f(N*dL,T-tau));
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
// [2*N,...,3*N-1] - n
// steady state - (l+2m)*d2 z /dz2 = rho g
class Z_solver: public basic_fr_solver {
public:
	// inputs
	double rho_s,lambda,mu,B_;
	H_solver *Hs; // solver for H equation

	FILE *out2;
	FILE *out3;

	int N1,O1; // last index for z, first for v
	int N2,O2; // last index for v, first for n

	double *pr_D2=NULL;
	double *pr_D3=NULL;
	double *pr_D4=NULL;
	double *pr_D5=NULL;

	// testing
	const double mult=10;
	double testing_f(double z,double t) // solution for testing mode
	{
		return mult*t*(2-z*z)/(30.0*86400*30.0*86400);
	}
	double testing_dfdt(double z,double t)
	{
		return mult*(2-z*z)/(30.0*86400*30.0*86400);
	}
	double testing_d2fdt2(double z,double t)
	{
		return  0;
	}
	double testing_dfdz(double z,double t)
	{
		return -2*mult*t*z/(30.0*86400*30.0*86400);
	}
	double testing_dafdz(double z,double t)
	{
		return -(2/(2-alpha))*mult*t*z/(30.0*86400*30.0*86400);
	}
	double testing_d2fdzdt(double z,double t)
	{
		return -2*mult*z/(30.0*86400*30.0*86400);
	}
	double testing_d2afdzdt(double z,double t)
	{
		return -(2/(2-alpha))*mult*z/(30.0*86400*30.0*86400);
	}
	double testing_d2fdz2(double z,double t)
	{
		return -2*mult*t/(30.0*86400*30.0*86400);
	}
	double testing_d2afdz2(double z,double t)
	{
		return -(2/(2-alpha))*mult*t/(30.0*86400*30.0*86400);
	}
	double testing_n(double z,double t) // solution of dn/dt=D_z((1-n)v)
	{
		double a=mult/(30.0*86400*30.0*86400);
		double b=(1-sigma)*2*a;
		return 1-b/testing_dfdt(z,t);
	}
	double testing_dndz(double z,double t)
	{
		double a=mult/(30.0*86400*30.0*86400);
		double b=(1-sigma)*2*a;
		return b*testing_d2fdzdt(z,t)/(testing_dfdt(z,t)*testing_dfdt(z,t));
	}
	double bf(double x,double t)  // d/dz theta*h/n
	{
		return ((Hs->dw_dh(-Hs->testing_f(x,t),0)*Hs->testing_dfdz(x,t)*Hs->testing_f(x,t)/testing_n(x,t))+
				Hs->testing_dfdz(x,t)*(Hs->wetness(0,-Hs->testing_f(x,t))/testing_n(x,t))-
				Hs->wetness(0,-Hs->testing_f(x,t))*Hs->testing_f(x,t)*testing_dndz(x,t)/(testing_n(x,t)*testing_n(x,t)));
	}
	double testing_F2(double z,double t) // for dz/dt
	{
		return ((1-testing_n(z,t))*rho_s+testing_n(z,t)*rho_w)*(testing_d2fdt2(z,t)-g);
	}
	double testing_F3(double z,double t) // for dz/dt
	{
		return -(lambda+2.0*mu)*testing_d2afdz2(z,t);
	}
	double testing_F4(double z,double t) // for dz/dt
	{
		if ((alpha==1)||(z==0))
			return B_*rho_w*g*bf(z,t);
		return B_*rho_w*g*pow(z,alpha-1)*integrate_bf(z,t);
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
		if (i<O2) // v
			return ((1-pickard_prevU[O2+i-O1])*rho_s+pickard_prevU[O2+i-O1]*rho_w)/tau;
		return 1.0/tau; // n
	}
	// right part 
	double Rp(int ii,double *b_Uold,double *pickard_prevU)
	{
		if (steady_state==1) 
				return ((ii<N1)?((1-sigma)*rho_s+sigma*rho_w)*g:((ii>=O2)?sigma:0.0));
		if ((ii<O1)||(ii>=O2)) // z or n
			return b_Uold[ii]/tau;
		// v
		double ret=((1-pickard_prevU[O2+ii-O1])*rho_s+pickard_prevU[O2+ii-O1]*rho_w)*((b_Uold[ii]/tau)+g);
		if (testing_mode==1) 
			ret+=testing_F2((ii-O1)*dL,T);
		return ret;
	}
	// linear solver
	// matrix multiplication on vector UU
	void mmult(double *UU)
	{
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
			MM[N2]=UU[N2];
			// upper boundary condition
			MM[O1]=UU[O1]-UU[O1+1];
			if (steady_state==1)
				MM[O1]=UU[O1];
		    // n
			// upper boundary condition
			if (testing_mode==0)
				MM[O2]=UU[O2]-UU[O2+1];
			else
				MM[O2]=UU[O2];
			// bottom boundary condition
			MM[N]=UU[N];
			
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
		    if (pr_D5==NULL) pr_D5=new double[N+2];
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
			pr_D4[i]=((Hs->pr_w[i+1]/pickard_prevU[O2+i+1])*Hs->b_U[i+1]-(Hs->pr_w[i-1]/pickard_prevU[O2+i-1])*Hs->b_U[i-1])/(2.0*dL);
		    else
			pr_D4[i]=((Hs->pr_w[i+1]/pickard_prevU[O2+i+1])*Hs->b_U[i+1]-(Hs->pr_w[i]/pickard_prevU[O2+i])*Hs->b_U[i])/dL;
		}
		// for n
#pragma omp barrier
#pragma omp for nowait
		for (int i = 0;i <=N1 ;i++)
		{
		    if ((i!=0)&&(i!=N1)) 
			pr_D5[i]= ((1-pickard_prevU[i+O2+1])*pickard_prevU[i+O1+1]-(1-pickard_prevU[i+O2-1])*pickard_prevU[i+O1-1])/(2.0*dL);
		    if (i==0)
			pr_D5[i]= ((1-pickard_prevU[i+O2+1])*pickard_prevU[i+O1+1]-(1-pickard_prevU[i+O2])*pickard_prevU[i+O1])/dL;
		    if (i==N1)
			pr_D5[i]= ((1-pickard_prevU[i+O2])*pickard_prevU[i+O1]-(1-pickard_prevU[i+O2-1])*pickard_prevU[i+O1-1])/dL;
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
		    if (testing_mode==1) 
				v+=testing_F3(i*dL,T);
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
			if (testing_mode==1) 
			    v-=testing_F4(i*dL,T);
			if (steady_state==0)
				RP[i+O1]-=((pr_R[i+O1]!=0)?v/pr_R[i+O1]:v);
			else
				RP[i]-=((pr_R[i]!=0)?v/pr_R[i]:v);
		}
		// n
#pragma omp barrier
		if (steady_state==0)
#pragma omp for nowait
		for (int i = 1;i < N1 ;i++)
		{
			double v=Dg*dL*dL*pow((i+1)*dL,alpha-1)*(SFD(i,pr_D5,alpha,sp_precalc,N1)+pr_D5[i]);
			RP[i+O2]+=((pr_R[i+O2]!=0)?v/pr_R[i+O2]:v);
		}
#pragma omp barrier
#pragma omp single
		{
			// z
			// bottom boundary condition
			RP[N1]=0.0;
			if (testing_mode==1)
			    RP[N1]=testing_f(N1*dL,T);
			// upper boundary condition
			RP[0]=0;
			if (testing_mode==1)
			    RP[0]=-dL*testing_dfdz(0,T);
			if (steady_state==1) RP[0]=0;
			// v
			// bottom boundary condition
			RP[N2]=0.0;
			if (testing_mode==1)
			    RP[N2]=testing_dfdt(N1*dL,T);
			// upper boundary condition
			RP[O1]=-dL*testing_d2fdzdt(0,T);
			if (steady_state==1)
				RP[O1]=0;
			// n
			RP[O2]=0.0;
			RP[N]=sigma;
			if (testing_mode==1)
			{
				RP[O2]=testing_n(0,T);
				RP[N]=testing_n(N1*dL,T);
			}
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
	Z_solver(int _N,double e,double _L,H_solver *_Vs,double _r,double _l,double _m,double _a,double _b) : basic_fr_solver(_L,_N,e), Hs(_Vs),rho_s(_r),lambda(_l),mu(_m),B_(_b)
	{
		out=fopen("out_Z.txt","wt");
		out2=fopen("out_Vz.txt","wt");
		out3=fopen("out_n.txt","wt");
		alpha=_a;
		N1=(N-2)/3;
		O1=N1+1;
		N2=O1+N1;
		O2=N2+1;
		dL = L/N1;
		Dg=pow(dL,-1-alpha);
		for (int i = 0;i < N + 1;i++)
		{
			b_U[i] = 0;
			if (testing_mode==1)
			{
			    if (i<=N1)
				 b_U[i]=testing_f(i*dL,0);
			    else
				if (i<=N2)
				 b_U[i]=testing_dfdt((i-O1)*dL,0);
			}
			if (i>=O2)
			    b_U[i]=((testing_mode==1)?testing_n((i-O2)*dL,0):sigma);
		}
	}
	// destructor
	~Z_solver()
	{
		if (pr_D2) delete [] pr_D2;
		if (pr_D3) delete [] pr_D3;
		if (pr_D4) delete [] pr_D4;
		if (pr_D5) delete [] pr_D5;
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
		    fprintf(out3,"%g ",b_U[i+O2]);
		}
		// balance of momentum
		if ((Hs->pr_w!=NULL)&&(testing_mode==0))
		{
		    double total_momentum=0,total_rhog=0;
		    double mult=1.0/Gamma(1+alpha);
		    double mult2=1.0/Gamma(2-alpha);
		    for (int i=0;i<N1;i++)
		    {
			double b=pow((N1-i)*dL,alpha)-pow((N1-(i+1))*dL,alpha);
			total_momentum+=0.5*(((1-b_U[O2+i])*rho_s+b_U[O2+i]*rho_w)*b_U[O1+i]+((1-b_U[O2+i+1])*rho_s+b_U[O2+i+1]*rho_w)*b_U[O1+i+1])*
					 b*mult2*pow((i+1)*dL,1-alpha);
			total_rhog+=0.5*(((1-b_U[O2+i])*rho_s+b_U[O2+i]*rho_w)*g+((1-b_U[O2+i+1])*rho_s+b_U[O2+i+1]*rho_w)*g)*b*mult2*pow((i+1)*dL,1-alpha);
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
		    pr_D4[0]=((Hs->pr_w[1]/b_U[O2+1])*Hs->b_U[1]-(Hs->pr_w[0]/b_U[O2+0])*Hs->b_U[0])/dL;
		    pr_D4[N1]=((Hs->pr_w[N1]/b_U[O2+N1])*Hs->b_U[N1]-(Hs->pr_w[N1-1]/b_U[O2+N1-1])*Hs->b_U[N1-1])/dL;
		    double flux2=-B_*rho_w*g*pr_D4[0];
		    double flux2a=-B_*rho_w*g*pr_D4[N1];
		    if (balance==1e30)
			 balance=total_momentum;
		    fprintf(out,"total_momentum %g balance %g flux_sigma %g %g fluxP %g %g total_rhog %g ",total_momentum,balance,flux1,flux1a,flux2,flux2a,total_rhog);
		    balance+=(flux1+flux1a+flux2+flux2a+total_rhog)*(T-tau-last_o);
		}
		if (testing_mode==1)
		{
		    double diff=0,sq1=0,sq2=0;
		    for (int i=0;i<=N1;i++)
		    {
			diff+=(b_U[i]-testing_f(i*dL,T-tau))*(b_U[i]-testing_f(i*dL,T-tau));
			sq1+=b_U[i]*b_U[i];
			sq2+=testing_f(i*dL,T-tau)*testing_f(i*dL,T-tau);
		    }
		    diff=sqrt(diff)/N1;
		    sq1=sqrt(sq1);
		    sq2=sqrt(sq2);
		    printf ("testing z diff %g rel %g at 0 %g %g at N %g %g\n",diff,fabs(sq1-sq2)/sq2,b_U[0],testing_f(0,T-tau),b_U[N1],testing_f(N1*dL,T-tau));
		    diff=sq1=sq2=0;
		    for (int i=0;i<=N1;i++)
		    {
			diff+=(b_U[i+O1]-testing_dfdt(i*dL,T-tau))*(b_U[i+O1]-testing_dfdt(i*dL,T-tau));
			sq1+=b_U[i+O1]*b_U[i+O1];
			sq2+=testing_dfdt(i*dL,T-tau)*testing_dfdt(i*dL,T-tau);
		    }
		    diff=sqrt(diff)/N1;
		    sq1=sqrt(sq1);
		    sq2=sqrt(sq2);
		    printf ("testing v diff %g rel %g at 0 %g %g at N %g %g\n",diff,fabs(sq1-sq2)/sq2,b_U[O1],testing_dfdt(0,T-tau),b_U[N2],testing_dfdt(N1*dL,T-tau));
		    diff=sq1=sq2=0;
		    for (int i=0;i<=N1;i++)
		    {
			diff+=(b_U[i+O2]-testing_n(i*dL,T-tau))*(b_U[i+O2]-testing_n(i*dL,T-tau));
			sq1+=b_U[i+O2]*b_U[i+O2];
			sq2+=testing_n(i*dL,T-tau)*testing_n(i*dL,T-tau);
		    }
		    diff=sqrt(diff)/N1;
		    sq1=sqrt(sq1);
		    sq2=sqrt(sq2);
		    printf ("testing n diff %g rel %g at 0 %g %g at N %g %g\n",diff,fabs(sq1-sq2)/sq2,b_U[O2],testing_n(0,T-tau),b_U[N],testing_n(N1*dL,T-tau));
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
		if (strcmp(argv[i],"testing_mode")==0)
			testing_mode=atoi(argv[i+1]);
		// (v,z) equation parameters
		if (strcmp(argv[i],"rho_s")==0)
			rho_s=atof(argv[i+1]);
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
	Z_solver *Zsolver=new Z_solver(3*N+2,ls_eps2,L,solver,rho_s,lambda,mu,alpha,B);
	solver->Zs=Zsolver;
	basic_solver *others[1]={Zsolver};
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
	solver->soil_calc_coefs();
	// solve steady state problem for H
	if (testing_mode==0)
	{
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
	}
	// run simulation
	solver->output_solution();
	Zsolver->output_solution();

	solver->solve_up_to_T(Tm,solver->tau,others,1,1);

	solver->output_solution();
	Zsolver->output_solution();
	return 0;
}

