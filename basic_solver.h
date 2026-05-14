// VGM parameters
int vgm_nlayers; // number of layers
double *vgm_ns,*vgm_s0s,*vgm_s1s,*vgm_as,*vgm_h0,*vgm_h1,*vgm_k,*vgm_power,*vgm_cw; // VGM coefs per layer of soil
double *soil_lambda,*soil_Ct,*soil_D[4]; // to read from input file, not used later
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

double Om; // output interval
double last_o=0; // last output time

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
