// Model from https://doi.org/10.1007/978-981-10-8773-8_1
// numerical scheme roughly from https://doi.org/10.1137/0913034
// d/dt(theta*c+s_a+s_s)=-duc/dx
// ds_a/dt=l_a*c*u-k_det*s_a
// ds_s/dt=l_s*c*u
// all concetrations >0
// c - concentration of suspended particles in water
// s_a - of attached particles
// s_s - of strained particles
// u - advection velocity
// theta - moisture content
// hyperbolic equation solver
//
// get vgm for (sand,silt,clay) from Rosetta python module
#include <Python.h>
class particle_transport:public basic_solver
{
public:
    double l_a,l_s,k_det; // attachment coef, straining coef, detachment coef
    double *b_v; // advection velocity
    double *theta; // moisture content
    double ds_s(int i, double *U,double *pickard_prevU) // ds_s/dt
    {
	double dss=l_s*pickard_prevU[i]*b_v[i];
	if ((U[2*(N/3)+i]+tau*dss)<0) dss=0;
	return dss;
    }
    double ds_a(int i, double *U,double *pickard_prevU) // ds_a/dt
    {
	double dsa=l_a*pickard_prevU[i]*b_v[i]-k_det*pickard_prevU[(N/3)+i];
	if ((U[(N/3)+i]+tau*dsa)<0) dsa=0;
	return dsa;
    }
    double dc(int i, double *U,double *pickard_prevU) // d(theta c)/dt
    {
	double dss=ds_s(i,U,pickard_prevU);
	double dsa=ds_a(i,U,pickard_prevU);
	double dcx=0; // d uc /dx
	double fi=b_v[i]*pickard_prevU[i];
	if (i==0) // boundary - uc[0]=0
	{
	    if (b_v[0]>0)
		dcx=fi;
	    else
		dcx=b_v[1]*pickard_prevU[1];
	}
	else
	{
	    if (i==((N/3)-1)) // boundary - uc[N]=0
	    {
		if (b_v[(N/3)-1]<=0)
		    dcx=-fi;
		else
		    dcx=-b_v[(N/3)-2]*pickard_prevU[(N/3)-2];
	    }
	    else
	    {
		double fl=b_v[i-1]*pickard_prevU[i-1];
		double fr=b_v[i+1]*pickard_prevU[i+1];
		if (b_v[i]>0)
		{
		    if (b_v[i-1]<=0) // shock/sonic points
			fl=(1.0/2.0)*(b_v[i]*pickard_prevU[i]+b_v[i-1]*pickard_prevU[i-1]);
		    if (b_v[i+1]<=0)
			fi=(1.0/2.0)*(b_v[i+1]*pickard_prevU[i+1]+b_v[i]*pickard_prevU[i]);
		    dcx=fi-fl; // general formula
		}
		else
		{
		    if (b_v[i+1]>0) // shock/sonic points
			fr=(1.0/2.0)*(b_v[i+1]*pickard_prevU[i+1]+b_v[i]*pickard_prevU[i]);
		    if (b_v[i-1]>0)
			fi=(1.0/2.0)*(b_v[i]*pickard_prevU[i]+b_v[i-1]*pickard_prevU[i-1]);
		    dcx=fr-fi; // general formula
		}
	    }
	}
	return -dss-dsa-dcx/dL;
    }
    // three-diagonal matrix
    double A_(int i, double *pickard_prevU)
    {
	return 0;
    }
    double B(int i, double *pickard_prevU)
    {
	return 0;
    }
    double R(int i, double *pickard_prevU)
    {
	return 1.0/tau;
    }
    double Rp(int i, double *U,double *pickard_prevU)
    {
	double ret=U[i]/tau;
	if ((i>=0)&&(i<N/3)) ret+=dc(i,U,pickard_prevU)/theta[i];
	if ((i>=N/3)&&(i<2*(N/3))) ret+=ds_a(i-N/3,U,pickard_prevU);
	if ((i>=2*(N/3))&&(i<N)) ret+=ds_s(i-2*N/3,U,pickard_prevU);
	return ret;
    }
    void mmult(double *UU)
    {
	if ((l_a==0)&&(l_s==0)&&(k_det==0)) return;
	mmult_main(UU);
#pragma omp barrier
#pragma omp single
	{
	    MM[0]=UU[0];
	    MM[N]=UU[N];
	}
#pragma omp barrier
    }
    void __precalc(double *pickard_prevU)
    {
	if ((l_a==0)&&(l_s==0)&&(k_det==0)) return;
	__precalc_main(pickard_prevU);
#pragma omp barrier
    }
    particle_transport(int _N, double _dL,double _tau,double _l_a,double _l_s,double _k_det,double *_b_v,double *_theta,double c0,int idx,double eps):
	    basic_solver(_dL*((_N/3)-1),_N,eps),l_a(_l_a),l_s(_l_s),k_det(_k_det),b_v(_b_v),theta(_theta)
    {
	char str[1024];
	sprintf(str,"out_P_%d.txt",idx);
	out=fopen(str,"wt");
	dL*=3;
	tau=_tau;
	pr_A=new double[N];
	pr_B=new double[N];
	pr_R=new double[N];
	for (int i=0;i<(N/3);i++)
	{
	    b_U[i]=0;
	    b_U[(N/3)+i]=c0;
	    b_U[2*(N/3)+i]=0;
	}
	printf("particle transport %d: l_a %g l_s %g k_det %g c0 %g\n",idx,l_a,l_s,k_det,c0);
    }
    ~particle_transport()
    {
	fclose(out);
    }
    void output_solution()
    {
	if (T<=tau)
	{
	    fprintf(out,"t(days) %g tau(seconds) %g - sumC: ",(T-tau)/86400,tau);
	    for (int i=0;i<(N/3);i++)
		fprintf(out,"%g ",i*dL);
	    fprintf(out,"c: ");
	    for (int i=0;i<(N/3);i++)
		fprintf(out,"%g ",i*dL);
	    fprintf(out,"s_a: ");
	    for (int i=0;i<(N/3);i++)
		fprintf(out,"%g ",i*dL);
	    fprintf(out,"s_s: ");
	    for (int i=0;i<(N/3);i++)
		fprintf(out,"%g ",i*dL);
	    fprintf(out,"\n");
	}
	fprintf(out,"t(days) %g tau(seconds) %g - sumC: ",(T-tau)/86400,tau);
	for (int i=0;i<(N/3);i++)
		fprintf(out,"%g ",theta[i]*b_U[i]+b_U[(N/3)+i]+b_U[2*(N/3)+i]);
	fprintf(out,"c: ");
	for (int i=0;i<(N/3);i++)
		fprintf(out,"%g ",b_U[i]);
	fprintf(out,"s_a: ");
	for (int i=0;i<(N/3);i++)
		fprintf(out,"%g ",b_U[(N/3)+i]);
	fprintf(out,"s_s: ");
	for (int i=0;i<(N/3);i++)
		fprintf(out,"%g ",b_U[2*(N/3)+i]);
	double sum=0;
	for (int i=0;i<(N/3);i++)
		sum+=dL*(theta[i]*b_U[i]+b_U[(N/3)+i]+b_U[2*(N/3)+i]);
	fprintf(out," totalC %g\n",sum);
    }
};
// call Rosetta python module to get VGM parameters for chaning (sand,silt,clay)
const char *p_code_import="from rosetta import rosetta, SoilData";
const char *p_code_run="data = [%s] \n\
mean, stdev, codes = rosetta(3, SoilData.from_array(data))\n\
mean=mean.tolist()";
char *data_string,*data_string2;
double *vgm_mean;
PyObject *mainModule;
PyObject *var1Py;
int *prev_pt_b_U=NULL;
int *cur_pt_b_U=NULL;
FILE *vgm_out;
void python_initialize(int N)
{
	Py_Initialize();
	PyRun_SimpleString(p_code_import);
	data_string=new char[N*100+500];
	data_string2=new char[N*100+500];
	vgm_mean=new double[5*N];
	mainModule = PyImport_AddModule("__main__");
	vgm_out=fopen("out_vgm.txt","wt");
}
void python_get_mean_vgm(particle_transport **pt,int N,double T)
{
	if (cur_pt_b_U==NULL)
	    cur_pt_b_U=new int[3*N];
	// form a string with (sand, silt, clay) for each cell
	data_string[0]=0;
	int first=1;
	for (int i=0;i<N;i++)
	{
		// calc
		double v0=pt[0]->theta[i]*pt[0]->b_U[i]+pt[0]->b_U[N+i]+pt[0]->b_U[2*N+i];
		double v1=pt[0]->theta[i]*pt[1]->b_U[i]+pt[1]->b_U[N+i]+pt[1]->b_U[2*N+i];
		double v2=pt[0]->theta[i]*pt[2]->b_U[i]+pt[2]->b_U[N+i]+pt[2]->b_U[2*N+i];
		int i0=int(100*v0);
		int i1=int(100*v1);
		int i2=int(100*v2);
		cur_pt_b_U[3*i+0]=i0;
		cur_pt_b_U[3*i+1]=i1;
		cur_pt_b_U[3*i+2]=i2;
		// check for changes
		int add=1;
		if (prev_pt_b_U)
		    if (abs(i0-prev_pt_b_U[3*i+0])==0)
		    if (abs(i1-prev_pt_b_U[3*i+1])==0)
		    if (abs(i2-prev_pt_b_U[3*i+2])==0)
			add=0;
		// add to string
		if (add)
		{
		    i0=int(i0/(v0+v1+v2));
		    i1=int(i1/(v0+v1+v2));
		    i2=int(i2/(v0+v1+v2));
		    if ((i0+i1+i2)!=100) i0+=100-i0-i1-i2;
		    if (first==0) strcat(data_string,",");
		    sprintf(data_string2,"[%d,%d,%d]",i0,i1,i2);
		    strcat(data_string,data_string2);
		    first=0;
		}
	}
	if (first)
	    return;
	sprintf(data_string2,p_code_run, data_string);
	// run python code
	PyRun_SimpleString(data_string2);
	// get back results
	var1Py = PyObject_GetAttrString(mainModule, "mean");
	int idx=0;
	for (int i=0;i<N;i++)
	{
		// check
		int get=1;
		if (prev_pt_b_U)
		    if (abs(cur_pt_b_U[3*i+0]-prev_pt_b_U[3*i+0])==0)
		    if (abs(cur_pt_b_U[3*i+1]-prev_pt_b_U[3*i+1])==0)
		    if (abs(cur_pt_b_U[3*i+2]-prev_pt_b_U[3*i+2])==0)
			get=0;
		if (get==0) continue;
		PyObject *l=PyList_GetItem(var1Py,idx++);
		for (int j=0;j<5;j++)
		{
		    PyObject *ll=PyList_GetItem(l,j);
		    vgm_mean[i*5+j]=PyFloat_AsDouble(ll);
		    if (j>=2) vgm_mean[i*5+j]=pow(10.0,vgm_mean[i*5+j]);
		    if (j==4) vgm_mean[i*5+j]/=100*86400; // cm/day to m/s
		}
		fprintf(vgm_out,"T %g cell %d new t_r %g t_s %g a %g n %g ks %g cur (%d %d %d) prev (%d %d %d)\n",T,i,
		    vgm_mean[i*5+0],vgm_mean[i*5+1],vgm_mean[i*5+2],vgm_mean[i*5+3],vgm_mean[i*5+4],
		    cur_pt_b_U[3*i+0],cur_pt_b_U[3*i+1],cur_pt_b_U[3*i+2],
		    (prev_pt_b_U?prev_pt_b_U[3*i+0]:0),(prev_pt_b_U?prev_pt_b_U[3*i+1]:0),(prev_pt_b_U?prev_pt_b_U[3*i+2]:0));
	}
	// save previous
	if (prev_pt_b_U==NULL)
	    prev_pt_b_U=new int[3*N];
	for (int i=0;i<3*N;i++)
	    prev_pt_b_U[i]=cur_pt_b_U[i];
}
void python_finalize()
{
    delete [] data_string;
    delete [] data_string2;
    delete [] vgm_mean;
    fclose(vgm_out);
    Py_Finalize();
}
