#include <fstream>
#include <cmath>
#include <cmath>


using namespace std;

void f (double *  y0, const double x);
void RKstep(double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,double *k7,  double*  y,  double x, double dx);
void RK4(double *k1, double *k2, double *k3, double *k4, double *k5, double *k6, double *k7,double *y4, double x,  double dx);
void RK5(double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,double *k7, double *y5,  double x,  double dx);
void max_bestimmen (double *yn2, double *yn, double &maximum);

int main()
{
	double x=0.0, dx=1e-3, tol=1e-5;
	const double L=17.065216560157;
	const int dim=4;
	double k1[dim], k2[dim], k3[dim], k4[dim], k5[dim], k6[dim], k7[dim];
	double  y0[dim]= { 0.994, 0.0 , 0.0 , -2.00158510637908};
	double  y02[dim]= { 0.994, 0.0 , 0.0 , -2.00158510637908};
	double  maximum;
	double y_alt[dim], dx_new;
	
	
	ofstream out("Runge_kutta.txt");
	while(x<L)
	{
		for(int i=0; i<dim; i++) y_alt[i]=y0[i];
		RK4(k1,k2,k3, k4,k5, k6, k7,y0, x,  dx);
		RK5(k1, k2, k3, k4, k5, k6,k7,y02, x,  dx);

		
		max_bestimmen(y0,y02, maximum);
		for(int i=0; i<dim; i++)
			y0[i]=y_alt[i];
		dx_new=dx*pow((tol/maximum), 1./5.);
		dx=dx_new;
		RK4(k1,k2,k3, k4,k5, k6, k7, y0, x,  dx);
		for(int i=0; i<dim; i++)
			y02[i]=y0[i];
		x+=dx;
		out<<x<<"\t "<<y0[0]<<"\t "<<y0[1]<<"\t "<<dx<<endl;

		
	}

	
    getche();
    return 0;
}
void f (double *  y0, double x)
{
	double const mu=0.012277471;
	double y[4]={y0[0],y0[1], y0[2], y0[3]}; // x , y, x_ableitung, y_ableitung
	
	
	const double r=sqrt(pow((y[0]+mu),2)+pow(y[1],2));
	const double s=sqrt(pow((y[0]-1+mu),2)+pow(y[1],2));

	y0[0]=y[2];
	y0[1]=y[3];
	y0[2]=y[0]+2*y[3]-((1-mu)*(y[0]+mu))/(pow(r,3))-mu*(y[0]-1+mu)/(pow(s,3));
	y0[3]=y[1]-2*y[2]-(y[1]*(1-mu)/(pow(r,3)))-mu*y[1]/(pow(s,3));
	
}

void RKstep(double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,double *k7,  double*  y,  double x, double dx)
{
	const int dim=4;
	
	for(int i=0; i<dim; i++)
		k1[i]=y[i];
	f(k1, x);
	for(int i=0; i<dim; i++)
		k2[i]=y[i]+1./5.*dx*k1[i];
	f(k2,x+dx*1./5.);
	for(int i=0; i<dim; i++)
		k3[i]=y[i]+3./40.*k1[i]*dx+9./40.*k2[i]*dx;
	f(k3,x+(3.0/10.0)*dx);
	for(int i=0; i<dim; i++)
		k4[i]=y[i]+44./45.*dx*k1[i]-56./15.*dx*k2[i]+32./9.*k3[i]*dx;
	f(k4,x+4./5.*dx);
	for(int i=0; i<dim; i++)
		k5[i]=y[i]+19372./6561.*dx*k1[i]-25360./2187.*dx*k2[i]+64448./6561.*dx*k3[i]-212./729.*k4[i]*dx;
	f(k5,x+8./9.*dx);
	for(int i=0; i<dim; i++)
		k6[i]=y[i]+dx*(9017./3168.*k1[i]-355./33.*k2[i]+46732./5247.*k3[i]+49./176.*k4[i]-5103./18656.*k5[i]);
	f(k6,x+dx);

	

	for(int i=0; i<dim; i++)
		k7[i]=y[i]+dx*(35./384.*k1[i]+500./1113*k3[i]+125./192*k4[i]-2187./6784.*k5[i]+11./84.*k6[i]);
	f(k7,x+dx);

	
}

void RK4( double *k1, double *k2, double *k3, double *k4, double *k5, double *k6, double *k7,double *y4, double x,  double dx)
{
	const int dim=4;
	RKstep(k1, k2, k3,k4, k5, k6, k7,  y4, x,dx);
	for(int i=0; i<dim; i++)
		y4[i]+=dx*(35./384.*k1[i]+500./1113*k3[i]+125./192*k4[i]-2187./6784.*k5[i]+11./84.*k6[i]);
}

void RK5(double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,double *k7, double *y5,  double x,  double dx)
{
	const int dim=4;
	RKstep(k1, k2, k3,k4, k5, k6, k7,  y5, x,dx);
	for(int i=0; i<dim; i++)
		y5[i]+=+dx*(5179./57600*k1[i]+7571./16695.*k3[i]+393./640.*k4[i]-92097./339200.*k5[i]+187./2100.*k6[i]+1./40.*k7[i] );

}
void max_bestimmen (double *yn2, double *yn, double &maximum)
{
	
	const int dim=4;
	double norm[dim];
	for(int i=0; i<dim; i++)
		norm[i]=abs(yn2[i]-yn[i]);

    maximum=norm[0];
	for(int i=0; i<dim; i++)
	{
		if(norm[i]>maximum)
			maximum=norm[i];
	}
	
}
