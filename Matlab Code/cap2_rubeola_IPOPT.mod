# Bioeconomic approach
#
# Jun/2012
#
# rubeola_IPOPT

#### PARAMETERS ###

param b=0.012;
param e=36.5;
param g=30.417;
param p=0.65;
param q=0.65;
param beta=527.59;
param A=100;




param N=999;					# steps number of Euler Scheme 
param tf=3;


param x1_0=0.0555;					# intial conditions
param x2_0=0.0003;
param x3_0=0.0004;
param x4_0=1;
param fc_0=0;

					 

#### VARIABLES ###
var x1 {0..N}, >=0;		# 
var x2 {0..N}, >=0;
var x3 {0..N}, >=0;		# 
var x4 {0..N}, >=0;
var fc {0..N}, >=0;		# functional

var u {0..N},	>=0 <=0.9;						# control variable 




#### OBJECTIVE FUNCTION ###
minimize cost: fc[N];

#### CONSTRAINTS ###
subject to
	i1: x1[0] = x1_0;
	i2: x2[0] = x2_0;
	i3: x3[0] = x3_0;
	i4: x4[0] = x4_0;
	i5: fc[0] = fc_0;
	
	f1 {i in 0..N-1}: x1[i+1] = x1[i] + (tf/N)*(b-b*(p*x2[i]+q*x3[i])-b*x1[i]-beta*x1[i]*x3[i]-u[i]*x1[i]);
	f2 {i in 0..N-1}: x2[i+1] = x2[i] + (tf/N)*(b*p*x2[i]+beta*x1[i]*x3[i]-(e+b)*x2[i]);
	f3 {i in 0..N-1}: x3[i+1] = x3[i] + (tf/N)*(e*x2[i]-(g+b)*x3[i]);
	f4 {i in 0..N-1}: x4[i+1] = x4[i] + (tf/N)*(b-b*x4[i]);
	f5 {i in 0..N-1}: fc[i+1] = fc[i] + (tf/N)*(A*x3[i]+u[i]^2);	


	
#### DATA ###



### SOLVE ###
#option solver ipopt;
#option ipopt_options "max_iter=10000";
solve;

### OUTPUT SCREEN ###
printf: " # cost = %24.16e\n", cost;
printf: " # N = %d\n", N;
printf: "Data\n";
printf{i in 0..N}: "%d %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n", 
i, i*tf/N, x1[i], x2[i], x3[i], x4[i], fc[i], u[i];
end;


