/*
2021/7/20 debug end

Network simulation of Kuramoto model

*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <time.h>

#include "math.c"		//random number (0~1) generation routine
#include "time.c"		//return date and time

//output folder
#define		OUT_FOLDER			"out"

//data label
#define		DLB_MAIN_BASE		"PRG"

#define 	PAI					(atan(1.e0) * 4.e0)
#define		EPS3				1.e-3
#define		EPS6				1.e-6
#define		EPS9				1.e-9
#define		EPS12				1.e-12
#define		BIG30				1.e30
#define		NBG3				1000
#define		NBG6				1000000
#define		NBG9				1000000000
#define		NSA					50
#define		NSB					100
#define		NSC					500
#define		NSD					1000
#define		VER					-2.e6
#define		NER					-20000

//oscillator model
#define		NUM_PHI				10					//number of oscillator
#define		NUM_X_ALL			(NUM_PHI)

#define		MAX_FILE_OPEN		500

#define		RANSUU_METHOD		1					//0: ran_fast_RCP
													//1: ran1_RCP
													//2: ran2_RCP
										
//calculation mode
#define		CALMODE			0						//network simulation(only once)

#define		MAX(a,b)		((a) >= (b) ? (a) : (b))
#define		MIN(a,b)		((a) <= (b) ? (a) : (b))
#define		EQ(a,b)			(fabs((a)-(b))<EPS6 ? 1 : 0)
#define		NE(a,b)			(fabs((a)-(b))>EPS6 ? 1 : 0)
#define		LE(a,b)			((a) < (b)+EPS6 ? 1 : 0)
#define		LT(a,b)			((a) < (b)-EPS6 ? 1 : 0)
#define		GE(a,b)			((a) > (b)-EPS6 ? 1 : 0)
#define		GT(a,b)			((a) > (b)+EPS6 ? 1 : 0)
#define		EQS(s1,s2)		(!strcmp((s1),(s2)))

void check_para_main();
void disp_start_main();
void disp_end_main();
void main_1kai();
void init_M1K();
void output_NW(char c[]);
void rand_NW();
void add_con(long p, long q);
void sim_NW();
void mod_phi(double phi_v[]);
void system_var(long ns, double t, double phi_v[], double *r_order);
void c_x_dt(double t, double phi_v[], double phi_dt_v[]);
void init_SNW(double phi_v[]);
void output_SNW(long ns, double t, double phi_v[], double r_order);
double c_r_order(double phi_v[]);
void output_jikei(long ns, double t, double phi_v[], double r_order);
void output_jikei_sub(long ns, long f_init, double t, double phi_v[], double r_order);
void ruku(double t, double x[], double h);
void my_fopen(FILE **fp, char fname[], long n_label);
void set_dlb_main(int argc, char *argv[]);
float ransuu(void);
void ransuu_init(long i_seed);
long n_ransuu(long n);
void disp_start_end_time(char md[], char c[]);
void disp_err(char c[]);

double T_step, T_init, T_fin;
double Sig_nw, P_rand_nw;

long Fout_jikei;
double T_jikei_min, T_jikei_max;
long Nstep_jikei;
long Id_phi_sel[200];

double T_order_ave_min, T_order_ave_max;

long Ransuu_seed; 

double A_nw[NUM_PHI][NUM_PHI];		//connection matrix of the oscillator network
double Phi_init_v[NUM_PHI];			//initial value of phase
double Omega_v[NUM_PHI];			//natural angular frequency

double R_order_ave;					//time average value of order parameter

char Dlb_main[200];
char Time_text[200];

long I_ransuu;

FILE *Fp_jikei;

int main(int argc, char *argv[]){
	//time
	T_step 		= 0.5e0;
	T_init		= -5.e0;
	T_fin		= 5000.e0;

	Sig_nw		= 0.25e0;
	P_rand_nw	= 0.5e0;
	
	Fout_jikei  	= 1;
	T_jikei_min 	= T_init;
	T_jikei_max 	= T_fin;
	Nstep_jikei 	= 20;
	Id_phi_sel[0] = 0; Id_phi_sel[1] = NUM_PHI/2; Id_phi_sel[2] = NUM_PHI-1; Id_phi_sel[3] = NER;

	//order parameter
	T_order_ave_min		= 0.e0;
	T_order_ave_max		= T_fin;
		
	Ransuu_seed  = -1;
	
    set_dlb_main(argc, argv);
	check_para_main();
	
	disp_start_main();
    disp_start_end_time("s", Dlb_main);
	
	//network simulation
	main_1kai();
	
	disp_start_end_time("e", Dlb_main);
	disp_end_main();
}

//parameter check
void check_para_main(){
	if (Fout_jikei == 1){
		if (!(GE(T_jikei_min, T_init) && LE(T_jikei_max, T_fin))) disp_err("err-CPM-TJKa");
	}
	
	if (!(GE(T_order_ave_min, T_init) && LE(T_order_ave_max, T_fin))) disp_err("err-COM-TOMa");
}

//display output
void disp_start_main(){
	printf("--start of %s--\n", Dlb_main);
}

void disp_end_main(){
	printf("--end of %s--\n", Dlb_main);
}

//network simulation
void main_1kai(){
	long i;
	init_M1K();
	sim_NW();
}

//initialization
void init_M1K(){
	long i;
	ransuu_init(Ransuu_seed); 
	rand_NW();
	output_NW("M1K");
	
	//Phi_init_v, Omega_v
	for (i = 0; i < NUM_PHI; i++) Phi_init_v[i] = ransuu() * 2.e0 * PAI;
	for (i = 0; i < NUM_PHI; i++) Omega_v[i] = ransuu() * 2.e0 - 1.e0;
}


//output network
void output_NW(char c[]){
	long i, j;
	FILE *fp1;
	char text1[200];
	sprintf(text1, "network(%s)", c);
	
	my_fopen(&fp1, text1, 1);
	fprintf(fp1, "ID_of_oscillator(%s)%s, ID_of_oscillator(%s)%s\n", c, Dlb_main, c, Dlb_main);
	
	for (i = 0; i < NUM_PHI; i++) for (j = 0; j < NUM_PHI; j++){
		if (A_nw[i][j] == 1) fprintf(fp1, "%ld, %ld\n", i, j);
	}
	
	fclose(fp1);
}

void rand_NW(){
	long i, j;
	for (i = 0; i < NUM_PHI; i++) for (j = 0; j < NUM_PHI; j++) A_nw[i][j] = 0.e0;
	
	for (i = 0; i < NUM_PHI; i++){
		for (j = i+1; j < NUM_PHI; j++){
			if (ransuu() < P_rand_nw) add_con(i,j);
		}
	}
}

void add_con(long p, long q){
	if (!(p >= 0 && p < NUM_PHI)) disp_err("err: add-con-a");
	if (!(q >= 0 && q < NUM_PHI)) disp_err("err: add-con-b");
	if (p == q) disp_err("err: add-con-c");

	A_nw[p][q] = 1.e0;
	A_nw[q][p] = 1.e0;
}

//network simulation
void sim_NW(){
	double t;
	double phi_v[NUM_PHI];
	long ns;
	double r_order;
	
	init_SNW(phi_v);
	
	for (ns = 0; ;ns++){
		t = T_init + T_step * ns;
		ruku(t, phi_v, T_step);
		mod_phi(phi_v);
		
		system_var(ns, t, phi_v, &r_order);
		output_SNW(ns, t, phi_v, r_order);
		if (GE(t, T_fin)) break;
	}
	printf("R_order_ave = %.7lf (end of SNW)\n", R_order_ave);
}

void mod_phi(double phi_v[]){
	long i;
	for (i = 0; i < NUM_PHI; i++){
		if (phi_v[i] < 0.e0) phi_v[i] += 2.e0 * PAI;
		if (phi_v[i] >= 2.e0 * PAI) phi_v[i] -= 2.e0 * PAI;
	}
}

//r_order, R_order_ave
void system_var(long ns, double t, double phi_v[], double *r_order){
	static double sum;
	if (ns == 0) sum = 0.e0;
	
	//r_order
	*r_order = c_r_order(phi_v);
	
	//R_order_ave
	if (GE(t, T_order_ave_min)&&LE(t, T_order_ave_max - T_step)){
		sum += (*r_order) * T_step;
		if (EQ(t, T_order_ave_max - T_step)) R_order_ave = sum / (T_order_ave_max - T_order_ave_min);
	}
}

//phi_v, phi_dt_v
void c_x_dt(double t, double phi_v[], double phi_dt_v[]){
	long i, j;
	double sum;
	for (i = 0; i < NUM_PHI; i++){
		sum = 0.e0;
		for (j = 0; j < NUM_PHI; j++) sum += A_nw[i][j] * sin(phi_v[j] - phi_v[i]);
		phi_dt_v[i] = Omega_v[i] + Sig_nw * sum;
	}
}

//initialization
void init_SNW(double phi_v[]){
	long i;
	for (i = 0; i < NUM_PHI; i++) phi_v[i] = Phi_init_v[i]; 
}

//output SNW
void output_SNW(long ns, double t, double phi_v[], double r_order){
	if (Fout_jikei == 1) output_jikei(ns, t, phi_v, r_order);
}

//order parameter
double c_r_order(double phi_v[]){
	long i;
	double wk1;
	double sum_re, sum_im;
	double re, im;
	
	//sum_re
	sum_re = 0.e0;
	for (i = 0; i < NUM_PHI; i++) sum_re += cos(phi_v[i]);
	
	//sum_im
	sum_im = 0.e0;
	for (i = 0; i < NUM_PHI; i++) sum_im += sin(phi_v[i]);
	
	//re, im
	re = sum_re / NUM_PHI;
	im = sum_im / NUM_PHI;
	
	wk1 = sqrt(re * re + im * im);
	return wk1;
}

//output of time series data
void output_jikei(long ns, double t, double phi_v[], double r_order){
	static long f_init;
	
	if (ns == 0) f_init = 0;
	
	if (GE(t, T_jikei_min)&&LE(t, T_jikei_max)){
	
		if (EQ(t, T_jikei_min)){ f_init = 1; }
		else if (EQ(t, T_jikei_max)){ f_init = -1;}
		else{ f_init = 0; }
		
		output_jikei_sub(ns, f_init, t, phi_v, r_order);
	}
}

void output_jikei_sub(long ns, long f_init, double t, double phi_v[], double r_order){
	long i;
	char text1[200];
	if (f_init == 1){ 
		my_fopen(&Fp_jikei, "out_jikei", 1); 
		
		fprintf(Fp_jikei, "t, ");
		for (i = 0; Id_phi_sel[i] != NER; i++) fprintf(Fp_jikei, "phi_v[%ld]%s, ", Id_phi_sel[i], Dlb_main);
		fprintf(Fp_jikei, "r_order%s\n", Dlb_main);
	}
	
	if (f_init == -1){ fclose(Fp_jikei); return; }
	
	if (ns % Nstep_jikei != 0) return;
	
	fprintf(Fp_jikei, "%.7le, ", t);
	for (i = 0; Id_phi_sel[i] != NER; i++) fprintf(Fp_jikei, "%.7le, ", phi_v[Id_phi_sel[i]]);
	fprintf(Fp_jikei, "%.7le\n", r_order);
}	

//Runge Kutta
void ruku(double t, double x[], double h){
	long i;
	double k1[NUM_X_ALL], k2[NUM_X_ALL], k3[NUM_X_ALL], k4[NUM_X_ALL];
	double wk1[NUM_X_ALL];
	
	//1
	c_x_dt(t, x, k1);
	
	//2
	for (i = 0; i < NUM_X_ALL; i++) wk1[i] = x[i] + h * k1[i]/2.e0;
	c_x_dt(t + h/2.e0, wk1, k2);
	
	//3
	for (i = 0; i < NUM_X_ALL; i++) wk1[i] = x[i] + h * k2[i]/2.e0;
	c_x_dt(t + h/2.e0, wk1, k3);
	
	//4
	for (i = 0; i < NUM_X_ALL; i++) wk1[i] = x[i] + h * k3[i];
	c_x_dt(t + h, wk1, k4);

	for (i = 0; i < NUM_X_ALL; i++){
		x[i] = x[i] + h * (k1[i]/6.e0 + k2[i]/3.e0 + k3[i]/3.e0 + k4[i]/6.e0);
	}
}

void my_fopen(FILE **fp, char fname[], long n_label){
	static long n_fop = 0;
	char text_wk1[200], text_folder[200];
	long i;	
	
	if (n_fop >= MAX_FILE_OPEN){ disp_err("err in my-fopen for opening too many files\n"); }
	if (n_label >= 2) disp_err("err my_fopen n-label\n");
	
	sprintf(text_wk1, "%s\\%s", OUT_FOLDER, fname);
	if (n_label == 1) strcat(text_wk1, Dlb_main);
	strcat(text_wk1, ".dat");
		
	if ( (*fp = fopen(text_wk1, "w") ) == NULL ){ 
		printf("err in fopcl for write mode open %s", text_wk1); 
		exit(1);
	}
	n_fop++;
}

void set_dlb_main(int argc, char *argv[]){
	if ( !strcmp(DLB_MAIN_BASE, "PRG") ){
		sprintf(Dlb_main, "(%s", argv[0]);
	}
	else{
		sprintf(Dlb_main, "(%s", DLB_MAIN_BASE);
	}
	
	strcat(Dlb_main, ")");
}

//random number (0~1) generation routine
float ransuu(void){
	if (RANSUU_METHOD == 2){
		return ran2_RCP(&I_ransuu);
	}
	else if (RANSUU_METHOD == 1){
		return ran1_RCP(&I_ransuu);
	}
	else if (RANSUU_METHOD == 0){
		return ran_fast_RCP(&I_ransuu);
	}
	else{
		printf("err in ran-suu\n"); 
		exit(1);
	}
}

//random number seed
void ransuu_init(long i_seed){
	if (i_seed >= 0){ printf("err i_seed >= 0\n"); exit(1); }
	I_ransuu = i_seed;
}

long n_ransuu(long n){
	long i_out;
	i_out = n * ransuu();
	if (i_out < 0 || i_out >= n) disp_err("err: n-ransuu");
	return i_out;
}

//return date and time
void disp_start_end_time(char md[], char c[]){
	FILE *fp;
	if ( (fp = fopen("time_now.dat", "a") ) == NULL){printf("err DPSET\n");exit(1); }
	
	c_time_now(Time_text);
	if (!strcmp(md, "s")){
		fprintf(fp, "start time of %s = %s\n", c, Time_text);
	}
	else{
		fprintf(fp, "end time of %s = %s\n", c, Time_text);
	}
	fclose(fp);
}

void disp_err(char c[]){
	printf("%s\n", c);
	exit(1);
}

