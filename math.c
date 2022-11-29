//random number (0~1) generation routine
float ran2_RCP(long *idum);
float ran1_RCP(long *idum);
float ran_fast_RCP(long *idum);


//constants used in ran2_RCP
#define IM1_RAN2_RCP	2147483563
#define IM2_RAN2_RCP	2147483399
#define AM_RAN2_RCP 	(1.0/IM1_RAN2_RCP)
#define IMM1_RAN2_RCP	(IM1_RAN2_RCP-1)
#define	IA1_RAN2_RCP	40014
#define IA2_RAN2_RCP	40692
#define IQ1_RAN2_RCP	53668
#define IQ2_RAN2_RCP	52774
#define IR1_RAN2_RCP	12211
#define IR2_RAN2_RCP	3791
#define NTAB_RAN2_RCP	32
#define NDIV_RAN2_RCP	(1+IMM1_RAN2_RCP/NTAB_RAN2_RCP)
#define EPS_RAN2_RCP	1.2e-7
#define RNMX_RAN2_RCP	(1.0-EPS_RAN2_RCP)

//constants used in ran1_RCP
#define	IA_RAN1_RCP		16807
#define	IM_RAN1_RCP		2147483647
#define	AM_RAN1_RCP 	(1.0/IM_RAN1_RCP)
#define	IQ_RAN1_RCP		127773
#define	IR_RAN1_RCP		2836
#define	NTAB_RAN1_RCP	32
#define	NDIV_RAN1_RCP	(1+(IM_RAN1_RCP-1)/NTAB_RAN1_RCP)
#define EPS_RAN1_RCP	1.2e-7
#define RNMX_RAN1_RCP	(1.0-EPS_RAN1_RCP)

//constants used in ran_fast_RCP
#define	IM_RAN_FAST_RCP 	714025
#define AM_RAN_FAST_RCP 	(1.0/IM_RAN_FAST_RCP)
#define IA_RAN_FAST_RCP 	1366
#define IC_RAN_FAST_RCP 	150889
#define EPS_RAN_FAST_RCP 	1.2e-7
#define RNMX_RAN_FAST_RCP 	(1.0-EPS_RAN_FAST_RCP)


//the number of random numbers is 2.3Å~10^18
float ran2_RCP(long *idum){
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB_RAN2_RCP];
	float temp;
	
	//init
	if (*idum <= 0){
		if ( -(*idum) < 1 ) *idum = 1;
		else *idum = -(*idum);
		idum2 = (*idum);
		for (j = NTAB_RAN2_RCP + 7; j >= 0; j--){
			k = (*idum) / IQ1_RAN2_RCP;
			*idum = IA1_RAN2_RCP * (*idum - k * IQ1_RAN2_RCP) - k * IR1_RAN2_RCP;
			if (*idum < 0) *idum += IM1_RAN2_RCP;
			if (j < NTAB_RAN2_RCP) iv[j] = *idum;
		}
		iy = iv[0];
	}
	
	k = (*idum) / IQ1_RAN2_RCP;
	*idum = IA1_RAN2_RCP * (*idum - k * IQ1_RAN2_RCP) - k * IR1_RAN2_RCP;
	if (*idum < 0) *idum += IM1_RAN2_RCP;
	
	k = idum2 / IQ2_RAN2_RCP;
	idum2 = IA2_RAN2_RCP * (idum2 - k * IQ2_RAN2_RCP) - k * IR2_RAN2_RCP;
	if (idum2 < 0) idum2 += IM2_RAN2_RCP;
	
	j = iy / NDIV_RAN2_RCP;
	iy = iv[j] - idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1_RAN2_RCP;
	
	if ( (temp = AM_RAN2_RCP * iy) > RNMX_RAN2_RCP ) return RNMX_RAN2_RCP;
	else return temp;
}

//the number of random numbers is 2.1Å~10^9
float ran1_RCP(long *idum){
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB_RAN1_RCP];
	float temp;
	
	//init
	if (*idum <= 0 || !iy){
		if ( -(*idum) < 1 ) *idum = 1;
		else *idum = -(*idum);
		for (j = NTAB_RAN1_RCP + 7; j >= 0; j--){
			k = (*idum) / IQ_RAN1_RCP;
			*idum = IA_RAN1_RCP * (*idum - k * IQ_RAN1_RCP) - IR_RAN1_RCP * k;
			if (*idum < 0) *idum += IM_RAN1_RCP;
			if (j < NTAB_RAN1_RCP) iv[j] = *idum;
		}
		iy = iv[0];
	}
	
	k = (*idum) / IQ_RAN1_RCP;
	*idum = IA_RAN1_RCP * (*idum - k * IQ_RAN1_RCP) - IR_RAN1_RCP * k;
	if (*idum < 0) *idum += IM_RAN1_RCP;
	j  = iy / NDIV_RAN1_RCP;
	iy = iv[j];
	iv[j] = *idum;
	if ((temp=AM_RAN1_RCP*iy) > RNMX_RAN1_RCP) return RNMX_RAN1_RCP;
	else return temp;
}

//the number of random numbers is 714025
float ran_fast_RCP(long *idum){
	float temp;
	if (*idum < 0){
		*idum = -(*idum);
	}
	*idum = ( (*idum) * IA_RAN_FAST_RCP + IC_RAN_FAST_RCP) % IM_RAN_FAST_RCP;
	temp= AM_RAN_FAST_RCP * (*idum);
	if (temp > RNMX_RAN_FAST_RCP) return RNMX_RAN_FAST_RCP;
	else return temp;
}



