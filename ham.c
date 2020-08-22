//needs cleaning

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <gsl_matrix.h>
#include <errno.h>
#include <gsl_errno.h>
#include <omp.h>
#include <math.h>

int main(int argc, char* argv[]){

	omp_set_num_threads(omp_get_max_threads());
	omp_set_dynamic(0);

	int rows = 9999, cols = 26, num_windows=200, iter=0;
	long double T = 310, KbT=0.0019872041*T, beta = 1/KbT, change = 1.0, diff;
	long double * Fx; Fx = (int *)malloc(sizeof(long double)*num_windows);
	long double * prevFx; prevFx = (int *)malloc(sizeof(long double)*num_windows);
	long double * Fx_old; Fx_old = (int *)malloc(sizeof(long double)*num_windows);
	long double * UB; UB = (int *)malloc(sizeof(long double)*num_windows*num_windows*rows);
	long double * UBk; UBk = (int *)malloc(sizeof(long double)*num_windows);

	if(Fx==NULL || Fx_old==NULL || UB==NULL){
		printf("malloc failed\n");
		exit(-1);
	}

	gsl_matrix * m[num_windows];			//create array of pointers that point to gsl_matrices
	gsl_matrix * centers; FILE * file; gsl_matrix * kmat;
	FILE * out;

	for(int i=0;i<num_windows;i++){
		m[i] = gsl_matrix_alloc(rows,cols);
		if(m[i]==NULL){
			GSL_ERROR("memory", GSL_ENOMEM);
		}
	}
	#pragma omp parallel for
	for(int i=0;i<num_windows;i++){
		FILE *f; char * filename;
		f = NULL;
		filename = malloc(sizeof(char)*200);		//dynamic allocation, otherwise buffer overflow
		sprintf(filename, "data/data_s100_string_%d.dat", i+1);
		f = fopen(filename, "r");
		gsl_matrix_fscanf(f,m[i]);
		printf("read window %d\n", i+1);
		free(filename);
	}
	//printf("%f\n", gsl_matrix_get(m[7], 0, 0));
	
	centers = gsl_matrix_alloc(num_windows, cols);
	file = fopen("../optstring/newconstr_3.dat", "r");
	int z = gsl_matrix_fscanf(file, centers);
	fclose(file);

	kmat = gsl_matrix_alloc(num_windows, cols);
	gsl_matrix_set_all(kmat, 20);
	int q = 0;

	//printf("%f\n", pow((gsl_matrix_get(centers, 155, 10) - gsl_matrix_get(m[2], 1000, 10)), 2));

	//#pragma omp parallel for
	for(int i=0;i<num_windows;i++){
		printf("eval. win. %d\n", i);
		for(int j=0;j<(m[i]->size1);j++){
			for(int k=0;k<num_windows;k++){
				long double Ubias = 0.0;
				for(int s=0;s<cols;s++){
					Ubias += 0.5*gsl_matrix_get(kmat, k, s)*pow((gsl_matrix_get(centers, k, s) - gsl_matrix_get(m[i], j, s)), 2);
				}
				UB[q] = exp(-Ubias*beta);
				q++;
			}
		}
	}
	printf("evaluated unbiasing factors\n");

	#pragma omp parallel for
	for(int i=0;i<num_windows;i++){
		Fx_old[i] = 0.0;
	}

	while(change>0.01){

		if(iter>0){
			#pragma omp parallel for
			for(int i=0;i<num_windows;i++){
				prevFx[i] = Fx_old[i];
			}
		}

		//#pragma omp parallel for
		for(int i=0;i<num_windows;i++){
			Fx_old[i] = (m[i]->size1)*exp(Fx_old[i]*beta);
			//printf("%d %f\n", i, Fx_old[i]);
		}
		//#pragma omp parallel for
		for(int i=0;i<num_windows;i++){
			Fx[i] = 0.0;
		}
		int q = 0;
	
		for(int i=0;i<num_windows;i++){
			for(int j=0;j<(m[i]->size1);j++){
				long double denom = 0;
				for(int k=0;k<num_windows;k++){
					UBk[k] = UB[q];
					denom += Fx_old[k]*UB[q];
					q++;
				}
				denom = 1/denom;	//check if this is inf
				for(int k=0;k<num_windows;k++){
					Fx[k] += UBk[k]*denom;
				}
			}
		}
		//#pragma omp parallel for
		for(int i=0;i<num_windows;i++){
			//printf("%d %f\n", i, Fx[i]);
			Fx[i] = -KbT*log(Fx[i]);
			printf("%d %f\n", i, Fx[i]);
		}
		//#pragma omp parallel for
		for(int i=1;i<num_windows;i++){
			Fx[i] = Fx[i] - Fx[0];
		}
		//#pragma omp parallel for
		for(int i=0;i<num_windows;i++){
			Fx_old[i] = Fx[i];
		}

		if(iter>0){
			long double diff1 = 0;
			for(int i=1;i<num_windows;i++){
				diff = abs(Fx[i]-prevFx[i]);
				if(diff > diff1){
					diff1 = diff;
				}
			}
			change = diff1;
		}
		printf("iteration %d, change = %Lf\n", iter, change);
		iter++;
	}
	printf("%d\n", q);

	out = fopen("Fx.dat", "w");
	for(int i=0;i<200;i++){
		fprintf(out, "%f\n", Fx_old[i]);
	}
	fclose(out);
	free(Fx); free(Fx_old); free(prevFx); free(UB); free(UBk);
}
