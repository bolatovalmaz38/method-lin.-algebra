#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include "rotation.h"

#define MAX_PRINT 10 

void multiply(double *A, double *b, double *c, int n);
double difference(double *a, int n);
void initialize(double *a, double *b, int n);
void initialize_m(double *a, int n);
void reverse(double *a, double *b, int n);                  
void *process_function(void *pa);
typedef struct Thread_data{
	pthread_t id;
	double *a, *x;
	int N;
	pthread_barrier_t *b;
}Thread_Data;
static pthread_barrier_t barrier;
static Thread_Data *threads;
static double *cos_phi;
static double *sin_phi;
int col;

static int lines_num;
static int columns_num;
static int coun;
static int N;
struct timespec diff(struct timespec start, struct timespec end);
struct timespec diff(struct timespec start, struct timespec end)
{
    struct timespec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}

int main(int argc, char *argv[])
{
  double l, *a, *x, *ax;  
  int i,j,n, p;
  FILE *f;
  struct timespec time_start, time_end;
  pthread_t *th;
        switch(argc) {

                case 2: {
                        f = fopen(argv[1], "r");
                        if (f == NULL){
                                printf("Error: 2 \n");
				fclose(f);
                                return -2;
                        }
                        if(fscanf(f, "%d", &n) == EOF) {
                                printf("Error: 3 \n");
				fclose(f);
                                return -3;
                        }
                        if(n <= 0) {
                                printf("Error: 4 \n");
				fclose(f);
                                return -4;
                        }
                        if(fscanf(f, "%d", &p) == EOF){
                                printf("Error: 5 \n");
				fclose(f);
                                return -5;
                        }
                        if(p <= 0 || p > 4) {
                                printf("Error: 6 \n");
				fclose(f);
                                return -6;
                        }
                        a = (double*)malloc((n * n)*sizeof(double));
                        for(i = 0; i < n; i++) 
                        {
                                for(j = 0; j < n; j++) {
                                        l = fscanf(f,"%lf", &a[i*n+j]);
                                        if(l == EOF || l == 0) {
                                                printf("Error: 7 \n");
        					free(a);
						fclose(f);
                                                return -7;
                                        }
                                }
                        }
                        fclose(f);
                        break;
                }

                case 3: {
                        sscanf(argv[1], "%d", &n);
                        sscanf(argv[2], "%d", &p);
                        a = (double*)malloc((n * n) * sizeof(double));
                        if(n <= 0){
                                printf("Error 1 \n");
        			free(a);
                                return -1;
                        }
                        if(p <= 0 || p > n){
                                printf("Error 3 \n");
        			free(a);
                                return -3;
                        }
                        break;
                }
                default:{
                	printf("Неверное количество аргументов\n");
                	return -1;
                }
	}
  x = (double *)malloc(n * n * sizeof(double));
  ax = (double *)malloc(n * n * sizeof(double));
  
  initialize(a, x, n); 

  if (!(threads = (Thread_Data*) malloc (p * sizeof(Thread_Data)))) {
  	free(a);
  	free(x);
  	free(ax);
  	return 1;
  }
  if (!(cos_phi = (double*) malloc((n-1)*sizeof(double)))) {
  	free(a);
  	free(x);
  	free(ax);
  	free(threads);
   	return 1;
  }
  if (!(sin_phi = (double*) malloc((n-1)*sizeof(double)))) {
  	free(a);
  	free(x);
  	free(ax);
  	free(threads);
    	free(cos_phi);
  	return 1;
  }
  if (!(th = (pthread_t*) malloc (p * sizeof(pthread_t)))){  
  	free(a);
  	free(x);
  	free(ax);
  	free(threads);
    	free(cos_phi);
    	free(sin_phi);
  	return 1;
  }
  N = n;
  coun = p;
  col = N / coun;
  clock_gettime(CLOCK_MONOTONIC, &time_start);
  pthread_barrier_init(&barrier, NULL, coun);
  for (i = 0; i < coun; i++){
    	threads[i].a = a;
    	threads[i].x = x;
	threads[i].N = i;
	threads[i].b = &barrier;
    	pthread_create(th + i, 0, process_function, threads + i);
   }
   for (i = 0; i < coun; i++)
     pthread_join(th[i], 0);
  pthread_barrier_destroy(&barrier);
  clock_gettime(CLOCK_MONOTONIC, &time_end);
  time_end = diff(time_start, time_end);
  printf("Time  = %lf ", time_end.tv_sec+(time_end.tv_nsec)/1000000000.0);
  printf("\n");

  initialize_m(a, n);
  multiply(a, x, ax, n);  
  printf("\n norma: %.10e\n", difference(ax, n));

  free(a);
  free(x);
  free(ax);
  free(threads);
  free(th);
  free(cos_phi);
  free(sin_phi);
  return 0;
}


void multiply(double *a, double *x, double *ax, int n)
{
  int i, j, k;

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++){
	for(k = 0; k < n; k++){
		ax[i*n+j] += a[i*n+k] * x[k*n + j];
	}
    }
  }

}

double difference(double *a, int n)
{
  int i, j;
  double max,m = 0;

  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
    	max += fabs(a[i*n+j]);
    }
    if (max > m)
	    m = max;
	max = 0;
  }
  return m-1;
}

void initialize(double *a, double *x, int n)
{
  int i, j;

  for (i = 0; i < n; i++)    
     for (j = 0; j < n; j++){
       //a[i*n + j] = 1.0/(i + j + 1);
       a[i*n + j] = fabs((double)(i - j));
       x[i*n+j] = (double)(i == j);
   }
}
void initialize_m(double *a, int n)
{
  int i, j;

  for (i = 0; i < n; i++)
     for (j = 0; j < n; j++){
       //a[i*n + j] = 1.0/(i + j + 1);
       a[i*n + j] = fabs((double)(i - j));
   }
}

void *process_function(void *pa)
{
  int k;
  double *A, *X;
  int i, j, l;
  double tmp;
  double *q, *pb;
  double *p, *p1, *p2, *cosp = cos_phi, *sinp = sin_phi;
  double cos_, sin_, x, y;
  int first_row, last_row;
  struct timespec time_start, time_end;
  double tmp1;
  Thread_Data *z = (Thread_Data*)pa;
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time_start);
  for (k = 1, A = z->a, X = z->x; k < N; k++, A += (N+1), X += N){
	if (z->N == 0){
		p = A + N;
		cosp = cos_phi;
		sinp = sin_phi;
		for (i = 1; i <= (N-k); i++, p += N, cosp++, sinp++) {
			x = *A;
			y = *p;
			tmp = sqrt(x*x + y*y);
			*cosp = cos_ = x / tmp;
			*sinp = sin_ = -y / tmp;
			*A = tmp;
			*p = 0;
		}
		lines_num = N - k;
		columns_num = (N - k) / coun;
		l = (N - k) % coun;
	}
	pthread_barrier_wait(z->b);
  	p = A + 1;
	cosp = cos_phi;
	sinp = sin_phi;
	if (columns_num > 0){
	for (i = 0, q = p + N; i < lines_num; i++, q += N, cosp++, sinp++) {
		cos_ = *cosp;
	   	sin_ = *sinp;
	    	for (j = 0, p1 = p + z->N*columns_num, p2 = q + z->N*columns_num; j < columns_num; j++, p1++, p2++) {
	      		x = *p1;
	      		y = *p2;
			*p1 = x*cos_ - y*sin_;
	      		*p2 = x*sin_ + y*cos_;
	    	}
	}
	}
	if (l > 0 && z->N == 0){
		p = A + 1 + columns_num*coun;
		cosp = cos_phi;
		sinp = sin_phi;
		for (i = 0, q = p + N; i < lines_num; i++, q += N, cosp++, sinp++) {
			cos_ = *cosp;
	    		sin_ = *sinp;
	    		for (j = 0, p1 = p, p2 = q; j < l; j++, p1++, p2++) {
	      			x = *p1;
	      			y = *p2;
	      			*p1 = x*cos_ - y*sin_;
	      			*p2 = x*sin_ + y*cos_;
	    		}
		}
  	}
	pthread_barrier_wait(z->b);
	p = X + col * z->N;
	cosp = cos_phi;
	sinp = sin_phi;
	if (z->N < coun-1){
	for(i = 0, pb = p + N; i < lines_num; pb += N, i++, cosp++, sinp++){
		cos_ = *cosp;
		sin_ = *sinp;
  		for (j = 0, p1 = p, p2 = pb; j < col; j++, p1++, p2++){
  			x = *p1;
  			y = *p2;
  			*p1 = x*cos_ - y*sin_;
  			*p2 = x*sin_ + y*cos_;
  		}
	}
	}
	else{
	l = col + N % coun;
	for(i = 0, pb = p + N; i < lines_num; pb += N, i++, cosp++, sinp++){
		cos_ = *cosp;
		sin_ = *sinp;
  		for (j = 0, p1 = p, p2 = pb; j < l; j++, p1++, p2++){
  			x = *p1;
  			y = *p2;
  			*p1 = x*cos_ - y*sin_;
  			*p2 = x*sin_ + y*cos_;
  		}
	}
	}
	pthread_barrier_wait(z->b);
}
	first_row = N * z->N;
	first_row /= coun;
	last_row = N * (z->N + 1);
	last_row /= coun;
	for (k = first_row; k < last_row; k++){
		for (i = N - 1; i >= 0; i--){
			tmp1 = z->x[i*N + k];
			for (j = i+1; j < N; j++)
				tmp1 -= z->a[i*N+j] * z->x[j*N+k];
			z->x[i*N+k] = tmp1/ z->a[i*N+i];
		}
	}
	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time_end);
        time_end = diff(time_start, time_end);
        printf("number %d Time  = %lf ",z->N, time_end.tv_sec+(time_end.tv_nsec)/1000000000.0);
        printf("\n");
	return 0;
}
