// João Vitor Oliveski Mesquita
// Março de 2022
// gcc [-O3] langevin.c -o exe -lm
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h> // RNG
#include <gsl/gsl_randist.h> // distribuição gaussiana

typedef struct{
	double x;
	double v;
}Particula;

typedef struct{
	double gamma;
	double T;
	double Kb;
	double k;
	double m;
	double tau;
	double Gamma;
}Parametros;

/*******************************
*			       *
*     Declarações de funções   *
*			       *
********************************/

void iniciaParticulas(Particula *particulas, int N){
	int i;

	for(i = 0; i < N; i++){
		particulas[i].x = 0;
		particulas[i].v = 0;
	}
}

// Força externa
double forca(double x, double k){
	return -k*x;
}

// recebe uma partícula e evolui ela
// em um passo no tempo
void passoSDE(Particula *particula, Parametros p, gsl_rng *r){
	double fx;	// força externa
	double v, x;
	v = particula->v;	// apelidos
	x = particula->x;

	fx = forca(x, p.k) / p.m;
	v = (1 - p.gamma*p.tau)*v + (fx + sqrt(p.Gamma/p.tau)*gsl_ran_gaussian(r,1))*p.tau;
	x = x + v*p.tau;
	particula->v = v;	// atualizo velocidade
	particula->x = x;	// e posição da partícula
}

// Media sobre todas particulas
void mediaPart(Particula *particulas, int N, Parametros p, double *Ec, double *Ep){
	double v2, x2;
	double acum_Ec, acum_Ep;
	int i;

	acum_Ec = 0;
	acum_Ep = 0;
	for(i = 0; i < N; i++){
		v2 = particulas[i].v * particulas[i].v;
		x2 = particulas[i].x * particulas[i].x;
		acum_Ec += 0.5 * p.m * v2;
		acum_Ep += 0.5 * p.k * x2;
	}

	*Ec = acum_Ec / N;
	*Ep = acum_Ep / N;
}

// autocorrelação da velocidade de uma partícula
void autocorrVel(double *vel_i, double dt, double t_max, double *corr_i){
	// t_max é o tempo máximo da medida da correlação
	int t, passo, passo_max = t_max/dt;
	int n;

	for(passo = 0; passo < passo_max; passo++){
		t = 0;
		n = 0;
		corr_i[passo] = 0;
		while(t + passo <= passo_max){
			corr_i[passo] = corr_i[passo] + vel_i[t]*vel_i[passo + t];
			n++;
			t++;
		}
		t = 0;
		corr_i[passo] /= n;
	}
}

// media normalizada das autocorrelações sobre todas partículas
void mediaCorr(double **corr, int N, double dt, double t_max, double *corr_media){
	int passo_max = t_max/dt;
	int i, j;
	double acum;

	// media sobre as partículas
	for(i = 0; i < passo_max; i++){
		acum = 0;
		for(j = 0; j < N; j++) acum = acum + corr[j][i];
		corr_media[i] = acum / N;
	}

	// normalização
	for(i = 0; i < passo_max; i++){
		corr_media[i] /= corr_media[0];
	}
}

int main(){

	FILE *file;
	const int N = 1000;
	double t, dt, t_max; // tempo e seus amigos
	int din_size;
	Particula *particulas;
	Parametros par;
	int i, passo = 0;
	double Ec, Ep;
	double *vel[N];
	double *corr[N];
	double *corr_media;
	
	Ec = 0;	// cond inicial
	Ep = 0;

	par.gamma = 2.0;
	par.T = 0.8;
	par.Kb = 1;
	par.k = 3;
	par.m = 1;
	par.tau = 0.01;
	par.Gamma = 2*par.gamma*par.Kb*par.T/par.m;

	t = 0;
	dt = par.tau;
	t_max = 100;
	din_size = t_max / dt;	// número de passos na dinâmica

	particulas = (Particula *) malloc(N*sizeof(Particula));
	corr_media = (double *) malloc(din_size*sizeof(double));
	for(i = 0; i < N; i++){
		vel[i] = (double *) malloc(din_size*sizeof(double));
		vel[i][0] = 0;
		corr[i] = (double *) malloc(din_size*sizeof(double));
	}
	
	gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, 11111);	// semente qualquer

	// Todas partículas com x(0) = 0, v(0) = 0;
	iniciaParticulas(particulas, N);
	file = fopen("medidas.dat", "w");
	fprintf(file, "%lf %lf %lf\n", t, Ec, Ep); // escreve cond inicial

	// faz as realizações da SDE pra cada partícula
	while(t < t_max){
		passo += 1;
		for(i = 0; i < N; i++){
			passoSDE(&particulas[i], par, r);
			vel[i][passo] = particulas[i].v;
		}
		mediaPart(particulas, N, par, &Ec, &Ep);
		t += dt;
		fprintf(file, "%lf %lf %lf\n", t, Ec, Ep);
	}
	fclose(file);

	// mede a autocorrelação pra cada partícula
	for(i = 0; i < N; i++){
		autocorrVel(vel[i], dt, 10, corr[i]);
	}

	// autocorrelação média
	// tempo máximo 10 arbitrário
	mediaCorr(corr, N, dt, 10, corr_media);

	// escreve a autocorrelação média
	file = fopen("autocorrelação.dat", "w");
	for(i = 0; i < (10/dt); i++){
		fprintf(file, "%lf %lf\n", i*dt, corr_media[i]);
	}
	fclose(file);

	// cleaning up
	free(particulas);
	free(corr_media);
	for(i = 0; i < N; i++){
		free(vel[i]);
		free(corr[i]);
	}

	return 0;
}
