#ifndef IO_H
#define IO_H

#include "pint.h"
#include "tipsy.h"

void print_queue_helper(struct particle *p, int lvl);
void print_queue(struct particle *p);
int loadascii(char *filename, struct env *env);
int loadtipsy(TCTX tctx, struct env *env);
int judge_and_load_file(char *filename, struct env *env);
void write_tipsy(char *filename, int seqno, struct env *env);
void write_ascii_fp(FILE *fp, struct env *env);
void write_ascii(char *filename, int seqno, struct env *env);
void write_density(FILE *fp, struct env *env);
void write_timesteps(FILE *fp, struct env *env);
void write_energy(FILE *fp, double time, struct env *env);
void write_mass(FILE *fp, double time, struct env *env);
void ic_binarystar(struct env *env);
void ic_earthsun(struct env *env);
void ic_randomsphere(struct env *env);
void ic_threebody(struct env *env);

void write_state_every(double timeNext, int maxTime);
void write_state(double timeNext);

#endif

