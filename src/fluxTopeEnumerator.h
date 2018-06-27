///////////////////////////////////////////////////////////////////////////////
// Author: Matthias P. Gerstl
// Email: matthias.gerstl@acib.at
// Company: Austrian Centre of Industrial Biotechnology (ACIB)
// Web: http://www.acib.at Copyright
// (C) 2016 Published unter GNU Public License V3
///////////////////////////////////////////////////////////////////////////////
//Basic Permissions.
//
// All rights granted under this License are granted for the term of copyright
// on the Program, and are irrevocable provided the stated conditions are met.
// This License explicitly affirms your unlimited permission to run the
// unmodified Program. The output from running a covered work is covered by
// this License only if the output, given its content, constitutes a covered
// work. This License acknowledges your rights of fair use or other equivalent,
// as provided by copyright law.
//
// You may make, run and propagate covered works that you do not convey,
// without conditions so long as your license otherwise remains in force. You
// may convey covered works to others for the sole purpose of having them make
// modifications exclusively for you, or provide you with facilities for
// running those works, provided that you comply with the terms of this License
// in conveying all material for which you do not control copyright. Those thus
// making or running the covered works for you must do so exclusively on your
// behalf, under your direction and control, on terms that prohibit them from
// making any copies of your copyrighted material outside their relationship
// with you.
//
// Disclaimer of Warranty.
//
// THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE
// LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR
// OTHER PARTIES PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND,
// EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
// ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.
// SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY
// SERVICING, REPAIR OR CORRECTION.
//
// Limitation of Liability.
//
// IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL
// ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS THE
// PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY
// GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE
// OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA
// OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD
// PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
// EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGES.
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <ilcplex/cplex.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>

#define BITSIZE CHAR_BIT
#define BITMASK(b) (1 << ((b) % BITSIZE))
#define BITSLOT(b) ((b) / BITSIZE)
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a, b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a, b) ((a)[BITSLOT(b)] & BITMASK(b))

#define MAX_ARGS        19
#define ARG_FT_OUT       0
#define ARG_OUTFORMAT    1 
#define ARG_KERNEL       2
#define ARG_RM_LIST      3
#define ARG_SFILE        4
#define ARG_RFILE        5
#define ARG_RVFILE       6
#define ARG_MFILE        7
#define ARG_ZERO         8
#define ARG_THREADS      9
#define ARG_LP_OUT      10
#define ARG_LOG         11
#define ARG_MINFLUX     12
#define ARG_MAXFLUX     13
#define ARG_NET         14
#define ARG_CONFIGFILE  15
#define ARG_FT_IN       16
#define ARG_FT_START    17
#define ARG_FT_END      18

#define ERROR_ARGS     1
#define ERROR_TREE     2
#define ERROR_FILE     3
#define ERROR_RAM      4
#define ERROR_CPLEX    5

#define OUT_BIN 1
#define OUT_TXT 2

unsigned long long tree_exist;
unsigned long long tree_size;
unsigned long long act_tree_exist;
unsigned long long act_tree_size;
pthread_mutex_t tree_mutex;

struct thread_args
{
    char *bv;
    int bv_size;
    char *rev;
    int k_rev_len;
    int cols;
    char ***adj_bv;
    int *adj_count;
    int thread_id;
    int *finished;
    int net;
    int **map;
    int *maps;
    int all_rx_len;
    double *lb;
    double *ub;
    double min_flux;
    CPXENVptr env;
    CPXLPptr lp;
};

struct metabolite
{
    int use;
    int ignore;
    double dfg;
    double min;
    double max;
    char* name;
};

struct node
{
    char* value;
    struct node *left;
    struct node *right;
};

struct node *root;

// tree
void tree_destroy(struct node *leaf);
int tree_insert_not_safe(char* val, struct node **leaf, int size);
int tree_insert(char* val, struct node **leaf, int size);

// Bit
unsigned long getBitsize(unsigned long count);
void copyBitvector(char *from, char *to, int len);

// input output
void open_file(FILE **file, char *name, char *handle);
char* getTime();
void printHeader(FILE *f, char **names, int *rmlist, int len, char* rev, int
        rev_len, int format);
void printVector(FILE *f, char *v, int len, int **map, int *maps,
        int all_rx_len, int outformat, char* rev, int rev_len);
void readReversibleFile(char *filename, char *reversible_reactions,
        unsigned int rx_count, int *rev_rx_count);
void readRemoveList(char *filename, int *list, unsigned int rx_count);
void readInputList(char *filename, char*** rx_names, unsigned int rx_count);
void readMatrix(int ncol, int* nrow, double*** k_mat, char *filename);
void readFluxtopes(int ncol, unsigned long long* nrow, char*** ft, char
        *filename, int step, int** rx_map, int map_len, int all_rx_ct, char*
        rev, int rev_len, int log);
int addToFluxtopes(char*** ft, unsigned long* ft_ct, int step, int*
        act_step, int* step_ct, unsigned long long* act_size, int* end, int
        log);

// for reading fluxtopes 
// searches only first appearance in the reaction map
int getBitVectorPosition(int** map, int map_len, int orig_rx_ix);

// memory
void reallocDoubleMatrix (double*** matrix, unsigned long index);
void refreshDoubleMatrix (double*** matrix, unsigned long rows);
void allocateDoubleRow (double** vector, unsigned long reactions);

// arguments
void handle_arguments(int argc, char *argv[], char **ftout, int *outformat,
        char **kernel, char **rmlist, char **sfile, char **rfile, char
        **rvfile, char **mfile, double *eprhs, int *threads, char **lp_out, int
        *net, int *log, double *min_flux, double *max_flux, char **configfile,
        char **ft_in, int *ft_start, int *ft_end);

// LP
void init_cplex(CPXENVptr *env, double eprhs);
void createNetworkLpProblem (CPXENVptr env, CPXLPptr *lp, double **stoich,
        int s_rows, int s_cols, struct metabolite *metabolites,
        char** rx_names, double* lower_bounds,
        double* upper_bounds, double min_flux, char* name);
void createKernelLpProblem (CPXENVptr env, CPXLPptr *lp, double **kernel,
        int k_rows, int k_cols, char* name);
int checkFeasibility(CPXENVptr env, CPXLPptr lp);
int setNetworkBounds(char *bv, int bv_size, int ix, int **map, int *maps,
        int all_rx_len, CPXENVptr env, CPXLPptr lp, double *lb, double *ub,
        double min_flux);
void setKernelBounds(char *bv, int bv_size, int ix, int cols, CPXENVptr env,
        CPXLPptr lp);
void defineBounds(double **lower_bounds, double **upper_bounds, int s_cols,
        char** reaction_names, char* reversibles, double min_flux, double
        max_flux, char* configfile);

// general
char* getArg (int argc, char *argv[], char *opt);
void readArgs (int argc, char *argv[], int optc, char *optv[], char *optr[]);
void usage (char* description, char* usage, int optc, char* optv[],
        char *optd[]);
int get_col_count(char *filename);
int getRxCount(FILE* file);
void createReactionMapping(int *rmlist, char *reversible, int ***rx_map,
        char **kern_rev, int *map_size, int **map_lengths, int rx_count,
        int bitsize);
void quitError (char *message, int rv);
void quitErrorArg (char *message, char *arg, int rv);

// threads
int allThreadsFree(int threads, int *thread_lock);
int getFreeThread(int threads, int *thread_lock);
int getFinishedThread(int threads, int *thread_finished);

// flux topes
void* findAdj(void *pointer_thread_args);
void findFluxTopes(unsigned long long *foundCounter, char** bv, unsigned long
        long bv_len, int bv_size, char* rev, int k_rev_len, int cols, int
        threads, CPXENVptr *envs, CPXLPptr *lps, int **map, int *maps, int
        all_rx_len, double min_flux, int net, FILE *f, int log, double
        *lower_bounds, double *upper_bounds, int step, int endstep, int
        outformat);

void defineMetabolites(struct metabolite *metabolites, char **met_names, int
        met_count);

// combinatorics
double binom(int n, int k);
int min(int a, int b);
