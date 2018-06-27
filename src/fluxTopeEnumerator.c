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

#include "fluxTopeEnumerator.h"

int main (int argc, char *argv[])
{

    // handle arguments
    //==================================================
    char *a_kernel, *a_rmlist, *a_sfile, *a_rfile, *a_rvfile, *a_mfile,
         *a_ftout, *lp_out, *a_configfile, *a_ft_in;
    double eprhs, min_flux, max_flux;
    int threads, a_net, a_log, a_ft_start, a_ft_end, a_outformat;
    handle_arguments(argc, argv, &a_ftout, &a_outformat, &a_kernel, &a_rmlist,
            &a_sfile, &a_rfile, &a_rvfile, &a_mfile, &eprhs, &threads, &lp_out,
            &a_net, &a_log, &min_flux, &max_flux, &a_configfile, &a_ft_in,
            &a_ft_start, &a_ft_end);

    if (a_log > 0)
    {
        fprintf(stdout, "Start: %s\n", getTime());
    }

    // read column count of matrices
    //==================================================
    int s_cols = get_col_count(a_sfile);
    int k_cols = get_col_count(a_kernel);

    // define reversible reactions
    //==================================================
    unsigned long all_rx_bitsize = getBitsize(s_cols);
    char *reversible_reactions   = calloc(1, all_rx_bitsize);
    int rev_rx_count;
    readReversibleFile(a_rvfile, reversible_reactions, s_cols, &rev_rx_count);

    // read reaction names
    //==================================================
    char **reaction_names = NULL;
    readInputList(a_rfile, &reaction_names, s_cols);

    // set default lower and upper bounds
    //==================================================
    double *lower_bounds, *upper_bounds;
    defineBounds(&lower_bounds, &upper_bounds, s_cols, reaction_names,
            reversible_reactions, min_flux, max_flux, a_configfile);

    // read kernel
    //==================================================
    if (a_log > 0)
    {
        fprintf(stdout, "%s: read Kernel\n", getTime());
    }
    double **kernel;
    char *kern_rev;
    int k_rows, map_size, *map_lengths, **rx_map;
    int *rmlist = calloc(sizeof(int), s_cols);
    readMatrix(k_cols, &k_rows, &kernel, a_kernel);
    readRemoveList(a_rmlist, rmlist, s_cols);
    createReactionMapping(rmlist, reversible_reactions, &rx_map,
            &kern_rev, &map_size, &map_lengths, s_cols, all_rx_bitsize);
    
    // print header
    //==================================================
    FILE *f_ftout;
    open_file(&f_ftout, a_ftout, "w");
    printHeader(f_ftout, reaction_names, rmlist, s_cols, kern_rev, k_rows,
            a_outformat);

    // read stoichiometric matrix
    //==================================================
    double **stoich;
    int s_rows;
    readMatrix(s_cols, &s_rows, &stoich, a_sfile);

    // get metabolites 
    //==================================================
    char **met_names = NULL;
    readInputList(a_mfile, &met_names, s_rows);
    struct metabolite *metabolites = calloc(s_rows, sizeof(struct metabolite));
    defineMetabolites(metabolites, met_names, s_rows);

    // create CPLEX objects
    //==================================================
    CPXENVptr *envs = calloc(threads, sizeof(CPXENVptr));
    CPXLPptr *lps = calloc(threads, sizeof(CPXLPptr));
    if (a_log > 0)
    {
        fprintf(stdout, "%s: init CPLEX\n", getTime());
    }
    int i;
    for (i = 0; i < threads; i++)
    {
        init_cplex(&envs[i], eprhs);
        if (i < 1 && a_log > 0)
        {
            fprintf(stdout, "%s: create LP \n", getTime());
        }
        if (a_net)
        {
            createNetworkLpProblem(envs[i], &lps[i], stoich, s_rows, s_cols,
                    metabolites, reaction_names, lower_bounds, upper_bounds,
                    min_flux, "lp");
        }
        else
        {
            createKernelLpProblem(envs[i], &lps[i], kernel, k_rows, k_cols,
                    "lp");
        }
        // print LP file if filename is given
        if (i < 1)
        {
            if (strcmp(lp_out, "not_available") != 0)
            {
                if (a_log > 0)
                {
                    fprintf(stdout, "%s: write LP to %s\n", getTime(), lp_out);
                }
                CPXwriteprob (envs[0], lps[0], lp_out, "LP");
            }
            if (a_log > 0)
            {
                fprintf(stdout, "%s: check feasibility \n", getTime());
            }
            fflush(stdout);
            int feas = checkFeasibility(envs[i], lps[i]);
            if (feas != 0)
            {
                fprintf(stderr, "LP is not feasible\n");
                fprintf(stderr, "CPLEX returned status: %d\n", feas);
                exit(ERROR_CPLEX);
            }
        }
    }

    // free memory
    //==================================================
    for (i=0; i<k_rows; i++)
    {
        free(kernel[i]);
    }
    free(kernel);
    if (a_net)
    {
        for (i=0; i<s_rows; i++)
        {
            free(stoich[i]);
        }
    }
    free(stoich);

    // create full bitvector
    //==================================================
    unsigned long k_bitsize = getBitsize(k_rows);

    // create bitvector of reversible reactions
    //==================================================
    int k_rev_len = 0;
    for (i = 0; i < k_rows; i++)
    {
        if (BITTEST(kern_rev, i))
        {
            k_rev_len++;
        }
    }
    unsigned long k_rev_bitsize = getBitsize(k_rev_len);
    char *rev_bv = calloc(1, k_rev_bitsize);
    for (i = 0; i < k_rev_len; i++)
    {
        BITSET(rev_bv, i);
    }

    // LOG information
    //==================================================
    if (a_log > 0)
    {
        fprintf(stdout, "\nReversible reactions in kernel: %d\n", k_rev_len);
    }
    unsigned long long poss_tree_size = 0;
    for (i = 1; i < k_rev_len; i++)
    {
        poss_tree_size += binom(k_rev_len, i);
    }
    if (a_log > 0)
    {
        fprintf(stdout, "Possible tree entries: %llu\n\n", poss_tree_size);
    }

    // initialize tree
    //==================================================
    tree_size = 0;
    tree_exist = 0;
    int inserted = tree_insert(rev_bv, &root, k_rev_len);
    if (!inserted)
    {
        fprintf(stderr, "Could not insert into tree\n");
        exit(ERROR_TREE);
    }

    // run algorithm
    //==================================================
    int step = 1;
    unsigned long long foundCounter = 0;
    unsigned long long rows = 1;
    char **bvs;
    if (a_ft_start > 0){
        readFluxtopes(k_rows, &rows, &bvs, a_ft_in, a_ft_start - 1, rx_map,
                map_size, s_cols, kern_rev, k_rev_len, a_log);
        step = a_ft_start;
        fprintf(stdout, "%s: imported %llu fluxtopes\n", getTime(), rows);
        fflush(stdout);
    } else {
        char *row_bv = calloc(1, k_bitsize);
        for (i = 0; i < k_rows; i++)
        {
            BITSET(row_bv, i);
        }
        bvs = calloc(1, sizeof(char*));
        bvs[0] = row_bv;
    }
    if (a_log > 0)
    {
        fprintf(stdout, "%s: start algorithm to find flux topes\n", getTime());
        fflush(stdout);
    }
    findFluxTopes(&foundCounter, bvs, rows, k_rows, kern_rev, k_rev_len,
            k_cols, threads, envs, lps, rx_map, map_lengths, s_cols, min_flux,
            a_net, f_ftout, a_log, lower_bounds, upper_bounds, step, a_ft_end,
            a_outformat);
    if (a_log > 0)
    {
        fprintf(stdout, "Found flux topes: %llu\n", foundCounter);
    }

    fclose(f_ftout);

    // set memory free
    //==================================================
    free(lower_bounds);
    free(upper_bounds);
    tree_destroy(root);
    for (i = 0; i < s_cols; i++)
    {
        free(reaction_names[i]);
    }
    free(reaction_names);
    for (i = 0; i < s_rows; i++)
    {
        free(met_names[i]);
    }
    free(met_names);
    free(reversible_reactions);
    free(rmlist);
    for (i = 0; i < map_size; i++)
    {
        free(rx_map[i]);
    }
    free(rx_map);
    free(metabolites);
    free(map_lengths);
    free(kern_rev);
    for (i = 0; i < threads; i++)
    {
        CPXfreeprob(envs[i], &(lps[i]));
        CPXcloseCPLEX(&(envs[i]));
    }
    free(lps);
    free(envs);
    if (a_log > 0)
    {
        fprintf(stdout, "Tree size: %llu\n", tree_size);
        fprintf(stdout, "Found in tree: %llu\n", tree_exist);
        fprintf(stdout, "Finished: %s\n", getTime());
    }
    return EXIT_SUCCESS;
}

void tree_destroy(struct node *leaf)
{
    if( leaf != 0 )
    {
        tree_destroy(leaf->left);
        tree_destroy(leaf->right);
        free( leaf->value );
        free( leaf );
    }
}

int tree_insert_not_safe(char* val, struct node **leaf, int size)
{
    int return_val = 0;
    if( *leaf == 0 )
    {
        *leaf = (struct node*) malloc( sizeof( struct node ) );
        if (leaf == NULL)
        {
            fprintf(stderr, "Not enough memory for adding to the tree\n");
            exit(ERROR_RAM);
        }
        (*leaf)->value = val;
        (*leaf)->left = 0;
        (*leaf)->right = 0;
        tree_size++;
        act_tree_size++;
        return_val = 1;
    }
    else
    {
        int cmp = 0;
        int i;
        char* l_val = (*leaf)->value;
        for (i = 0; i < size; i++)
        {
            if (BITTEST(l_val, i) && !BITTEST(val, i))
            {
                cmp = -1;
                i = size;
            }
            else if (!BITTEST(l_val, i) && BITTEST(val, i))
            {
                cmp = 1;
                i = size;
            }
        }
        if (cmp < 0)
        {
            return_val = tree_insert_not_safe( val, &(*leaf)->left, size );
        }
        else if (cmp > 0)
        {
            return_val = tree_insert_not_safe( val, &(*leaf)->right, size );
        }
        else
        {
            tree_exist++;
            act_tree_exist++;
            return_val = 0;
        }
    }
    return return_val;
}

int tree_insert(char* val, struct node **leaf, int size)
{
    pthread_mutex_lock(&tree_mutex);
    int return_val = 0;
    if( *leaf == 0 )
    {
        *leaf = (struct node*) malloc( sizeof( struct node ) );
        if (leaf == NULL)
        {
            fprintf(stderr, "Not enough memory for adding to the tree\n");
            exit(ERROR_RAM);
        }
        (*leaf)->value = val;
        (*leaf)->left = 0;
        (*leaf)->right = 0;
        tree_size++;
        act_tree_size++;
        return_val = 1;
    }
    else
    {
        int cmp = 0;
        int i;
        char* l_val = (*leaf)->value;
        for (i = 0; i < size; i++)
        {
            if (BITTEST(l_val, i) && !BITTEST(val, i))
            {
                cmp = -1;
                i = size;
            }
            else if (!BITTEST(l_val, i) && BITTEST(val, i))
            {
                cmp = 1;
                i = size;
            }
        }
        if (cmp < 0)
        {
            return_val = tree_insert_not_safe( val, &(*leaf)->left, size );
        }
        else if (cmp > 0)
        {
            return_val = tree_insert_not_safe( val, &(*leaf)->right, size );
        }
        else
        {
            tree_exist++;
            act_tree_exist++;
            return_val = 0;
        }
    }
    pthread_mutex_unlock(&tree_mutex);
    return return_val;
}

/**
 * calculates needed number of chars for bit support of ltcs
 */
unsigned long getBitsize(unsigned long count)
{
    unsigned long bitsize = count / BITSIZE;
    int modulo  = count % BITSIZE;
    if (modulo > 0)
    {
        bitsize++;
    }
    return bitsize;
}

void copyBitvector(char *from, char *to, int len)
{
    int i;
    for (i = 0; i < len; i++)
    {
        if (BITTEST(from, i))
        {
            BITSET(to,i);
        }
        else
        {
            BITCLEAR(to,i);
        }
    }
}

void open_file(FILE **file, char *name, char *handle)
{
   *file = fopen(name, handle);
   if (*file == NULL)
   {
       quitErrorArg("Error in opening file: ", name, ERROR_FILE);
   }
}

void printHeader(FILE *f, char **names, int* rmlist, int len, char* rev, int
        rev_len, int format)
{
    int i;
    if (format == OUT_BIN)
    {
        fprintf(f, "compressedfluxtopes:%d:%d\n", len, rev_len);
        for (i = 0; i < len; i++)
        {
            if (i > 0)
            {
                fprintf(f, ",");
            }
            fprintf(f, "%d", rmlist[i]);
        }
        fprintf(f, "\n");
        for (i = 0; i < rev_len; i++) 
        {
            if (i > 0)
            {
                fprintf(f, ",");
            }
            if (BITTEST(rev,i))
            {
                fprintf(f, "1");
            }
            else
            {
                fprintf(f, "0");
            }
        }
        fprintf(f, "\n");
    }
    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            fprintf(f, ",");
        }
        fprintf(f, "%s", names[i]);
    }
    fprintf(f, "\n");
}

void printVector(FILE *f, char *v, int len, int **map, int *maps, int
        all_rx_len, int format, char* rev, int rev_len)
{
    if (format == OUT_BIN)
    {
        int bv_size = getBitsize(rev_len);
        char* cv = calloc(1, bv_size);
        int i;
        int t_ix = -1;
        for (i = 0; i < len; i++) 
        {
            if (BITTEST(rev, i))
            {
                t_ix++;
                if (BITTEST(v, i))
                {
                    BITSET(cv, t_ix);
                }
                else
                {
                    BITCLEAR(cv, t_ix);
                }
            }
        }
        fwrite(cv, 1, bv_size, f);
        free(cv);
    }
    else
    {
        char *exp = calloc(1, getBitsize(all_rx_len));
        if (exp == NULL)
        {
            fprintf(stderr, "Not enough RAM left\n");
            exit(ERROR_RAM);
        }
        int i;
        for (i = 0; i < len; i++)
        {
            if (BITTEST(v, i))
            {
                int j;
                for (j = 0; j < maps[i]; j++)
                {
                    BITSET(exp, map[i][j]);
                }
            }
        }
        for (i = 0; i < all_rx_len; i++)
        {
            if (i > 0)
            {
                fprintf(f, ",");
            }
            if (BITTEST(exp, i))
            {
                fprintf(f, "1");
            }
            else
            {
                fprintf(f, "-1");
            }
        }
        fprintf(f, "\n");
        free(exp);
    }
}


/*
 * return actual time as string
 */
char* getTime()
{
    static char time_string[20];
    time_t now = time (0);
    strftime (time_string, 100, "%Y-%m-%d %H:%M:%S", localtime (&now));
    return time_string;
}

/*
 * increase memory size of stoichiometric matrix
 */
void reallocDoubleMatrix (double*** matrix, unsigned long index)
{
    int alloc_size = 100;
    unsigned long new_size = 0;
    if (index % alloc_size == 0)
    {
        new_size = index + alloc_size;
        *matrix = (double**) realloc(*matrix, new_size * sizeof(double*));
        if (NULL == matrix)
        {
            quitError("Not enough free memory in reallocDoubleMatrix\n",
                    ERROR_RAM);
        }
    }
}

/*
 * remove not used allocations of a double matrix
 */
void refreshDoubleMatrix (double*** matrix, unsigned long rows)
{
    *matrix = (double**) realloc(*matrix, rows * sizeof(double*));
    if (NULL == matrix)
    {
        quitError("Not enough free memory in refreshDoubleMatrix\n",
                ERROR_RAM);
    }
}


/*
 * allocate memory for an stoichiometric row
 */
void allocateDoubleRow (double** vector, unsigned long reactions)
{
    *vector = calloc(sizeof(double), reactions);
    if (NULL == vector)
    {
        quitError("Not enough free memory\n", ERROR_RAM);
    }
}

/*
 * read file containing reversibilty of reactions in form of:
 * 1 0 0 1 1
 * where 1 describes a reversible and 0 a non reversible reaction
 */
void readReversibleFile(char *filename, char *reversible_reactions, unsigned
        int rx_count, int *rev_rx_count)
{
    char* line = NULL;
    size_t len = 0;
    FILE *file;
    open_file(&file, filename, "r");
    int counter = 0;
    while ( getline(&line, &len, file) != -1)
    {
        char *ptr;
        ptr = strtok(line, "\n\t ");
        int i = 0;
        while(ptr != NULL)
        {
            int x = atoi(ptr);
            if (x == 1)
            {
                BITSET(reversible_reactions, i);
                counter++;
            }
            else if (x == 0)
            {
                BITCLEAR(reversible_reactions, i);
            }
            else
            {
                quitError("not a correct reversibility file\n", ERROR_FILE);
            }
            i++;
            ptr = strtok(NULL, "\n\t ");
        }
        if (i != rx_count)
        {
            quitError("not a correct reversibility file\n", ERROR_FILE);
        }
    }
    fclose(file);
    free(line);
    *rev_rx_count = counter;
}

/*
 * read removed reaction list file
 */
void readRemoveList(char *filename, int *list, unsigned int rx_count)
{
    char* line = NULL;
    size_t len = 0;
    FILE *file;
    open_file(&file, filename, "r");
    int i = 0;
    while ( getline(&line, &len, file) != -1)
    {
        char *ptr;
        ptr = strtok(line, "\n\t ");
        if(ptr != NULL)
        {
            list[i] = atoi(ptr);
            i++;
        }
        ptr = strtok(NULL, "\n\t ");
    }
    if (i != rx_count)
    {
        fprintf(stderr, "rx_count = %d, rmlist count = %d\n", rx_count, i);
        quitError("not a correct rmlist file\n", ERROR_FILE);
    }
    fclose(file);
    free(line);
}

void readInputList(char *filename, char*** rx_names, unsigned int rx_count)
{
    char* line = NULL;
    size_t len = 0;
    FILE *file;
    open_file(&file, filename, "r");
    char** m_rx_names = calloc(rx_count, sizeof(char*));
    if ( getline(&line, &len, file) != -1)
    {
        char *ptr;
        ptr = strtok(line, "\"\n\t ");
        int i = 0;
        while(ptr != NULL)
        {
            if (i >= rx_count)
            {
                quitError("not correct reaction file: more reactions than in "
                        "EFM file\n", ERROR_FILE);
            }
            size_t len = strlen(ptr);
            m_rx_names[i] = calloc(len + 1, sizeof(char));
            strcpy(m_rx_names[i], ptr);
            ptr = strtok(NULL, "\"\n\t ");
            i++;
        }
        if (i != rx_count)
        {
            fprintf(stderr, "Number of reactions in sfile: %d\n", rx_count);
            fprintf(stderr, "Number of reactions in reaction file: %d\n", i);
            quitError("not a correct reaction file\n", ERROR_FILE);
        }
    }
    fclose(file);
    free(line);
    *rx_names = m_rx_names;
}

/*
 * read matrix
 */
void readMatrix(int ncol, int* nrow, double*** k_mat, char *filename)
{
    double** m_kernel = NULL;
    char* line = NULL;
    size_t len = 0;
    unsigned long mat_ix = 0;
    FILE *file;
    open_file(&file, filename, "r");
    while ( getline(&line, &len, file) != -1)
    {
        reallocDoubleMatrix(&m_kernel, mat_ix);
        allocateDoubleRow(&m_kernel[mat_ix], ncol);
        char *ptr;
        ptr = strtok(line, "\n\t ");
        int i = 0;
        while(ptr != NULL)
        {
            m_kernel[mat_ix][i] = atof(ptr);
            i++;
            ptr = strtok(NULL, "\n\t ");
        }
        mat_ix++;
    }
    refreshDoubleMatrix(&m_kernel, mat_ix);
    free(line);
    line = NULL;
    fclose(file);
    *nrow = mat_ix;
    *k_mat = m_kernel;
}

int getBitVectorPosition(int** map, int map_len, int orig_rx_ix)
{
    int i;
    for (i = 0; i < map_len; i++) {
        if (map[i][0] == orig_rx_ix)
        {
            return i;
        }
    }
    return -1;
}

int addToFluxtopes(char*** ft, unsigned long* ft_ct, int step, int*
        act_step, int* step_ct, unsigned long long* act_size, int* end, int
        log) 
{
    int result = 0;
    if (*step_ct == step)
    {
        *ft_ct = *ft_ct + 1;
        *ft = (char**) realloc(*ft, *ft_ct * sizeof(char*));
        if (ft == NULL)
        {
            fprintf(stderr, "Not enough RAM left\n");
            exit(ERROR_RAM);
        }
        result = 1;
    } 
    else 
    {
        if (*step_ct > *act_step)
        {
            if (log > 1)
            {
                fprintf(stdout, "Found %llu fluxtopes for step %d\n", *act_size, *act_step);
            }
            *act_step = *step_ct;
            *act_size = 1;
        }
        else
        {
            *act_size = *act_size + 1;
        }
        if (*step_ct > step)
        {
            *end = 1;
        }
    }
    return result;
}

void readFluxtopes(int ncol, unsigned long long* nrow, char*** ft, char*
        filename, int step, int** rx_map, int map_len, int all_rx_ct, char*
        rev, int rev_len, int log)
{
    char** m_ft = NULL;
    char* line = NULL;
    size_t len = 0;
    int bv_size = getBitsize(ncol);
    int rev_bv_size = getBitsize(rev_len);
    unsigned long ft_ct = 0;
    FILE* file;
    // store compressed positions of reactions
    int* positions = calloc(all_rx_ct, sizeof(int));
    int i;
    for (i = 0; i < all_rx_ct; i++) {
        positions[i] = getBitVectorPosition(rx_map, map_len, i);
    }
    // read file
    open_file(&file, filename, "r");
    int act_step = 0;
    unsigned long long act_size = 0;
    int started = 0;
    int parse = 0;
    int end = 0;
    int format = OUT_TXT;
    while (getline(&line, &len, file) != -1 && !end )
    {
        if (!started)
        {
            format = (strncmp(line, "compressedfluxtopes:", 20) == 0) ? OUT_BIN : OUT_TXT;
            if (format == OUT_BIN)
            {
                getline(&line, &len, file);
                getline(&line, &len, file);
                getline(&line, &len, file);
                char buffer[rev_len];
                while ( (fread( &buffer, 1, rev_bv_size, file) >= rev_bv_size) && !end ){
                    int step_ct = 0;
                    char* bv = calloc(1, bv_size);
                    int t_ix = -1;
                    for (i = 0; i < ncol; i++){
                        if (BITTEST(rev, i))
                        {
                            t_ix++;
                            if (BITTEST(buffer, t_ix))
                            {
                                BITSET(bv, i);
                            }
                            else
                            {
                                BITCLEAR(bv, i);
                                step_ct++;
                            }
                        }
                        else
                        {
                            BITSET(bv, i);
                        }
                    }
                    if (addToFluxtopes(&m_ft, &ft_ct, step, &act_step,
                                &step_ct, &act_size, &end, log))
                    {
                         m_ft[ft_ct - 1] = bv;
                    }
                    else 
                    {
                        free(bv);
                    }
                }
                started = 1;
            }
            else
            {
                parse = ((strncmp(line, "1,", 2) == 0) || (strncmp(line, "-1,", 3) == 0)) ? 1 : 0;
                started = 1;
            }
        }
        if (format == OUT_TXT)
        {
            if (parse)
            {
                char* bv = calloc(1, bv_size);
                int step_ct = 0;
                char* ptr = strtok(line, ",\n");
                int orig_rx_ix = 0;
                while(ptr != NULL)
                {
                    int pos = positions[orig_rx_ix];
                    if (pos > -1)
                    {
                        if (atoi(ptr) > 0)
                        {
                            BITSET(bv, pos);
                        } 
                        else
                        {
                            BITCLEAR(bv, pos);
                            step_ct++;
                        }
                    }
                    orig_rx_ix++;
                    ptr = strtok(NULL, ",\n");
                }
                if (addToFluxtopes(&m_ft, &ft_ct, step, &act_step,
                            &step_ct, &act_size, &end, log))
                {
                    m_ft[ft_ct - 1] = bv;
                }
                else 
                {
                    free(bv);
                }
            }
            else
            {
                parse = 1;
            }
        }
    }
    free(positions);
    positions = NULL;
    free(line);
    line = NULL;
    fclose(file);
    *nrow = ft_ct;
    *ft = m_ft;
}

void handle_arguments(int argc, char *argv[], char **ftout, int *outformat,
        char **kernel, char **rmlist, char **sfile, char **rfile, char
        **rvfile, char **mfile, double *eprhs, int *threads, char **lp_out, int
        *net, int *log, double *min_flux, double *max_flux, char **configfile,
        char** ft_in, int *ft_start, int *ft_end)
{
    //==================================================
    // define arguments and usage
    char *optv[MAX_ARGS] = { "-o", "--format", "-k", "--list", "-s", "-r",
        "-v", "-m", "-e", "-t", "-l", "--log", "--min", "--max", "-n",
        "--config", "--ftin", "--start", "--end"};
    char *optd[MAX_ARGS] = { 
        "output file [default: fluxtopes.out]",
        "format [txt,bin; default: bin]",
        "kernel matrix file",
        "removed list file",
        "stoichiometric file",
        "reaction file",
        "reversibility file",
        "metabolite file",
        "eprhs [default: 1e-6]",
        "number of threads [default: 1]",
        "lp output file [optional]",
        "log level (0-2) [default: 1]",
        "minimum flux [default: 1e-6]",
        "maximum flux [default: 1e3]",
        "use network LP (yes/no) [default: yes]",
        "config file [optional];\n\t\tl = lower bound; u = upper bound;"
            "\n\t\texample:\n\t\t  Reaction_id;l;0\n\t\t  Reaction_id;u;50",
        "input fluxtope file [optional, file will not be validated!]",
        "start step [optional, needs --ftin]",
        "end step [optional]"
    };
    char *optr[MAX_ARGS];
    char *description = "Calculates flux topes of an metabolic network";
    char *usg = "\n   fluxTopeEnumerator -o fluxtopes.out -k kernel --list "
        "rmlist -s sfile -r rfile -v rvfile -m mfile -e 1e-6 -t 8 -l debug.lp "
        "-n no --min 1e-5 --max 1e4 --ftin fluxtopes.out.1 --start 4";
    // end define arguments and usage
    //==================================================

    //==================================================
    // read arguments
    readArgs(argc, argv, MAX_ARGS, optv, optr);
    if ( !optr[ARG_KERNEL] || !optr[ARG_RM_LIST] || !optr[ARG_SFILE] ||
            !optr[ARG_RFILE] || !optr[ARG_RVFILE] || !optr[ARG_MFILE] )
    {
        usage(description, usg, MAX_ARGS, optv, optd);
        quitError("Missing argument\n", ERROR_ARGS);
    }
    if ( optr[ARG_FT_START] && !optr[ARG_FT_IN] )
    {
        usage(description, usg, MAX_ARGS, optv, optd);
        quitError("--start needs --ftin argument\n", ERROR_ARGS);
    }

    *ftout      = optr[ARG_FT_OUT] ? optr[ARG_FT_OUT] : "fluxtopes.out";
    *outformat  = OUT_BIN;
    if (optr[ARG_OUTFORMAT])
    {
        if (strcmp(optr[ARG_OUTFORMAT], "txt") == 0){
            *outformat = OUT_TXT;
        }
    }
    *kernel     = optr[ARG_KERNEL];
    *rmlist     = optr[ARG_RM_LIST];
    *sfile      = optr[ARG_SFILE];
    *rfile      = optr[ARG_RFILE];
    *rvfile     = optr[ARG_RVFILE];
    *mfile      = optr[ARG_MFILE];
    *min_flux   = optr[ARG_MINFLUX] ? atof(optr[ARG_MINFLUX]) : 1e-6;
    *max_flux   = optr[ARG_MAXFLUX] ? atof(optr[ARG_MAXFLUX]) : 1e3;
    *eprhs      = optr[ARG_ZERO] ? atof(optr[ARG_ZERO]) : 1e-6;
    *threads    = optr[ARG_THREADS] ? atoi(optr[ARG_THREADS]) : 1;
    *log        = optr[ARG_LOG] ? atoi(optr[ARG_LOG]) : 1;
    if (*log < 0 || *log > 2)
    {
        *log = 1;
    }
    *configfile = optr[ARG_CONFIGFILE] ? optr[ARG_CONFIGFILE] :
        "not_available";
    *lp_out     = optr[ARG_LP_OUT] ? optr[ARG_LP_OUT] : "not_available";
    char *s_net = "yes";
    s_net = optr[ARG_NET] ? optr[ARG_NET] : "yes";
    if (strcmp(s_net, "no") == 0)
    {
        *net = 0;
    }
    else
    {
        *net = 1;
    }
    *ft_start = optr[ARG_FT_START] ? atoi(optr[ARG_FT_START]) : 0;
    *ft_end = optr[ARG_FT_END] ? atoi(optr[ARG_FT_END]) : -1;
    if (optr[ARG_FT_IN])
    {
        *ft_in = optr[ARG_FT_IN];
    }
    // end read arguments
    //==================================================

    //==================================================
    // print arguments summary
    if (*log > 0)
    {
        fprintf(stdout, "\n");
        fprintf(stdout, "Output:             %s\n", *ftout);
        char *temp = "bin";
        if (*outformat == OUT_TXT)
        {
            temp = "txt";
        }
        fprintf(stdout, "Format:             %s\n", temp);
        fprintf(stdout, "kernel:             %s\n", *kernel);
        fprintf(stdout, "rmlist file:        %s\n", *rmlist);
        fprintf(stdout, "sfile:              %s\n", *sfile);
        fprintf(stdout, "rfile:              %s\n", *rfile);
        fprintf(stdout, "rvfile:             %s\n", *rvfile);
        fprintf(stdout, "mfile:              %s\n", *mfile);
        fprintf(stdout, "eprhs:              %.2e\n", *eprhs);
        fprintf(stdout, "Threads:            %d\n", *threads);
        fprintf(stdout, "LP outfile:         %s\n", *lp_out);
        fprintf(stdout, "Use network LP:     %s\n", s_net);
        fprintf(stdout, "Log level:          %d\n", *log);
        fprintf(stdout, "minimum flux:       %.2e\n", *min_flux);
        fprintf(stdout, "maximum flux:       %.2e\n", *max_flux);
        fprintf(stdout, "config file:        %s\n", *configfile);
        if (*ft_start > 0)
        {
            fprintf(stdout, "Fluxtope inputfile: %s\n", *ft_in);
            fprintf(stdout, "Start step:         %d\n", *ft_start);
        }
        if (*ft_end > -1)
        {
            fprintf(stdout, "End step:           %d\n", *ft_end);
        }
        fprintf(stdout, "\n");
    }
    // end print arguments summary
    //==================================================
}

int get_col_count(char *filename)
{
    FILE *file;
    open_file(&file, filename, "r");
    int col_count = getRxCount(file);
    fclose(file);
    if (col_count < 1)
    {
        quitError("Error in file format; number of columns < 1\n",
                ERROR_FILE);
    }
    return col_count;
}

void createReactionMapping(int *rmlist, char *reversible, int
        ***rx_map, char **kern_rev, int *map_size, int **map_lengths, int
        rx_count, int bitsize)
{
    int i;
    int ix = 0;
    int **map = calloc(rx_count, sizeof(int**));
    int *maps = calloc(rx_count, sizeof(int*));
    char *rev = calloc(1, bitsize);
    for (i = 0; i < rx_count; i++)
    {
        if (rmlist[i] < 0)
        {
            map[ix] = calloc(sizeof(int*), rx_count);
            map[ix][0] = i;
            int counter = 1;
            int j;
            for (j = 0; j < rx_count; j++)
            {
                if (rmlist[j] == i)
                {
                    map[ix][counter] = j;
                    counter++;
                }
            }
            map[ix] = (int*) realloc(map[ix], sizeof(int) * counter);
            maps[ix] = counter;
            int t_rev = 1;
            for (j = 0; j < counter; j++)
            {
                if (!BITTEST(reversible, map[ix][j]))
                {
                    t_rev = 0;
                    j = counter;
                }
            }
            if (t_rev)
            {
                BITSET(rev, ix);
            }
            else
            {
                BITCLEAR(rev, ix);
            }
            ix++;
        }
    }
    *rx_map = map;
    *kern_rev = rev;
    *map_size = ix;
    *map_lengths = maps;
}

void init_cplex(CPXENVptr *env, double eprhs)
{
    int status = 0;
    CPXENVptr m_env = CPXopenCPLEX (&status);
    if ( m_env == NULL )
    {
        char  errmsg[CPXMESSAGEBUFSIZE];
        fprintf (stderr, "Could not open CPLEX environment.\n");
        CPXgeterrorstring (m_env, status, errmsg);
        fprintf (stderr, "%s", errmsg);
        exit(ERROR_CPLEX);
    }
    status = CPXsetdblparam(m_env, CPX_PARAM_EPRHS, eprhs);
    *env = m_env;
}

void createNetworkLpProblem (CPXENVptr env, CPXLPptr *lp, double **stoich, int
        s_rows, int s_cols, struct metabolite *metabolites, 
        char** rx_names, double* lower_bounds, double* upper_bounds, double
        min_flux, char* name)
{
    int status = 0;
    // create LP problem
    CPXLPptr m_lp = CPXcreateprob (env, &status, name);
    if ( m_lp == NULL )
    {
        fprintf (stderr, "Failed to create LP.\n");
        exit(ERROR_CPLEX);
    }
    // LP Columns
    // ==================================================
    int col_counter = 0;
    int i;
    // allocate bound arrays (rx, DrG, conc, DfG)
    int col_len = 2 * s_cols + 2 * s_rows;
    double lb[col_len];
    double ub[col_len];
    char **colnames = calloc(col_len, sizeof(char*));
    // define columns for reactions
    for (i = 0; i < s_cols; i++)
    {
        lb[col_counter] = lower_bounds[col_counter] > min_flux ?
            lower_bounds[col_counter] : min_flux;
        ub[col_counter] = upper_bounds[col_counter];
        colnames[col_counter] = calloc(strlen(rx_names[i]) + 1, sizeof(char));
        strcpy(colnames[col_counter], rx_names[i]);
        col_counter++;
    }
    // add columns to LP
    status = CPXnewcols(env, m_lp, col_counter, NULL, lb, ub, NULL, colnames);
    // LP Rows
    // ==================================================
    // define non zero count
    // ---------------------
    int numNewNonZero = 0;
    // reaction constraints
    for (i = 0; i < s_rows; i++)
    {
        int j;
        for (j = 0; j < s_cols; j++)
        {
           if (stoich[i][j] != 0.0)
           {
               numNewNonZero++;
           }
        }
    }
    int use_rx = 0;
    int constraint_count = s_rows + use_rx;
    double rhs[constraint_count];
    char sense[constraint_count];
    int rmatbeg[constraint_count];
    char** rownames = calloc(constraint_count, sizeof(char*));
    int rmatidx[numNewNonZero];
    double rmatval[numNewNonZero];
    // reaction constraints
    int nonZeroCnt = 0;
    for (i = 0; i < s_rows; i++)
    {
        rownames[i] = calloc(strlen(metabolites[i].name) + 1, sizeof(char));
        strcpy(rownames[i], metabolites[i].name);
        rhs[i] = 0;
        sense[i] = 'E';
        rmatbeg[i] = nonZeroCnt;
        int j;
        for (j = 0; j < s_cols; j++)
        {
            if (stoich[i][j] != 0.0)
            {
                rmatidx[nonZeroCnt] = j;
                rmatval[nonZeroCnt] = stoich[i][j];
                nonZeroCnt++;
            }
        }
    }
    int ix = s_rows - 1;
    // add constraints to LP
    status = CPXaddrows(env, m_lp, 0, constraint_count, numNewNonZero, rhs,
            sense, rmatbeg, rmatidx, rmatval, NULL, rownames);
    for (i = 0; i < col_counter; i++)
    {
        free(colnames[i]);
    }
    free(colnames);
    for (i = 0; i <= ix; i++)
    {
        free(rownames[i]);
    }
    free(rownames);
    *lp = m_lp;
}

void createKernelLpProblem (CPXENVptr env, CPXLPptr *lp, double **kernel, int
        k_rows, int k_cols, char* name)
{
    int status = 0;
    CPXLPptr m_lp = CPXcreateprob (env, &status, name);
    if ( m_lp == NULL )
    {
        fprintf (stderr, "Failed to create LP.\n");
        exit(ERROR_CPLEX);
    }
    int col_counter = 0;
    int i;
    double lb[k_cols + k_rows];
    double ub[k_cols + k_rows];
    for (i = 0; i < k_cols; i++)
    {
        lb[i] = -1000;
        ub[i] = 1000;
        col_counter++;
    }
    for (i = 0; i < k_rows; i++)
    {
        lb[col_counter] = 1e-8;
        ub[col_counter] = 1000;
        col_counter++;
    }
    status = CPXnewcols(env, m_lp, col_counter, NULL, lb, ub, NULL, NULL);
    double rhs[k_rows];
    char sense[k_rows];
    int numNewNonZero = 0;
    for (i = 0; i < k_rows; i++)
    {
        rhs[i] = 0;
        sense[i] = 'E';
        int j;
        for (j = 0; j < k_cols; j++)
        {
           if (kernel[i][j] != 0.0)
           {
               numNewNonZero++;
           }
        }
        numNewNonZero++;
    }
    int rmatbeg[k_rows + 1];
    int rmatidx[numNewNonZero];
    double rmatval[numNewNonZero];
    int nonZeroCnt = 0;
    for (i = 0; i < k_rows; i++)
    {
        rmatbeg[i] = nonZeroCnt;
        int j;
        for (j = 0; j < k_cols; j++)
        {
            if (kernel[i][j] != 0.0)
            {
                rmatidx[nonZeroCnt] = j;
                rmatval[nonZeroCnt] = kernel[i][j];
                nonZeroCnt++;
            }
        }
        rmatidx[nonZeroCnt] = k_cols + i;
        rmatval[nonZeroCnt] = -1;
        nonZeroCnt++;
    }
    rmatbeg[i] = nonZeroCnt;
    status = CPXaddrows(env, m_lp, 0, k_rows, numNewNonZero, rhs, sense,
            rmatbeg, rmatidx, rmatval, NULL, NULL);
    *lp = m_lp;
}

int checkFeasibility(CPXENVptr env, CPXLPptr lp)
{
    int status = CPXlpopt (env, lp);
    if (status != 0)
    {
        return status;
    }
    int lpstat = CPXgetstat (env, lp);
    if (lpstat == CPX_STAT_OPTIMAL)
    {
        return 0;
    }
    return lpstat;
}

int setNetworkBounds(char *bv, int bv_size, int ix, int **map, int *maps, int
        all_rx_len, CPXENVptr env, CPXLPptr lp, double *lb, double *ub, double
        min_flux)
{
    int m_size = all_rx_len;
    double m_lb[m_size];
    double m_ub[m_size];
    int ind[m_size];
    char llu[m_size];
    char lub[m_size];
    int i;
    int failed = 0;
    for (i = 0; i < bv_size; i++)
    {
        int j;
        for (j = 0; j < maps[i]; j++)
        {
            int constr = map[i][j];
            ind[constr] = constr;
            llu[constr] = 'L';
            lub[constr] = 'U';
            if (i == ix || !BITTEST(bv, i))
            {
                m_lb[constr] = lb[constr]; 
                m_ub[constr] = -1 * min_flux; 
                if (m_lb[constr] > m_ub[constr])
                {
                    failed = 1;
                }
            }
            else
            {
                m_lb[constr] = lb[constr] > min_flux ? lb[constr] : min_flux;
                m_ub[constr] = ub[constr];
                if (m_ub[constr] < m_lb[constr])
                {
                    failed = 1;
                }
            }
        }
    }
    int status = CPXchgbds(env, lp, m_size, ind, llu, m_lb);
    if (status != 0)
    {
        fprintf(stderr, 
                "Could not change lower bound of LP problem. status = %d\n",
                status);
        exit(ERROR_CPLEX);
    }
    status = CPXchgbds(env, lp, m_size, ind, lub, m_ub);
    if (status != 0)
    {
        fprintf(stderr, 
                "Could not change upper bound of LP problem. status = %d\n",
                status);
        exit(ERROR_CPLEX);
    }
    return failed;
}

void setKernelBounds(char *bv, int bv_size, int ix, int cols, CPXENVptr env,
        CPXLPptr lp)
{
    double lb[bv_size];
    double ub[bv_size];
    int ind[bv_size];
    char llu[bv_size];
    char lub[bv_size];
    int i;
    for (i = 0; i < bv_size; i++)
    {
        ind[i] = i + cols;
        llu[i] = 'L';
        lub[i] = 'U';
        if (i == ix || !BITTEST(bv, i))
        {
            lb[i] = -1000;
            ub[i] = -1e-8;
        }
        else
        {
            lb[i] = 1e-8;
            ub[i] = 1000;
        }
    }
    int status = CPXchgbds(env, lp, bv_size, ind, llu, lb);
    if (status != 0)
    {
        fprintf(stderr, 
                "Could not change lower bound of LP problem. status = %d\n",
                status);
        exit(ERROR_CPLEX);
    }
    status = CPXchgbds(env, lp, bv_size, ind, lub, ub);
    if (status != 0)
    {
        fprintf(stderr, 
                "Could not change upper bound of LP problem. status = %d\n",
                status);
        exit(ERROR_CPLEX);
    }
}

void* findAdj(void *pointer_thread_args)
{
    struct thread_args *thread_args = (struct thread_args*)
        pointer_thread_args;
    char *bv       = thread_args->bv;
    int bv_size    = thread_args->bv_size;
    char *rev      = thread_args->rev;
    int k_rev_len  = thread_args->k_rev_len;
    int cols       = thread_args->cols;
    char ***adj_bv = thread_args->adj_bv;
    int *adj_count = thread_args->adj_count;
    int *finished  = thread_args->finished;
    int net        = thread_args->net;
    int **map      = thread_args->map;
    int *maps      = thread_args->maps;
    int all_rx_len = thread_args->all_rx_len;
    double *lb     = thread_args->lb;
    double *ub     = thread_args->ub;
    double minflux = thread_args->min_flux;
    CPXENVptr env  = thread_args->env;
    CPXLPptr lp    = thread_args->lp;
    int m_adj_ct = 0;
    char **m_adj_bv = NULL;
    int i;
    for (i = 0; i < bv_size; i++)
    {
        if (BITTEST(rev, i) && BITTEST(bv,i))
        {
            char *nv = calloc(1, getBitsize(k_rev_len));
            if (nv == NULL)
            {
                fprintf(stderr, "Not enough RAM left\n");
                exit(ERROR_RAM);
            }
            int j;
            int t_ix = -1;
            for (j = 0; j < bv_size; j++)
            {
                if (j == i)
                {
                    t_ix++;
                }
                else if (BITTEST(rev, j))
                {
                    t_ix++;
                    if (BITTEST(bv, j))
                    {
                        BITSET(nv, t_ix);
                    }
                }
            }
            int tree_check = -1;
            while (tree_check < 0)
            {
                tree_check = tree_insert(nv, &root, k_rev_len);
            }
            if (tree_check > 0)
            {
                int feasible = 0;
                if (net)
                {
                    feasible = setNetworkBounds(bv, bv_size, i, map, maps, all_rx_len,
                            env, lp, lb, ub, minflux);
                }
                else
                {
                    setKernelBounds(bv, bv_size, i, cols, env, lp);
                }
                if (feasible == 0)
                {
                    feasible = checkFeasibility(env, lp);
                }
                if (feasible == 0)
                {
                    char *n_bv = calloc(1, getBitsize(bv_size));
                    if (n_bv == NULL)
                    {
                        fprintf(stderr, "Not enough RAM left\n");
                        exit(ERROR_RAM);
                    }
                    int j;
                    for (j = 0; j < bv_size; j++)
                    {
                        if (j != i)
                        {
                            if (BITTEST(bv, j))
                            {
                                BITSET(n_bv, j);
                            }
                        }
                    }
                    m_adj_ct++;
                    m_adj_bv = (char**) realloc(m_adj_bv, m_adj_ct *
                            sizeof(char*));
                    if (m_adj_bv == NULL)
                    {
                        fprintf(stderr, "Not enough RAM left\n");
                        exit(ERROR_RAM);
                    }
                    m_adj_bv[m_adj_ct - 1] = n_bv;
                }
            }
            else
            {
                free(nv);
            }
        }
    }
    free(bv);
    *adj_count = m_adj_ct;
    if (m_adj_ct > 0)
    {
        *adj_bv = m_adj_bv;
    }
    *finished = 1;
    return((void*)NULL);
}

int allThreadsFree(int threads, int *thread_lock)
{
    int i;
    for (i = 0; i < threads; i++)
    {
        if (thread_lock[i] > 0)
        {
            return 0;
        }
    }
    return 1;
}

int getFreeThread(int threads, int *thread_lock)
{
    int i;
    for (i = 0; i < threads; i++)
    {
        if (thread_lock[i] < 1)
        {
            return i;
        }
    }
    return -1;
}

int getFinishedThread(int threads, int *thread_finished)
{
    if (threads < 2)
    {
        usleep(1);
    }
    int i;
    for (i = 0; i < threads; i++)
    {
        if (thread_finished[i] > 0)
        {
            return i;
        }
    }
    return -1;
}

void findFluxTopes(unsigned long long *foundCounter, char** bv, unsigned long
        long bv_len, int bv_size, char* rev, int k_rev_len, int cols, int
        threads, CPXENVptr *envs, CPXLPptr *lps, int **map, int *maps, int
        all_rx_len, double min_flux, int net, FILE *f, int log, double
        *lower_bounds, double *upper_bounds, int step, int end_step, int
        outformat)
{
    int *thread_lock                = calloc(threads, sizeof(int));
    int *thread_finished            = calloc(threads, sizeof(int));
    int *adj_ct                     = calloc(threads, sizeof(int));
    char ***adj_bv                  = calloc(threads, sizeof(char**));
    char **m_adj_bv                 = NULL;
    unsigned long long m_adj_bv_len = 0;
    struct thread_args *thread_args = calloc(threads, 
            sizeof(struct thread_args));
    pthread_t thread[threads];
    int alloc_size = getBitsize(bv_size);
    int i;
    // initialize variables
    for (i = 0; i < threads; i++)
    {
        thread_lock[i] = 0;
        thread_finished[i] = 0;
        adj_ct[i] = 0;
        adj_bv[i] = NULL;
        thread_args[i].net        = net;
        thread_args[i].bv_size    = bv_size;
        thread_args[i].rev        = malloc(getBitsize(bv_size));
        copyBitvector(rev, thread_args[i].rev, bv_size);
        thread_args[i].k_rev_len  = k_rev_len;
        thread_args[i].cols       = cols;
        thread_args[i].thread_id  = i;
        thread_args[i].adj_count  = &adj_ct[i];
        thread_args[i].adj_bv     = &adj_bv[i];
        thread_args[i].finished   = &thread_finished[i];
        thread_args[i].all_rx_len = all_rx_len;
        thread_args[i].map        = map;
        thread_args[i].maps       = maps;
        thread_args[i].env        = envs[i];
        thread_args[i].lp         = lps[i];
        thread_args[i].lb         = lower_bounds;
        thread_args[i].ub         = upper_bounds;
        thread_args[i].min_flux   = min_flux;
    }
    if (log > 0)
    {
        fprintf(stdout,
                "--------------------------------------\n");
        unsigned long long possibles = bv_len * (k_rev_len - step);
        fprintf(stdout, "Step: %d, possible flux topes: %llu\n", step,
                possibles);
        fflush(stdout);
    }
    if (step == 1)
    {
        (*foundCounter)++;
        printVector(f, bv[0], bv_size, map, maps, all_rx_len, outformat, rev,
                k_rev_len);
    }
    // as long as adjacent flux topes need to be searched
    int run = 1;
    while (run) 
    {
        if (bv_len < 1)
        {
            if (allThreadsFree(threads, thread_lock))
            {
                if (m_adj_bv_len > 0)
                {
                    fflush(f);
                    if (log > 1)
                    {
                        fprintf(stdout, "%s   "
                                "Found flux topes: %llu, "
                                "tree size: %llu, "
                                "found in current tree: %llu, "
                                "flux topes for this step: %llu, "
                                "flux topes for next step: %llu, "
                                "all tree entries: %llu, "
                                "all found in tree: %llu\n",
                                getTime(), (*foundCounter), act_tree_size,
                                act_tree_exist, bv_len, m_adj_bv_len,
                                tree_size, tree_exist);
                        fprintf(stdout,
                                "--------------------------------------\n");
                    }
                    step++;
                    if (end_step > -1)
                    {
                        if (step > end_step)
                        {
                            run = 0;
                            int abi;
                            for (abi = 0; abi < m_adj_bv_len; abi++) {
                                free(m_adj_bv[abi]);
                            }
                        }
                    }
                    if (run)
                    {
                        free(bv);
                        bv = m_adj_bv;
                        m_adj_bv = NULL;
                        bv_len = m_adj_bv_len;
                        m_adj_bv_len = 0;
                        tree_destroy(root);
                        act_tree_exist = 0;
                        act_tree_size = 0;
                        root = 0;
                        if (log > 0)
                        {
                            unsigned long long possibles = 
                                bv_len * (k_rev_len - step);
                            fprintf(stdout, 
                                    "Step: %d, possible flux topes: %llu\n", step,
                                    possibles);
                            fflush(stdout);
                        }
                    }
                }
                else
                {
                    run = 0;
                }
            }
        }
        int free_thread = -1;
        // join finished threads and check free thread
        while (free_thread < 0)
        {
            int finished = getFinishedThread(threads, thread_finished);
            if (finished > -1)
            {
                pthread_join(thread[finished], NULL);
                if (adj_ct[finished] > 0)
                {
                    unsigned long long old_last = m_adj_bv_len - 1;
                    m_adj_bv_len += adj_ct[finished];
                    m_adj_bv = (char**) realloc(m_adj_bv, m_adj_bv_len *
                            sizeof(char*));
                    if (m_adj_bv == NULL)
                    {
                        fprintf(stderr, "Not enough RAM left\n");
                        exit(ERROR_RAM);
                    }
                    for (i = 0; i < adj_ct[finished]; i++)
                    {
                        old_last++;
                        m_adj_bv[old_last] = malloc(alloc_size);
                        if (m_adj_bv[old_last] == NULL)
                        {
                            fprintf(stderr, "Not enough RAM left\n");
                            exit(ERROR_RAM);
                        }
                        (*foundCounter)++;
                        if ((*foundCounter)%100 == 0)
                        {
                            fflush(f);
                            if (log > 1)
                            {
                                fprintf(stdout, "%s   "
                                        "Found flux topes: %llu, "
                                        "tree size: %llu, "
                                        "found in current tree: %llu, "
                                        "flux topes for this step: %llu, "
                                        "flux topes for next step: %llu, "
                                        "all tree entries: %llu, "
                                        "all found in tree: %llu\n",
                                        getTime(), (*foundCounter),
                                        act_tree_size, act_tree_exist, bv_len,
                                        m_adj_bv_len, tree_size, tree_exist);
                                fflush(stdout);
                            }
                        }
                        printVector(f, adj_bv[finished][i], bv_size, map, maps,
                                all_rx_len, outformat, rev, k_rev_len);
                        copyBitvector(adj_bv[finished][i], m_adj_bv[old_last],
                                bv_size);
                        free(adj_bv[finished][i]);
                    }
                    free(adj_bv[finished]);
                }
                thread_lock[finished] = 0;
                thread_finished[finished] = 0;
                free_thread = finished;
            }
            else
            {
                free_thread = getFreeThread(threads, thread_lock);
            }
        }
        if (free_thread > -1 && bv_len > 0)
        {
            thread_lock[free_thread] = 1;
            unsigned long long arr_ix = bv_len - 1;
            thread_args[free_thread].bv = bv[arr_ix];
            int status = pthread_create(&thread[free_thread], NULL,
                    findAdj, (void*)&thread_args[free_thread]);
            if (status != 0)
            {
                fprintf(stderr, "ERROR: thread status (%d) != 0\n", status);
            }
            bv_len--;
        }
    }
    for (i = 0; i < threads; i++)
    {
        free(thread_args[i].rev);
    }
    free(thread_lock);
    free(thread_finished);
    free(thread_args);
    free(adj_ct);
    free(adj_bv);
    free(m_adj_bv);
    free(bv);
}

void defineMetabolites(struct metabolite *metabolites, char **met_names, int
        met_count) 
{
    int i;
    for (i = 0; i < met_count; i++)
    {
        metabolites[i].name = met_names[i];
    }
}

void defineBounds(double **lower_bounds, double **upper_bounds, int s_cols,
        char** reaction_names, char* reversibles, double min_flux, double
        max_flux, char* configfile)
{
    int i;
    double *m_lb = calloc(s_cols, sizeof(double));
    double *m_ub = calloc(s_cols, sizeof(double));
    for (i = 0; i < s_cols; i++)
    {
        if (BITTEST(reversibles, i))
        {
            m_lb[i] = -1 * max_flux;
        }
        else
        {
            m_lb[i] = min_flux;
        }
        m_ub[i] = max_flux;
    }
    if (strcmp(configfile, "not_available") != 0)
    {
        char *line = NULL;
        size_t len = 0;
        FILE *file;
        open_file(&file, configfile, "r");
        while (getline(&line, &len, file) != -1)
        {
            char *ptr;
            ptr = strtok(line, "\n ");
            if (ptr != NULL)
            {
                int index = -1;
                double value = -1;
                int is_lower = -1;
                char *ptr2;
                ptr2 = strtok(ptr, ";");
                int c = 0;
                while (ptr2 != NULL)
                {
                    c++;
                    if (c == 1)
                    {
                        int j;
                        for (j = 0; j < s_cols; j++)
                        {
                            if (strcmp(reaction_names[j], ptr2) == 0)
                            {
                                index = j;
                                j = s_cols;
                            }
                        }
                    }
                    else if (c == 2)
                    {
                        if (strcmp(ptr2, "l") == 0)
                        {
                            is_lower = 1;
                        }
                        else if (strcmp(ptr2, "u") == 0)
                        {
                            is_lower = 0;
                        }
                    }
                    else if (c == 3)
                    {
                        value = atof(ptr2);
                    }
                    ptr2 = strtok(NULL, ";");
                }
                if (index >= 0 && index < s_cols && value >= 0 && 
                        is_lower >= 0)
                {
                    if (is_lower)
                    {
                        m_lb[index] = value;
                    }
                    else
                    {
                        m_ub[index] = value;
                    }
                }
            }
        }
        free(line);
        line = NULL;
        fclose(file);
    }
    *lower_bounds = m_lb;
    *upper_bounds = m_ub;
}


void quitError ( char *message, int rv )
{
    fprintf( stderr, "%s", message );
    exit(rv);
}

void quitErrorArg ( char *message, char *arg, int rv )
{
    fprintf( stderr, "%s %s\n", message, arg );
    exit(rv);
}

int getRxCount(FILE* file)
{
    char*  line = NULL;
    size_t len  = 0;
    int rx_count = 0;
    if ( getline(&line, &len, file) != -1)
    {
        char *ptr;
        ptr = strtok(line, "\n\t ");
        while(ptr != NULL)
        {
            rx_count++;
            ptr = strtok(NULL, "\n\t ");
        }
    }
    free(line);
    line = NULL;
    return rx_count;
}

double binom(int n, int k)
{
    if (k < 0)
    {
        return 0;
    }
    double* C = (double*)calloc(k+1, sizeof(double));
    int i, j;
    double res;
    C[0] = 1;
    for(i = 1; i <= n; i++)
    {
        for(j = min(i, k); j > 0; j--)
        {
            C[j] += C[j-1];
        }
    }
    res = C[k];  // Store the result before freeing memory
    free(C);  // free dynamically allocated memory to avoid memory leak
    return res;
}

int min(int a, int b)
{
    return (a < b) ? a : b;
}

void readArgs ( int argc, char *argv[], int optc, char *optv[], char *optr[] )
{
    unsigned int i;
    for (i = 0; i < optc; i++)
    {
        optr[i] = getArg(argc, argv, optv[i]);
    }
}

char* getArg ( int argc, char *argv[], char *opt )
{
    char *ret = NULL;
    unsigned int i;
    for (i = 0; i < (argc-1); i++)
    {
        if (!strcmp(argv[i], opt))
        {
            ret = argv[i+1];
            break;
        }
    }
    return ret;
}

void usage (char* description, char* usage, int optc, char* optv[], 
        char *optd[])
{
    printf("\n%s\n", description);
    printf("\nusage: %s\n", usage);
    printf("arguments:\n");
    unsigned int i;
    for (i = 0; i < optc; i++)
    {
        printf("   %s\t%s\n", optv[i], optd[i]);
    }
    printf("\n");
}
