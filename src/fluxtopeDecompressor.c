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

#include "fluxtopeDecompressor.h"

int main (int argc, char *argv[])
{

    char *ft_in, *ft_out, *sep, *fwd_str, *rev_str;
    int ft_start, ft_end, use_sep, log;

    handle_arguments(argc, argv, &ft_in, &ft_out, &ft_start, &ft_end, &use_sep,
            &sep, &fwd_str, &rev_str, &log);

    if (log > 0)
    {
        fprintf(stdout, "Start: %s\n", getTime());
    }
    decompressFluxtopes(ft_in, ft_out, use_sep, sep, ft_start, ft_end, fwd_str,
            rev_str, log);
    return EXIT_SUCCESS;
}


void decompressFluxtopes(char* in, char* out, int use_sep, char* sep, int
        ft_start, int ft_end, char* fwd_str, char* rev_str, int log)
{
    FILE *ftout, *ftin;
    open_file(&ftout, out, "w");
    open_file(&ftin, in, "r");
    char* line = NULL;
    size_t len = 0;
    int reactions = 0;
    int kernel_len = 0;
    int map_size;
    int *map_lengths;
    int **rx_map;
    char *rev;
    int rev_ct = 0;
    // read first line
    int linelength = getline(&line, &len, ftin);
    parseFirstLine(linelength, line, &reactions, &kernel_len);
    // define end step
    if (ft_end < 0)
    {
        ft_end = reactions + 1;
    }
    // read second line (reaction map)
    linelength = getline(&line, &len, ftin);
    createReactionMapping(linelength, line, reactions, &rx_map, &map_size,
            &map_lengths); 
    // read third line (reversible kernel reactions)
    linelength = getline(&line, &len, ftin);
    parseReversibles(linelength, line, kernel_len, &rev, &rev_ct);
    // read fourth line (reaction names)
    linelength = getline(&line, &len, ftin);
    printReactionNames(linelength, ftout, line, reactions);
    int bitsize = getBitsize(reactions);
    int bv_size = getBitsize(map_size);
    int rev_size = getBitsize(rev_ct);
    char buffer[kernel_len];
    int end = 0;
    int act_step = 0;
    unsigned long long act_size = 0;
    unsigned long long tot_size = 0;
    while ( (fread( &buffer, 1, rev_size, ftin) >= rev_size) && !end ){
        int step_ct = 0;
        char* bv = calloc(1, bv_size);
        int i;
        int t_ix = -1;
        for (i = 0; i < kernel_len; i++){
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
        int val = isPrintableFluxtope(&act_step, &step_ct, &act_size,
                &tot_size, ft_start, ft_end, log);
        if (val == 1)
        {
            printVector(ftout, bv, map_size, rx_map, map_lengths, reactions,
                    use_sep, sep, bitsize, fwd_str, rev_str);
        }
        else if (val == 2)
        {
            end = 1;
        }
        free(bv);
    }
    if (log > 0)
    {
        fprintf(stdout, "Total printed fluxtopes: %llu\n", tot_size);
    }
    fclose(ftout);
    fclose(ftin);
    free(line);
    free(rev);
    int i;
    for (i = 0; i < map_size; i++) {
        free(rx_map[i]);
    }
    free(rx_map);
    free(map_lengths);
}

char* getTime()
{
    static char time_string[20];
    time_t now = time (0);
    strftime (time_string, 100, "%Y-%m-%d %H:%M:%S", localtime (&now));
    return time_string;
}

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

void open_file(FILE **file, char *name, char *handle)
{
   *file = fopen(name, handle);
   if (*file == NULL)
   {
       quitErrorArg("Error in opening file: ", name, ERROR_FILE);
   }
}

void handle_arguments(int argc, char *argv[], char** ft_in, char **ft_out,
        int *ft_start, int *ft_end, int* use_sep, char** sep, char** fwd_str,
        char** rev_str, int *log)
{
    //==================================================
    // define arguments and usage
    char *optv[MAX_ARGS] = { "-i", "-o", "--start", "--end", "--sep", "--fwd",
        "--rev", "--log" };
    char *optd[MAX_ARGS] = { 
        "input compressed fluxtopefile",
        "output file",
        "start step [default: 0]",
        "end step [optional]",
        "separator [default: no separator]",
        "forward string [default: 1]",
        "reverse string [default: -1]",
        "log level (0-2) [default: 1]"
    };
    char *optr[MAX_ARGS];
    char *description = "Decompresses fluxtopes";
    char *usg = "\n   fluxtopeDecompressor -i fluxtopes.in -o fluxtopes.out "
        "--start 4 --end 4 --sep ; --fwd + --rev - --log 2";
    // end define arguments and usage
    //==================================================

    //==================================================
    // read arguments
    readArgs(argc, argv, MAX_ARGS, optv, optr);
    if ( !optr[ARG_FT_IN] || !optr[ARG_FT_OUT] )
    {
        usage(description, usg, MAX_ARGS, optv, optd);
        quitError("Missing argument\n", ERROR_ARGS);
    }

    *ft_in  = optr[ARG_FT_IN];
    *ft_out = optr[ARG_FT_OUT];
    *sep = ",";
    *use_sep = 0;
    if (optr[ARG_SEP])
    {
        *sep = optr[ARG_SEP];
        *use_sep = 1;
    }
    *fwd_str = optr[ARG_FWD] ? optr[ARG_FWD] : "1";
    *rev_str = optr[ARG_REV] ? optr[ARG_REV] : "-1";
    *log = optr[ARG_LOG] ? atoi(optr[ARG_LOG]) : 1;
    if (*log < 0 || *log > 2)
    {
        *log = 1;
    }
    *ft_start = optr[ARG_FT_START] ? atoi(optr[ARG_FT_START]) : 0;
    *ft_end = optr[ARG_FT_END] ? atoi(optr[ARG_FT_END]) : -1;
    // end read arguments
    //==================================================

    //==================================================
    // print arguments summary
    if (*log > 0)
    {
        fprintf(stdout, "\n");
        fprintf(stdout, "Input:              %s\n", *ft_in);
        fprintf(stdout, "Output:             %s\n", *ft_out);
        fprintf(stdout, "Separator:          ");
        if (*use_sep)
        {
            fprintf(stdout, "%s", *sep);
        }
        fprintf(stdout, "\n");
        fprintf(stdout, "Log level:          %d\n", *log);
        fprintf(stdout, "Start step:         %d\n", *ft_start);
        fprintf(stdout, "End step:           ");
        if (*ft_end > -1)
        {
            fprintf(stdout, "%d\n", *ft_end);
        }
        else
        {
            fprintf(stdout, "last\n");
        }
        fprintf(stdout, "Forward string:     %s\n", *fwd_str);
        fprintf(stdout, "Reverse string:     %s\n", *rev_str);
    }
    // end print arguments summary
    //==================================================
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

void parseReversibles(int linelength, char* line, int size, char** rev, int* rev_ct)
{
    if (linelength < 1)
    {
        quitError("Wrong format of input file: reversible kernel reactions are missing\n", 
                ERROR_FORMAT);
    }
    char* m_rev = calloc(1, getBitsize(size));
    int ct = 0;
    int ix = -1;
    char* ptr = strtok(line, ",\n");
    while(ptr != NULL)
    {
        ix++;
        int val = (int) strtol(ptr, NULL, 10);
        if (val == 1)
        {
            BITSET(m_rev, ix);
            ct++;
        }
        ptr = strtok(NULL, ",\n");
    }
    if (size != (ix + 1))
    {
        quitError("Input file is corrupt: kernel size does not match\n", ERROR_FORMAT);
    }
    *rev = m_rev;
    *rev_ct = ct;
}

void createReactionMapping(int linelength, char *line, int rx_count, int
        ***rx_map, int *map_size, int **map_lengths)
{
    if (linelength < 1)
    {
        quitError("Wrong format of input file: linear dependent reactions are missing\n", 
                ERROR_FORMAT);
    }
    int i;
    int *in_array = calloc(rx_count, sizeof(int*));
    int ix = -1;
    char* ptr = strtok(line, ",\n");
    while(ptr != NULL)
    {
        ix++;
        int val = (int) strtol(ptr, NULL, 10);
        in_array[ix] = val;
        ptr = strtok(NULL, ",\n");
    }
    if (rx_count != (ix + 1))
    {
        quitError("Input file is corrupt: reaction count does not match\n", ERROR_FORMAT);
    }
    int **map = calloc(rx_count, sizeof(int**));
    int *maps = calloc(rx_count, sizeof(int*));
    ix = 0;
    for (i = 0; i < rx_count; i++)
    {
        if (in_array[i] < 0)
        {
            map[ix] = calloc(sizeof(int*), rx_count);
            map[ix][0] = i;
            int counter = 1;
            int j;
            for (j = 0; j < rx_count; j++)
            {
                if (in_array[j] == i)
                {
                    map[ix][counter] = j;
                    counter++;
                }
            }
            map[ix] = (int*) realloc(map[ix], sizeof(int) * counter);
            maps[ix] = counter;
            ix++;
        }
    }
    map = (int**) realloc(map, sizeof(int*) * ix);
    maps = (int*) realloc(maps, sizeof(int) * ix);
    *rx_map = map;
    *map_size = ix;
    *map_lengths = maps;
    free(in_array);
}


void printReactionNames(int linelength, FILE *f, char *line, int rx_count)
{
    if (linelength < 1)
    {
        quitError("Wrong format of input file: reaction names are missing\n", 
                ERROR_FORMAT);
    }
    char* copy = malloc(strlen(line) + 1);
    strcpy(copy, line);
    int ix = 0;
    char* ptr = strtok(line, ",\n");
    while(ptr != NULL)
    {
        ix++;
        ptr = strtok(NULL, ",\n");
    }
    if (ix != rx_count)
    {
        quitError("Wrong format of input file: wrong number of reaction names\n", 
                ERROR_FORMAT);
    }
    fprintf(f, "%s", copy);
    free(copy);
}

void printVector(FILE *f, char *v, int len, int **map, int *maps, int
        all_rx_len, int use_sep, char* sep, int bitsize, char* fwd_str, char*
        rev_str)
{
    // char *exp = calloc(1, bitsize);
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
        if (use_sep && i > 0)
        {
            fprintf(f, "%s", sep);
        }
        if (BITTEST(exp, i))
        {
            fprintf(f, "%s", fwd_str);
        }
        else
        {
            fprintf(f, "%s", rev_str);
        }
    }
    fprintf(f, "\n");
    free(exp);
}

// return values
// 0: not print
// 1: print
// 2: end
int isPrintableFluxtope(int* act_step, int* step_ct, unsigned long long*
        act_size, unsigned long long* tot_size, int start, int end, int log) 
{
    if (*step_ct > *act_step)
    {
        if (log > 1)
        {
            if (*act_step < start)
            {
                fprintf(stdout, "Found ");
            }
            else
            {
                fprintf(stdout, "Printed ");
            }
            fprintf(stdout, "%llu fluxtopes for step %d\n", *act_size, *act_step);
        }
        *act_step = *step_ct;
        *act_size = 1;
    }
    else
    {
        *act_size = *act_size + 1;
    }
    if (*act_step >= start && *act_step <= end)
    {
        *tot_size = *tot_size + 1;
        return 1;
    }
    if (*step_ct > end)
    {
        return 2;
    }
    return 0;
}

void parseFirstLine(int linelength, char* line, int* reactions, int* kernel_len)
{
    if (linelength > 20)
    {
        if (strncmp(line, "compressedfluxtopes:", 20) == 0)
        {
            char* ptr = strtok(line, ":\n");
            if (ptr != NULL)
            {
                ptr = strtok(NULL, ":\n");
                if (ptr != NULL)
                {
                    *reactions = (int) strtol(ptr, NULL, 10);
                    ptr = strtok(NULL, ":\n");
                    if (ptr != NULL)
                    {
                        *kernel_len = (int) strtol(ptr, NULL, 10);
                    }
                }
            }
        }
    }
    if (*reactions == 0 || *kernel_len == 0)
    {
        quitError("Wrong format of input file\n", ERROR_FORMAT);
    }
}
