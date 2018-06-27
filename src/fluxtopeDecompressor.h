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
#include <limits.h>
#include <time.h>

#define BITSIZE CHAR_BIT
#define BITMASK(b) (1 << ((b) % BITSIZE))
#define BITSLOT(b) ((b) / BITSIZE)
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a, b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a, b) ((a)[BITSLOT(b)] & BITMASK(b))

#define MAX_ARGS     8
#define ARG_FT_IN    0
#define ARG_FT_OUT   1
#define ARG_FT_START 2
#define ARG_FT_END   3
#define ARG_SEP      4
#define ARG_FWD      5
#define ARG_REV      6
#define ARG_LOG      7

#define ERROR_ARGS   1
#define ERROR_FILE   2
#define ERROR_FORMAT 3
#define ERROR_RAM    4

unsigned long getBitsize(unsigned long count);
void open_file(FILE **file, char *name, char *handle);
char* getTime();
void handle_arguments(int argc, char *argv[], char** ft_in, char** ft_out,
        int *ft_start, int *ft_end, int* use_sep, char** sep, char** fwd_str,
        char** rev_str, int* log);
char* getArg (int argc, char *argv[], char *opt);
void readArgs (int argc, char *argv[], int optc, char *optv[], char *optr[]);
void usage (char* description, char* usage, int optc, char* optv[],
        char *optd[]);
void quitError (char *message, int rv);
void quitErrorArg (char *message, char *arg, int rv);

void printVector(FILE *f, char *v, int len, int **map, int *maps, int
        all_rx_len, int use_sep, char* sep, int bitsize, char* fwd_str, char*
        rev_str);
void printReactionNames(int linelength, FILE *f, char *line, int rx_count);
void decompressFluxtopes(char* in, char* out, int use_sep, char* sep, int
        start, int end, char* fwd_str, char* rev_str, int log);

void parseFirstLine(int linelength, char* line, int* reactions, int*
        kernel_len);
void parseReversibles(int linelength, char* line, int kernel_len, char** rev,
        int* rev_ct);
void createReactionMapping(int linelength, char *line, int rx_count, int
        ***rx_map, int *map_size, int **map_lengths);
int isPrintableFluxtope(int* act_step, int* step_ct, unsigned long long*
        act_size, unsigned long long* tot_size, int start, int end, int log);
