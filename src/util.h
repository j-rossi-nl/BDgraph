#ifndef UTIL_H
#define UTIL_H

#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T // For Fortran character strings
#endif

#ifndef FCONE
#define FCONE
#endif

// Helper: for a 2-D matrix stored as flat array in COLUMN-MAJOR fashion
// Retrives the index in the flat array of the element at (row, col)
inline int index_in_array(int row, int col, int nrows)
{
    return col * nrows + row;
}

#endif
