#define R_NO_REMAP
#include <Rinternals.h>
#include <string.h>

/* Finite times and a stable R order by (patient, time) are supplied by
   the R validator. Return duplicate flags in the original record order.
   Only the output vector is allocated during the scan. */
SEXP gpc_mark_measurement_repeats(SEXP person, SEXP time, SEXP order) {
    R_xlen_t n = XLENGTH(order);
    if (TYPEOF(person) != INTSXP || TYPEOF(time) != REALSXP ||
        TYPEOF(order) != INTSXP || XLENGTH(person) != n || XLENGTH(time) != n)
        Rf_error("Invalid measurement key inputs.");
    SEXP result = PROTECT(Rf_allocVector(LGLSXP, n));
    int *flags = LOGICAL(result);
    memset(flags, 0, n * sizeof(*flags));
    const int *ids = INTEGER(person), *rows = INTEGER(order);
    const double *times = REAL(time);
    int previous = 0;
    for (R_xlen_t i = 0; i < n; ++i) {
        if (rows[i] < 1 || rows[i] > n) Rf_error("Invalid measurement row index.");
        int current = rows[i] - 1;
        if (i && ids[current] == ids[previous] && times[current] == times[previous])
            flags[current] = 1;
        previous = current;
    }
    UNPROTECT(1);
    return result;
}
