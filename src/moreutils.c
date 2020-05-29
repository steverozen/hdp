#include <stdlib.h>
#include <R.h>

/*
 * More utilities added in package hdpx.
 */

/* Not sure which of these will be best.... */

void* malloc_and_check(size_t x) { /* We raise an error  */
  void *r = malloc(x);
  if (NULL == r) Rf_error("malloc failed");
  return r;
}

void* mallocR(size_t x) {
  /* R checks and raises the error, and frees the memory at the end of
    .C, .Call, or .External. But an attempt to use this caused problems
    because it is incompatible with normal free (not to be confused with
    Free == Rf_Free)*/
  return(void *) R_alloc(x, 1);
}


/* Another option is Rf_Calloc, Rf_Free, and set an on.exit to Rf_free; see
 * 6.1.2 in
 * https://cran.r-project.org/doc/manuals/R-exts.html#Memory-allocation
 */

