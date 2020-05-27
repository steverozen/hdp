#include <stdlib.h>
#include <R.h>

/*
 * More utilities added in hdpx.
 */

/* Not sure which of these will be best.... */

static void* x_safe_malloc(size_t x) { /* We raise there error  */
  void *r = malloc(x);
  if (NULL == r) Rf_error("malloc failed");
  return r;
}

static void* R_safe_malloc(size_t x) {
  /* R checks and raises the error, and frees the memory at the end of
    .C, .Call, or .External */
  return(void *) R_alloc(x, 1);
}

void *hdp_malloc(size_t x) {
  return R_safe_malloc(x);
}

/* Another option is Rf_Calloc, Rf_Free, and set an on.exit to free; see
 * 6.1.2 in
 * https://cran.r-project.org/doc/manuals/R-exts.html#Memory-allocation
 */

