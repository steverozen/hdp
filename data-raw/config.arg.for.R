./configure \
    --prefix=/opt/R/${R_VERSION} \
    --enable-memory-profiling \
    --enable-R-shlib \
    --with-blas \
    --with-lapack \
    --with-system-valgrind-headers \
    --with-valgrind-instrumentation \
    --with-pcre1

    export CFLAGS="-g -O1"
# export MAIN_FFLAGS="-g -02"
