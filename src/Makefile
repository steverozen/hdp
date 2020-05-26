
objects = R-base.o R-conparam.o R-dp.o R-hdp.o R-hdpMultinomial_iterate.o R-multinomial.o R-utils.o randutils.o

hdpx.so : $(objects)
	gcc -std=gnu99 -shared -L/opt/R/R-4.0.0/lib64/R/lib -L/usr/local/lib64 -o hdpx.so $(objects) -L/opt/R/R-4.0.0/lib64/R/lib -lR


$(objects): %.o: %.c
	gcc -std=gnu99 -I"/opt/R/R-4.0.0/lib64/R/include" -DNDEBUG   -I/usr/local/include   -fpic  -g -O0  -c $< -o $@

clean:
	rm -f hdpx.so *.o

# hdpx.so: R-base.o R-conparam.o R-dp.o R-hdp.o R-hdpMultinomial_iterate.o R-multinomial.o R-utils.o randutils.o
#	gcc -std=gnu99 -shared -L/opt/R/R-4.0.0/lib64/R/lib -L/usr/local/lib64 -o hdpx.so R-base.o R-conparam.o R-dp.o R-hdp.o R-hdpMultinomial_iterate.o R-multinomial.o R-utils.o randutils.o -L/opt/R/R-4.0.0/lib64/R/lib -lR

