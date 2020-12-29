CC = icc
#CFLAGS= -pedantic -Wall -Weffc++ -Wextra -I/usr/include/openblas -lopenblas

succprob: Essential_new.cpp Read.cpp Parameters.cpp Product.cpp QA.cpp Lanczos_succ.cpp
	$(CC) -qopenmp-O3 Essential_new.cpp Read.cpp Parameters.cpp QA.cpp Lanczos_succ.cpp Product.cpp 
-o succ_prob -mkl

instover: Essential.cpp Read.cpp Parameters_Nonstoq.cpp Product_Nonstoq.cpp StateOverlap.cpp QA_StateOverlap.cpp 
	$(CC) -O3 Essential.cpp Read.cpp Parameters_Nonstoq.cpp Product_Nonstoq.cpp StateOverlap.cpp QA_StateOverlap.cpp -o lanczos_instoverlap -mkl
clean: 
	rm -f *.o lanczos_prob lanczos_instoverlap 
