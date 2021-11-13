rcp:
	gcc -Wall rcp.c brute.c -o rcp -lm
so:
	gcc -fPIC -c -Wall rcp.c
	gcc -fPIC -c -Wall brute.c
	ld -shared rcp.o brute.o -o librcp.so
