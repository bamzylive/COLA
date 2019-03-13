simple:
	gcc -c -g *.c && gcc -o test *.o -lcrypto

clean:
	rm -f *.o test	

          
