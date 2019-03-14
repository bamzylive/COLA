simple:
	gcc -c -g *.c -I/usr/local/Cellar/openssl/1.0.2q/include && gcc -o test *.o -L/usr/local/Cellar/openssl/1.0.2q/lib -lcrypto

clean:
	rm -f *.o test	

          
