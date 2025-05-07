gcc -c -o cf.o src/CustardFlow.c
gcc -c -o main.o src/Main.c
gcc -o build/mm cf.o main.o