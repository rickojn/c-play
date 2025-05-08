gcc -c -o build/ma.o src/mul-add.s
gcc -c -o build/main.o src/Main.c
gcc -o build/mamu  build/main.o build/ma.o