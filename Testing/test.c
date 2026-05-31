#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include "LibPrimeNum.h"

int main() {
    printf("Тестирование библиотеки простых чисел\n");
    // Обязательная инициализация глобального контекста таблицы простых чисел библиотеки
    prime_lib_init();

    FILE *f = fopen("test.txt", "r");
    if (f == NULL) {
        printf("Ошибка открытия файла\n");
        return 1;
    }

    FILE *f_prost = fopen("testeas.txt", "r");
    if (f_prost == NULL) {
        printf("Ошибка открытия файла\n");
        return 1;
    }

    // Подсчет количества чисел в файле
    int n = 0;
    int n_prost = 0;
    mpz_t x;
    mpz_init(x);
    while (gmp_fscanf(f, "%Zd", x) == 1) n++;
    while (gmp_fscanf(f_prost, "%Zd", x) == 1) n_prost++;
    mpz_clear(x);
    printf("\nКоличество чисел в файле: %d\n", n);
    fseek(f, 0, SEEK_SET);
    fseek(f_prost, 0, SEEK_SET);


    // Чтение чисел из файла в массив
    mpz_t *X = malloc(n * sizeof(mpz_t));
    for (int i=0; i<n; i++) mpz_init(X[i]);

    mpz_t *X_prost = malloc(n_prost * sizeof(mpz_t));
    for (int i=0; i<n_prost; i++) mpz_init(X_prost[i]);

    int i = 0;
    while (gmp_fscanf(f, "%Zd", X[i]) == 1) i++;
    printf("В массив X успешно записано %d чисел\n", n);
    fclose(f);

    i = 0;
    while (gmp_fscanf(f_prost, "%Zd", X_prost[i]) == 1) i++;
    printf("В массив X_prost успешно записано %d чисел\n", n_prost);
    fclose(f_prost);

    printf("Тестирование скорости и точности вероятностных алгоритмов на больших числах (от 617 до 1234 знаков):\n");

    printf("\n1. Алгоритмом Миллера-Рабина:\n");
    for (int i=0; i<n; i++) {

        clock_t start = clock();
        char turn = test_miller_rabin(X[i]);
        clock_t end = clock();
        double time = (double)(end - start) / CLOCKS_PER_SEC;
        if (turn==0) printf("Число №%hd из файла составное\n", i+1);
        else if (turn==1) printf("Число №%hd из файла вероятно простое\n", i+1);
        else if (turn==2) printf("Число №%hd из файла точно простое\n", i+1);

        printf("Время рассчета для числа %hd: %f секунд\n", i+1, time);
    }

    printf("\n2. Алгоритмом Ферма:\n");
    for (int i=0; i<n; i++) {

        clock_t start = clock();
        char turn = test_fermat(X[i], 10);
        clock_t end = clock();
        double time = (double)(end - start) / CLOCKS_PER_SEC;
        if (turn==0) printf("Число №%hd из файла составное\n", i+1);
        else if (turn==1) printf("Число №%hd из файла вероятно простое\n", i+1);
        else if (turn==2) printf("Число №%hd из файла точно простое\n", i+1);

        printf("Время рассчета для числа %hd: %f секунд\n", i+1, time);
    }



    printf("\n3. Алгоритмом Соловея-Штрассена:\n");
    for (int i=0; i<n; i++) {

        clock_t start = clock();
        char turn = test_solovay_strassen(X[i], 10);
        clock_t end = clock();
        double time = (double)(end - start) / CLOCKS_PER_SEC;
        if (turn==0) printf("Число №%hd из файла составное\n", i+1);
        else if (turn==1) printf("Число №%hd из файла вероятно простое\n", i+1);
        else if (turn==2) printf("Число №%hd из файла точно простое\n", i+1);

        printf("Время рассчета для числа %hd: %f секунд\n", i+1, time);
    }

    printf("\nТестирование скорости и точности вероятностных и детерминированных алгоритмов на малых числах (от X до Y знаков):\n");

    printf("\n1. Алгоритмом Миллера-Рабина:\n");
    for (int i=0; i<n_prost; i++) {

        clock_t start = clock();
        char turn = test_miller_rabin(X_prost[i]);
        clock_t end = clock();
        double time = (double)(end - start) / CLOCKS_PER_SEC;
        if (turn==0) printf("Число №%hd из файла составное\n", i+1);
        else if (turn==1) printf("Число №%hd из файла вероятно простое\n", i+1);
        else if (turn==2) printf("Число №%hd из файла точно простое\n", i+1);

        printf("Время рассчета для числа %hd: %f секунд\n", i+1, time);
    }

    printf("\n2. Алгоритмом Ферма:\n");
    for (int i=0; i<n_prost; i++) {

        clock_t start = clock();
        char turn = test_fermat(X_prost[i], 10);
        clock_t end = clock();
        double time = (double)(end - start) / CLOCKS_PER_SEC;
        if (turn==0) printf("Число №%hd из файла составное\n", i+1);
        else if (turn==1) printf("Число №%hd из файла вероятно простое\n", i+1);
        else if (turn==2) printf("Число №%hd из файла точно простое\n", i+1);

        printf("Время рассчета для числа %hd: %f секунд\n", i+1, time);
    }

    printf("\n3. Алгоритмом Соловея-Штрассена:\n");
    for (int i=0; i<n_prost; i++) {

        clock_t start = clock();
        char turn = test_solovay_strassen(X_prost[i], 10);
        clock_t end = clock();
        double time = (double)(end - start) / CLOCKS_PER_SEC;
        if (turn==0) printf("Число №%hd из файла составное\n", i+1);
        else if (turn==1) printf("Число №%hd из файла вероятно простое\n", i+1);
        else if (turn==2) printf("Число №%hd из файла точно простое\n", i+1);

        printf("Время рассчета для числа %hd: %f секунд\n", i+1, time);
    }

//    printf("\n4. Алгоритмом Агравал-Каян-Саксена:\n");
//    for (int i=0; i<n; i++) {

//        clock_t start = clock();
//        char turn = test_aks(X_prost[i]);
//        clock_t end = clock();
//        double time = (double)(end - start) / CLOCKS_PER_SEC;
//        if (turn==0) printf("Число №%hd из файла составное\n", i+1);
//        else if (turn==1) printf("Число №%hd из файла вероятно простое\n", i+1);
//        else if (turn==2) printf("Число №%hd из файла точно простое\n", i+1);

//        printf("Время рассчета для числа %hd: %f секунд\n", i+1, time);
//    }

    printf("\n5. Алгоритмом перебора делителей:\n");
    for (int i=0; i<n_prost; i++) {

        clock_t start = clock();
        char turn = test_trivial_division(X_prost[i]);
        clock_t end = clock();
        double time = (double)(end - start) / CLOCKS_PER_SEC;
        if (turn==0) printf("Число №%hd из файла составное\n", i+1);
        else if (turn==1) printf("Число №%hd из файла вероятно простое\n", i+1);
        else if (turn==2) printf("Число №%hd из файла точно простое\n", i+1);

        printf("Время рассчета для числа %hd: %f секунд\n", i+1, time);
    }

    printf("\nГенераторы простых чисел: \n");

    printf("1. Решето Эратосфена (1000000 чисел): \n");
    FILE *E = fopen("Eratosthenes.txt", "w");
    if (E == NULL) {
        printf("Ошибка открытия файла\n");
        return 1;
    }
    unsigned long *Eratosthenes = malloc(1000000 * sizeof(unsigned long));
    generate_primes_eratosthenes(1000000, Eratosthenes);
    for (int i=0; i<1000000; i++) {
        fprintf(E, "%lu\n", Eratosthenes[i]);
    }
    free(Eratosthenes);
    fclose(E);


    printf("2. Решето Аткина (1000000 чисел): \n");
    FILE *A = fopen("Atkin.txt", "w");
    if (A == NULL) {
        printf("Ошибка открытия файла\n");
        return 1;
    }
    unsigned long *Atkin = malloc(1000000 * sizeof(unsigned long));
    generate_primes_atkin(1000000, Atkin);
    for (int i=0; i<1000000; i++) {
        fprintf(A, "%lu\n", Atkin[i]);
    }
    free(Atkin);
    fclose(A);

    // Освобождение памяти контекста
    prime_lib_cleanup();
    free(X); free(X_prost);
    printf("Тестирование завершено успешно\n");
    return 0;
}
