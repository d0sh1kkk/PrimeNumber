#ifndef PRIME_LIB_H
#define PRIME_LIB_H

#include <gmp.h>
#include <stdbool.h>

// Статусы простоты
typedef enum {
    PRIME_COMPOSITE = 0, // Составное
    PRIME_PROBABLY = 1, // Вероятно простое
    PRIME_DEFINITELY = 2 // Определенно простое
} PrimeStatus;

// Инициализация и очистка (для таблицы малых простых чисел)
void prime_lib_init();
void prime_lib_cleanup();

// Детерминированные тесты
PrimeStatus test_trivial_division(mpz_t n);
PrimeStatus test_aks(mpz_t n);
unsigned long test_mersenne_lucas_lehmer(unsigned long p, mpz_t *out_primes);

// Вероятностные тесты
PrimeStatus test_miller_rabin(mpz_t n);
PrimeStatus test_fermat(mpz_t n, int iterations);
PrimeStatus test_solovay_strassen(mpz_t n, int iterations);

// Генераторы (Решета)
void generate_primes_eratosthenes(int n, unsigned long *out_primes);
void generate_primes_atkin(int n, unsigned long *out_count);

#endif // PRIME_LIB_H
