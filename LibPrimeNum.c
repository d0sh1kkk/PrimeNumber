#include "LibPrimeNum.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <omp.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#define NUM_PROS_TABL 512
unsigned long *prostie;

// РЕШЕТО ЭРАТОСФЕНА
void generate_primes_eratosthenes(int n, unsigned long* temp_prostie) {
    if (n <= 0 || temp_prostie == NULL) return; // Если пользователь ввел некорректное значение

    temp_prostie[0] = 2;
    if (n == 1) return;

    // Определение границ
    long long limit;
    if (n < 6) limit = 15; // Фиксированная граница для малых n
    else limit = (long long)(n * (log(n) + log(log(n)))); // верхняя оценка для больших n

    // Выделение памяти
    bool *is_prime = malloc((limit + 1) * sizeof(bool));
    if (!is_prime) {
        printf("Ошибка памяти\n");
        return;
    }


    memset(is_prime, true, (limit + 1) * sizeof(bool)); // Сначала все числа считаем простыми
    is_prime[0] = is_prime[1] = false; // 0 и 1 не являются простыми

    // Решето Эратосфена
    for (long long p = 3; p * p <= limit; p += 2) {
        if (is_prime[p]) {
            // Отмечаем все кратные p как составные
            for (long long i = p * p; i <= limit; i += p)
                is_prime[i] = false;
        }
    }

    int count = 1;
    for (long long p = 3; count < n && p <= limit; p += 2) {
        if (is_prime[p]) {
            temp_prostie[count] = (long int)p; // Записываем простое число в файл
            count++;
        }
    }

    free(is_prime); // Освобождаем память
}

// РЕШЕТО АТКИНА
void generate_primes_atkin(int n, unsigned long* temp_prostie) {
    if (n <= 0) return;

    // Определение границ
    long long limit;
    if (n < 6) limit = 15; // Фиксированная граница для малых n
    else limit = (long long)(n * (log(n) + log(log(n)))); // Верхняя оценка для больших n

    // Выделение памяти
    bool *is_prime = calloc(limit + 1, sizeof(bool));
    if (!is_prime) {
        return;
    }

    // Предварительно обрабатываем 2, 3 и 5 (база для модуля 60)
    if (limit >= 2) is_prime[2] = true;
    if (limit >= 3) is_prime[3] = true;
    if (limit >= 5) is_prime[5] = true;

    // Предвычисляем квадраты, чтобы не множить в циклах
    for (long long x = 1; x * x <= limit; x++) {
        long long x2 = x * x;
        long long fourX2 = 4 * x2;
        long long threeX2 = 3 * x2;

        // Основной алгоритм: перебор x и y для квадратичных форм
        for (long long y = 1; y * y <= limit; y++) { // Перебор всех возможных пар (x, y)
            long long y2 = y * y;
            // Формула 1: 4x^2 + y^2 = m
            long long m = fourX2 + y2;
            if (m <= limit && (m % 12 == 1 || m % 12 == 5)) {
                is_prime[m] ^= true; // Инвертируем состояние
            }

            // Формула 2: 3x^2 + y^2 = m
            m = threeX2 + y2;
            if (m <= limit && m % 12 == 7) {
                is_prime[m] ^= true;
            }

            // Формула 3: 3x^2 - y^2 = m
            m = threeX2 - y2;
            if (x > y && m <= limit && m % 12 == 11) {
                is_prime[m] ^= true;
            }
        }
    }

    // Отсеиваем квадраты простых чисел
    for (long long p = 5; p * p <= limit; p++) {
        if (is_prime[p]) {
            long long p2 = p * p;
            for (long long i = p2; i <= limit; i += p2)
                is_prime[i] = false; // Отсекаем числа, кратные квадрату p
        }
    }

    // Запись результатов
    int count = 0;
    for (long long p = 2; count < n && p <= limit; p++) {
        if (is_prime[p]) {
            temp_prostie[count] = (unsigned long)p;
            count++;
        }
    }

    free(is_prime); // Освобождаем память
}

// ИНИЦИАЛИЗАЦИЯ
void prime_lib_init () {
    prostie = malloc(NUM_PROS_TABL * sizeof(long int));
    if (prostie == NULL) {
        return;
    }
    generate_primes_eratosthenes(NUM_PROS_TABL, prostie);
}

// ТЕСТ ФЕРМА
PrimeStatus test_fermat(mpz_t n, int iterations) {
    // Обработка базовых случаев
    if (mpz_cmp_ui(n, 1) <= 0) return 0; // 0 и 1 не простые
    if (mpz_cmp_ui(n, 3) <= 0) return 2; // 2 и 3 простые
    if (mpz_even_p(n)) return 0; // Проверка на четность

    // Проверка по таблице малых простых чисел
    for (int i = 0; i < NUM_PROS_TABL; i++) {
        if (mpz_cmp_ui(n, prostie[i]) == 0) return 2;
        if (mpz_divisible_ui_p(n, prostie[i])) return 0;
    }

    int is_probably_prime = 1;
    // Параллельный блок
    #pragma omp parallel
    {
        // Инициализация переменных для каждого потока
        mpz_t n_minus_1, res;
        mpz_init(n_minus_1);
        mpz_init(res);

        mpz_sub_ui(n_minus_1, n, 1);

        // Распределение итераций между потоками
        #pragma omp for schedule(static)
        for (int i = 0; i < iterations; i++) {
            if (!is_probably_prime) continue; // Если уже составное, пропускаем
            mpz_set_ui (res, prostie[i]);
            mpz_powm(res, res, n_minus_1, n); // Вычисляем res = a^(n-1) mod n

            // Если res != 1, число точно составное
            if (mpz_cmp_ui(res, 1) != 0) {
                is_probably_prime = 0;
            }
        }
        mpz_clears(n_minus_1, res, NULL);
    }
    return is_probably_prime;
}


// АЛГОРИТМ МИЛЛЕРА-РАБИНА
// Оптимальный рассчет количества итераций для теста Миллера-Рабина
int optimal_iteration(mpz_t n) {
    int k = mpz_sizeinbase(n, 2); // Количество бит в n

    // Если число условно мало
    if (k < prostie[NUM_PROS_TABL - 1]) return 1;

    // Если число условно большое
    double target_log = -80.0; // log10 от 10^-80
    double log_p = log10(4.0 * k) - sqrt((double)k) * log10(2.0);

    if (log_p > -0.602) log_p = -0.602; // Ограничиваем сверху log10(1/4)

    int iteration = (int)ceil(target_log / log_p);

    // printf("\nИтераций предв: %d\n", iteration);

    if (iteration > 40) iteration = 40;
    if (iteration < 2) iteration = 2;

    return iteration;
}

// Функция теста Миллера-Рабина
PrimeStatus test_miller_rabin(mpz_t n) {
    // Базовые проверки
    if (mpz_cmp_ui(n, 2) < 0) return 0;
    if (mpz_cmp_ui(n, 2) == 0) return 2;
    if (mpz_even_p(n)) return 0;

    // Проверка по таблице малых простых чисел
    for (int i = 0; i < NUM_PROS_TABL; i++) {
        if (mpz_cmp_ui(n, prostie[i]) == 0) return 2;
        if (mpz_divisible_ui_p(n, prostie[i])) return 0;
    }

    // Подготовка переменных: n - 1 = 2^s * d
    mpz_t d, n_minus_1, x;
    mpz_inits(d, n_minus_1, x, NULL);
    mpz_sub_ui(n_minus_1, n, 1); // n - 1

    long int s = mpz_scan1(n_minus_1, 0); // Поиск первого бита = 1
    mpz_tdiv_q_2exp(d, n_minus_1, s); // d = (n - 1) / 2^s

    char is_prime = 1;
    int iteration = optimal_iteration(n);
    // Основной цикл Миллера-Рабина
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < iteration; i++) {
        if (mpz_cmp_ui(n, prostie[i]) <= 0) continue; // Если a >= n, выходим
        mpz_set_ui(x, prostie[i]); // x = a
        mpz_powm(x, x, d, n); // x = a^d mod n
        // Если x = 1 или n-1, то число проходит тест для a
        if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, n_minus_1) == 0) continue;

        // Если первый этап не пройден, проверяем второй
        char composite_flag = 1; // Флаг для составного числа
        // Повторяем возведение в квадрат s-1 раз
        for (unsigned long int r = 1; r < s; r++) {
            mpz_powm_ui(x, x, 2, n); // x = x^2 mod n
            // Если x = n-1, то число проходит тест для a
            if (mpz_cmp(x, n_minus_1) == 0) {
                composite_flag = 0;
                break;
            }
        }
        // Если composite_flag все еще 1, то число составное
        if (composite_flag) {
            is_prime = 0;
            continue;
        }
    }

    mpz_clears(d, n_minus_1, x, NULL); // Освобождаем память
    return is_prime;
}

// АЛГОРИТМ AKS
// Структура для представления многочлена Z_n[x]/(x^r - 1)
typedef struct {
    mpz_t *coeffs;
    unsigned long r;
} Poly;

// Инициализация многочлена (массив mpz_t)
void poly_init(Poly *p, unsigned long r) {
    p->r = r;
    p->coeffs = malloc(r * sizeof(mpz_t));
    for (unsigned long i = 0; i < r; i++) {
        mpz_init_set_ui(p->coeffs[i], 0);
    }
}

// Освобождение памяти многочлена (mpz_t)
void poly_clear(Poly *p) {
    for (unsigned long i = 0; i < p->r; i++) {
        mpz_clear(p->coeffs[i]);
    }
    free(p->coeffs);
}

// Перемножение двух многочленов
void poly_mul(Poly *res, const Poly *p1, const Poly *p2, mpz_t n, Poly *temp) {
    unsigned long r = p1->r;

    for (unsigned long i = 0; i < r; i++) mpz_set_ui(temp->coeffs[i], 0);

    for (unsigned long i = 0; i < r; i++) {
        if (mpz_sgn(p1->coeffs[i]) == 0) continue;
        for (unsigned long j = 0; j < r; j++) {
            if (mpz_sgn(p2->coeffs[j]) == 0) continue;

            unsigned long k = (i + j) % r;
            mpz_addmul(temp->coeffs[k], p1->coeffs[i], p2->coeffs[j]);

            // Если коэффициент слишком большой, берем остаток по модулю n
            if (mpz_sizeinbase(temp->coeffs[k],2) > mpz_sizeinbase(n,2) + 64) {
                mpz_mod(temp->coeffs[k], temp->coeffs[k], n);
            }
        }
    }

    for (unsigned long i = 0; i < r; i++) {
        mpz_mod(res->coeffs[i], temp->coeffs[i], n);
    }
}

// Специальная версия возведения в квадрат для полиномов
void poly_sqr(Poly *res, const Poly *p, mpz_t n, Poly *temp) {
    unsigned long r = p->r;
    for (unsigned long i = 0; i < r; i++) mpz_set_ui(temp->coeffs[i], 0);

    for (unsigned long i = 0; i < r; i++) {
        if (mpz_sgn(p->coeffs[i]) == 0) continue;

        // Квадрат текущего члена: a_i^2 * x^(2i)
        unsigned long k2 = (i + i) % r;
        mpz_addmul(temp->coeffs[k2], p->coeffs[i], p->coeffs[i]);

        for (unsigned long j = i + 1; j < r; j++) {
            if (mpz_sgn(p->coeffs[j]) == 0) continue;

            // Удвоенное произведение: 2 * a_i * a_j * x^(i+j)
            unsigned long k = (i + j) % r;
            mpz_t prod; mpz_init(prod);
            mpz_mul(prod, p->coeffs[i], p->coeffs[j]);
            mpz_addmul_ui(temp->coeffs[k], prod, 2);
            mpz_clear(prod);
        }
        mpz_mod(temp->coeffs[i], temp->coeffs[i], n); // Частичный модуль
    }
    for (unsigned long i = 0; i < r; i++) mpz_mod(res->coeffs[i], temp->coeffs[i], n);
}


void poly_pow(Poly *res, const Poly *base, mpz_t exp, mpz_t n, Poly *t_res, Poly *t_base, Poly *backup) {
    unsigned long r = base->r;
    for (unsigned long i = 0; i < r; i++) {
        mpz_set_ui(t_res->coeffs[i], 0);
        mpz_set(t_base->coeffs[i], base->coeffs[i]);
    }
    mpz_set_ui(t_res->coeffs[0], 1);

    size_t num_bits = mpz_sizeinbase(exp, 2);
    for (size_t bit = 0; bit < num_bits; bit++) {
        if (mpz_tstbit(exp, bit)) {
            poly_mul(t_res, t_res, t_base, n, backup);
        }
        if (bit < num_bits - 1) {
            poly_sqr(t_base, t_base, n, backup); // Используем SQR вместо MUL
        }
    }
    for (unsigned long i = 0; i < r; i++) mpz_set(res->coeffs[i], t_res->coeffs[i]);
}

PrimeStatus test_aks(mpz_t n) {
    // Базовые проверки
    if (mpz_cmp_ui(n, 2) < 0) return 0;
    if (mpz_cmp_ui(n, 2) == 0) return 2;
    if (mpz_even_p(n)) return 0;

    for (int i = 0; i < NUM_PROS_TABL; i++) {
        if (mpz_cmp_ui(n, prostie[i]) == 0) return 2;
        if (mpz_divisible_ui_p(n, prostie[i])) return 0;
    }

    if (mpz_perfect_power_p(n)) return 0;

    unsigned long log2_n = mpz_sizeinbase(n, 2);
    unsigned long limit = log2_n * log2_n;
    unsigned long r = 2;
    mpz_t mod_res, mpz_r;
    mpz_inits(mod_res, mpz_r, NULL);

    while (1) {
        mpz_set_ui(mpz_r, r);
        mpz_gcd(mod_res, n, mpz_r);
        if (mpz_cmp_ui(mod_res, 1) > 0 && mpz_cmp(mod_res, n) < 0) {
            mpz_clears(mod_res, mpz_r, NULL);
            return 0;
        }
        char order_ok = 1;
        for (unsigned long k = 1; k <= limit; k++) {
            mpz_powm_ui(mod_res, n, k, mpz_r);
            if (mpz_cmp_ui(mod_res, 1) == 0) { order_ok = 0; break; }
        }
        if (order_ok) break;
        r++;
    }
    mpz_clears(mod_res, mpz_r, NULL);

    if (mpz_cmp_ui(n, r) <= 0) return 2;

    unsigned long a_limit = (unsigned long)(sqrt(r) * log2_n);

    // ВЫДЕЛЯЕМ ПАМЯТЬ ОДИН РАЗ НА ВЕСЬ ТЕСТ ЧИСЛА
    Poly base_poly, res_poly, t_res, t_base, backup;
    poly_init(&base_poly, r);
    poly_init(&res_poly, r);
    poly_init(&t_res, r);
    poly_init(&t_base, r);
    poly_init(&backup, r);

    char is_prime = 1;
    for (unsigned long a = 1; a <= a_limit; a++) {
        for (unsigned long i = 0; i < r; i++) mpz_set_ui(base_poly.coeffs[i], 0);
        mpz_set_ui(base_poly.coeffs[0], a);
        mpz_set_ui(base_poly.coeffs[1], 1);

        poly_pow(&res_poly, &base_poly, n, n, &t_res, &t_base, &backup);

        unsigned long n_mod_r = mpz_fdiv_ui(n, r);
        for (unsigned long i = 0; i < r; i++) {
            if (i == 0) {
                if (mpz_cmp_ui(res_poly.coeffs[i], a) != 0) { is_prime = 0; break; }
            } else if (i == n_mod_r) {
                if (mpz_cmp_ui(res_poly.coeffs[i], 1) != 0) { is_prime = 0; break; }
            } else {
                if (mpz_sgn(res_poly.coeffs[i]) != 0) { is_prime = 0; break; }
            }
        }
        if (!is_prime) break;
    }

    poly_clear(&base_poly); poly_clear(&res_poly);
    poly_clear(&t_res); poly_clear(&t_base); poly_clear(&backup);

    return is_prime ? 2 : 0;
}

// АЛГОРИТМ ЛЮКА-ЛЕМЕРА
// Упрощенное решето Эратосфена для поиска первых n простых чисел
bool* first_primes_ML(int n) {
    if (n <= 0) return NULL; // Если пользователь ввел некорректное значение
    // Выделение памяти
    bool* is_prime = malloc((n + 1) * sizeof(bool));
    if (!is_prime) {
        printf("Ошибка памяти\n");
        return NULL;
    }

    memset(is_prime, true, (n + 1) * sizeof(bool)); // Сначала все числа считаем простыми
    is_prime[0] = is_prime[1] = false; // 0 и 1 не являются простыми

    // Решето Эратосфена
    for (long long p = 3; p * p <= n; p += 2) {
        if (is_prime[p]) {
            // Отмечаем все кратные p как составные
            for (long long i = p * p; i <= n; i += p)
                is_prime[i] = false;
        }
    }
    return is_prime;
}

// Оптимизированная функция остатка: s = s mod (2^p - 1)
void mersenne_mod(mpz_t s, unsigned long p, mpz_t temp_q) {
    // Пока количество битов в s больше p
    while (mpz_sizeinbase(s, 2) > p) {
        mpz_tdiv_q_2exp(temp_q, s, p); // Побитовый сдвиг вправо на p битов (верхние биты)
        mpz_tdiv_r_2exp(s, s, p);     // Оставляем только нижние биты
        mpz_add(s, s, temp_q);        // s = нижние + верхние
    }
}

// Проверка на простоту Мерсенна
unsigned long test_mersenne_lucas_lehmer(unsigned long p, mpz_t *out_primes) {

    if (p < 2) return 0;

    bool *is_prime = first_primes_ML(p);
    if (!is_prime) return 0;
    if (!is_prime[p]) {
        free(is_prime);
        return 0;
    }
    free(is_prime);

    // Инициализация переменных для проверки
    mpz_t s, M, temp_q;
    mpz_inits(s, M, temp_q, NULL);

    // M = 2^p - 1
    mpz_set_ui(s, 1);
    mpz_mul_2exp(M, s, p);
    mpz_sub_ui(M, M, 1);

    // Отсев ТОЧНО не чисел Мерсенна
    // Делители Mp всегда имеют вид 2*k*p + 1
    for (unsigned long k = 1; k <= 1000; k++) {
        unsigned long d = 2 * k * p + 1;
        if (mpz_divisible_ui_p(M, d)) return 0;
    }

    // Тест Люка-Лемера
    mpz_set_ui(s, 4); // Инициализация s = 4
    unsigned long iterations = p - 2;
    for (unsigned long i = 0; i < iterations; i++) { // Повторяем p-2 раз
        mpz_mul(s, s, s); // Возводим s в квадрат
        mpz_sub_ui(s, s, 2); // s = s^2 - 2
        mersenne_mod(s, p, temp_q); // Применяем оптимизированную функцию взятия остатка
    }

    if (mpz_cmp(s, M) >= 0) mpz_sub(s, s, M); // Нормализуем s, если оно больше или равно M

    // Если s == 0, то число Мерсенна простое
    unsigned long count = 0;
    if (mpz_cmp_ui(s, 0) == 0) {
        mpz_set(out_primes[0], M); // Записываем простое число в нулевой индекс
        count = 1;
    }
    mpz_clears(s, M, temp_q, NULL); // Освобождаем память
    return count;
}

// АЛГОРИТМ СОЛОВЕЯ ШТРАССЕНА
// Функция теста Соловея-Штрассена
PrimeStatus test_solovay_strassen(mpz_t n, int iterations) {
    // Базовые проверки
    if (mpz_cmp_ui(n, 1) <= 0) return 0;
    if (mpz_cmp_ui(n, 3) <= 0) return 2;
    if (mpz_even_p(n)) return 0;

    // Проверка по таблице малых простых чисел
    for (int i = 0; i < NUM_PROS_TABL; i++) {
        if (mpz_cmp_ui(n, prostie[i]) == 0) return 1;
        if (mpz_divisible_ui_p(n, prostie[i])) return 0;
    }

    int is_probably_prime = 1;

    // Параллельные вычисления
    #pragma omp parallel
    {
        // Подготовка переменных
        mpz_t n_minus_1, exp, res, mod_jacobi, g, a;
        mpz_inits(n_minus_1, exp, res, mod_jacobi, g, a, NULL);
        mpz_sub_ui(n_minus_1, n, 1);
        mpz_tdiv_q_2exp(exp, n_minus_1, 1); // Вычисляем exp = (n-1)/2

        // Распараллеливание итераций между потоками
        #pragma omp for schedule(static)
        for (int i = 0; i < iterations; i++) {
            if (!is_probably_prime) continue;

            mpz_set_ui(a, prostie[i]); // Выбираем простое число a из таблицы

            if (mpz_cmp(a, n) >= 0) continue;

            // Проверка НОД(a, n). Если > 1, число составное
            mpz_gcd(g, a, n);
            if (mpz_cmp_ui(g, 1) > 0) {
                is_probably_prime = 0;
                continue;
            }

            // Вычисляем левую часть: res = a^((n-1)/2) mod n
            mpz_powm(res, a, exp, n);

            // Вычисляем правую часть: символ Якоби (a/n)
            int jacobi = mpz_jacobi(a, n); // Может быть -1/0/1

            // Приводим символ Якоби к положительному остатку по модулю n
            // Если jacobi = -1, то mod_jacobi = n - 1
            if (jacobi == -1) mpz_set(mod_jacobi, n_minus_1);
            else mpz_set_si(mod_jacobi, jacobi);

            // Сравниваем: a^((n-1)/2) mod n != (a/n) mod n
            if (mpz_cmp(res, mod_jacobi) != 0) is_probably_prime = 0; // Если не равны, число составное
        }
        mpz_clears(n_minus_1, exp, res, mod_jacobi, g, a, NULL); // Освобождаем память
    }
    return is_probably_prime;
}

PrimeStatus test_trivial_division(mpz_t n) {
    if (mpz_cmp_ui(n, 2) < 0) return 0;
    if (mpz_cmp_ui(n, 2) == 0) return 2;
    if (mpz_even_p(n)) return 0;

    mpz_t limit, i;
    mpz_inits(limit, i, NULL);
    mpz_sqrt(limit, n);

    mpz_set_ui(i, 3);
    char is_prime = 1;

    while (mpz_cmp(i, limit) <= 0) {
        if (mpz_divisible_p(n, i)) {
            is_prime = 0;
            break;
        }
        mpz_add_ui(i, i, 2);
    }

    mpz_clears(limit, i, NULL);
    return is_prime;
}

void prime_lib_cleanup() {
    free(prostie);
}
