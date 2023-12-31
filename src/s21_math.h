#ifndef SRC_S21_MATH_H_
#define SRC_S21_MATH_H_

#include <ctype.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#define S21_INF 1.0 / 0.0
#define S21_NAN 0.0 / 0.0

#define is_fin(x) __builtin_isfinite(x)
#define is_nan(x) __builtin_isnan(x)
#define is_inf(x) __builtin_isinf(x)

#define S21_PI 3.14159265358979323846264338327950288
#define S21_EXP 2.7182818284590452353602874713526624
#define S21_EPS 1e-9

// FUNCTIONS
int s21_abs(int x);
long double s21_acos(double x);
long double s21_asin(double x);
long double s21_atan(double x);
long double s21_ceil(double x);
long double s21_cos(double x);
long double s21_exp(double x);
long double s21_fabs(double x);
long double s21_floor(double x);
long double s21_fmod(double x, double y);
long double s21_log(double x);
long double s21_pow(double base, double exp);
long double s21_sin(double x);
long double s21_sqrt(double x);
long double s21_tan(double x);

// UTILITIES
long double s21_fmax(double a, double b);
long double s21_factorial(int N);
long long int s21_abs_long_int(long long int x);
long double s21_pow_helper(double base, double exp);
long double s21_pow_int(long double base, long long int exp_int);

#endif
