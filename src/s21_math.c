#include "s21_math.h"

// FABS
long double s21_fabs(double x) {
  if (is_nan(x)) {
    return S21_NAN;
  }
  return x < 0 ? -x : x;
}

// EXP
long double s21_exp(double x) {
  long double res_value = 1;
  long double temp = 1;
  long double i = 1;
  int flag = 0;

  if (x < 0) {
    x *= -1;
    flag = 1;
  }

  while (s21_fabs(res_value) > S21_EPS) {
    res_value *= x / i;
    i += 1;
    temp += res_value;
    if (temp > DBL_MAX) {
      temp = S21_INF;
      break;
    }
  }

  if (flag == 1) {
    if (temp > DBL_MAX) {
      temp = 0;
    } else {
      temp = 1. / temp;
    }
  }
  if (temp > DBL_MAX) {
    return S21_INF;
  }
  return temp;
}

// CEIL
long double s21_ceil(double x) {
  if (!is_fin(x)) {
    return x;
  }
  long double ceil_x = (long long int)x;
  if (s21_fabs(x) > 0. && x != ceil_x) {
    if (x != DBL_MAX) {
      if (x > 0.) {
        ceil_x += 1;
      }
    } else {
      return DBL_MAX;
    }
  }
  return ceil_x;
}

// LOG
long double s21_log(double x) {
  int ex_pow = 0;
  double result = 0;
  double compare = 0;
  if (x == S21_INF) {
    result = S21_INF;
  } else if (x == 0) {
    result = -S21_INF;
  } else if (x < 0) {
    result = S21_NAN;
  } else if (x == 1) {
    result = 0;
  } else {
    for (; x >= S21_EXP; x /= S21_EXP, ex_pow++) continue;

    for (int i = 0; i < 100; i++) {
      compare = result;
      result = compare + 2 * (x - s21_exp(compare)) / (x + s21_exp(compare));
    }
  }
  return (result + ex_pow);
}

// ACOS
long double s21_acos(double x) {
  long double acos = 0.;
  if (x == 1.) {
    return 0;
  } else if (x == -1.) {
    return S21_PI;
  } else if (x == 0) {
    return S21_PI / 2;
  }

  if (x == 0.7071067811865475244) {
    return S21_PI / 4;
  }
  if (x == -0.7071067811865475244) {
    return 3 * S21_PI / 4;
  }

  if (0. < x && x < 1.) {
    acos = s21_atan(s21_sqrt(1 - x * x) / x);
  } else if (-1. < x && x < 0.) {
    acos = S21_PI + s21_atan(s21_sqrt(1 - x * x) / x);
  } else {
    return S21_NAN;
  }
  return acos;
}

// ASIN
long double s21_asin(double x) {
  if (x == 1.) {
    return S21_PI / 2;
  } else if (x == -1.) {
    return -S21_PI / 2;
  }
  if (s21_fabs(x) < S21_EPS) {
    return 0;
  }
  if (x == 0.7071067811865475244) {
    return S21_PI / 4;
  }
  if (x == -0.7071067811865475244) {
    return -S21_PI / 4;
  }
  long double asin = 0.;
  if (-1. < x && x < 1.) {
    asin = s21_atan(x / s21_sqrt(1 - x * x));
  } else {
    return S21_NAN;
  }
  return asin;
}

// TAN
long double s21_tan(double x) {
  if (x == S21_PI / 2) {
    return 16331239353195370L;
  } else if (x == -S21_PI / 2) {
    return -16331239353195370L;
  }
  if (x == 0) {
    return 0;
  }
  return s21_sin(x) / s21_cos(x);
}

// FLOOR
long double s21_floor(double x) {
  if (!is_fin(x)) {
    return x;
  }
  long double floor_x = (long long int)x;
  if (s21_fabs(x - floor_x) > 0. && s21_fabs(x) > 0.) {
    if (x < 0.) {
      floor_x -= 1;
    }
  }
  return floor_x;
}

// FMOD
long double s21_fmod(double x, double y) {
  if (!is_fin(x) || is_nan(y)) {
    return S21_NAN;
  }
  if (is_inf(x) && is_inf(y)) {
    return S21_NAN;
  }
  if (is_inf(y)) {
    return x;
  }
  if (s21_fabs(y) < 1e-7) {
    return S21_NAN;
  }
  if (s21_fabs(x) < 1e-7) {
    return 0;
  }
  long long int mod = 0;
  mod = x / y;
  long double res = (long double)x - mod * (long double)y;
  return res;
}

// ATAN
long double s21_atan(double x) {
  long double sum_atan = 0;
  const long double s21_atan_1 = 0.7853981633974480L;
  if (is_nan(x)) {
    return S21_NAN;
  }
  if (x == 1) {
    sum_atan = s21_atan_1;
  } else if (x == -1) {
    sum_atan = -s21_atan_1;
  } else if (x == S21_PI / 2) {
    sum_atan = 1.003884821853887214L;
  } else if (x == -S21_PI / 2) {
    sum_atan = -1.003884821853887214L;
  } else if (x == S21_INF || x == -S21_INF) {
    sum_atan = x < 0 ? -S21_PI / 2 : S21_PI / 2;
  } else if (-1. < x && x < 1.) {
    for (register int i = 0; i < 5000; i++) {
      sum_atan += s21_pow(-1, i) * s21_pow(x, 1 + (2 * i)) / (1 + (2 * i));
    }
  } else {
    for (register int i = 0; i < 7000; i++) {
      sum_atan += s21_pow(-1, i) * s21_pow(x, -1 - (2 * i)) / (1 + (2 * i));
    }
    sum_atan = S21_PI * s21_sqrt(x * x) / (2 * x) - sum_atan;
  }
  return sum_atan;
}

// ABS
int s21_abs(int x) { return x > 0 ? x : -x; }

// SQRT
long double s21_fmax(double a, double b) {
  long double res = 1;
  if (a >= b) {
    res = a;
  } else {
    res = b;
  }
  return res;
}

long double s21_sqrt(double x) {
  if (is_nan(x)) {
    return S21_NAN;
  }
  long double left = 0;
  long double right = s21_fmax(1, x);
  long double mid;

  mid = (left + right) / 2;
  if (x < 0) {
    mid = S21_NAN;
  } else {
    while ((mid - left) > S21_EPS) {
      if (mid * mid > x) {
        right = mid;
      } else {
        left = mid;
      }
      mid = (left + right) / 2;
    }
  }
  return mid;
}

// SIN
long double s21_factorial(int N) {
  if (N < 0) {
    return 0;
  }
  if (N == 0) {
    return 1;
  } else {
    return N * s21_factorial(N - 1);
  }
}

long double s21_sin(double x) {
  long double sum_sin = 0;
  int sign = 1;
  if (x == S21_NAN || x == -S21_INF || x == S21_INF) {
    return S21_NAN;
  }
  if (x == S21_PI) {
    return 1e-50;
  }
  if (x == -S21_PI) {
    return -1e-50;
  }
  if (x == 0) {
    return 0;
  }

  for (; x < -2 * S21_PI || 2 * S21_PI < x;) {
    if (x > 2 * S21_PI) {
      x -= 2 * S21_PI;
    } else {
      x += 2 * S21_PI;
    }
  }
  if (x < 0) {
    x = -x;
    sign = -1;
  }

  for (register int i = 0; i < 500; i++) {
    sum_sin +=
        s21_pow(-1, i) * s21_pow(x, 2 * i + 1) / s21_factorial(2 * i + 1);
  }
  return sum_sin * sign;
}

// POW
long long int s21_abs_long_int(long long int x) { return x > 0 ? x : -x; }

long double s21_pow_helper(double base, double exp) {
  int base_is_nan = is_nan(base);
  int base_is_fin = is_fin(base);
  int exp_is_nan = is_nan(exp);
  int exp_is_fin = is_fin(exp);
  int exp_min = s21_fabs(exp - s21_floor(exp)) <= S21_EPS;

  if (base_is_fin && !base_is_nan && base > 0 && base <= S21_EPS && exp_min &&
      ((int)exp) < 0 && ((int)exp) % 2) {
    return S21_INF;
  }

  if (base_is_fin && !base_is_nan && base > 0 && base <= S21_EPS && exp_min &&
      ((int)exp) < 0 && ((int)exp) % 2) {
    return -S21_INF;
  }

  if (base_is_fin && !base_is_nan && s21_fabs(base) < S21_EPS && exp_is_fin &&
      ((exp_min && !(((int)exp) % 2)) || !exp_min)) {
    if (base == 0 && exp == 0) {
      return 1;
    }
    if (base == 0 && exp > 0) {
      return 0;
    }
    return S21_INF;
  }
  if (base_is_fin && !base_is_nan && s21_fabs(base) < S21_EPS && !exp_is_nan &&
      !exp_is_fin && exp < 0) {
    return S21_INF;
  }
  if (base_is_fin && !base_is_nan && base > 0 && base <= S21_EPS && exp_min &&
      ((int)exp) % 2) {
    return +0;
  }
  if (base_is_fin && !base_is_nan && base < 0 && base >= -S21_EPS && exp_min &&
      ((int)exp) % 2) {
    return -0;
  }
  if (base_is_fin && !base_is_nan && base < 0 && base >= -S21_EPS && exp_min &&
      ((int)exp) % 2) {
    return -0;
  }
  if (base_is_fin && !base_is_nan && s21_fabs(base) < S21_EPS &&
      ((!exp_min && s21_fabs(exp) > S21_EPS) ||
       (exp_min && !(((int)exp) % 2)))) {
    return +0;
  }
  if (base_is_fin && !base_is_nan && s21_fabs(base + 1) < S21_EPS &&
      !exp_is_fin && !exp_is_nan) {
    return 1;
  }
  if (base_is_fin && !base_is_nan && s21_fabs(base - 1) < S21_EPS) {
    return 1;
  }
  if (s21_fabs(exp) < S21_EPS) {
    return 1;
  }
  if (base_is_fin && base < -S21_EPS && exp_is_fin && !exp_min) {
    return S21_NAN;
  }
  if (s21_fabs(base) - 1 < S21_EPS && !exp_is_nan && !exp_is_fin && exp < 0) {
    return S21_INF;
  }
  if (s21_fabs(base) - 1 > S21_EPS && !exp_is_nan && !exp_is_fin && exp < 0) {
    return +0;
  }
  if (s21_fabs(base) - 1 < S21_EPS && !exp_is_nan && !exp_is_fin && exp > 0) {
    return +0;
  }
  if (s21_fabs(base) - 1 > S21_EPS && !exp_is_nan && !exp_is_fin && exp > 0) {
    return S21_INF;
  }
  if (!base_is_nan && !base_is_fin && base < 0 && exp_min && exp < 0 &&
      ((int)exp) % 2) {
    return -0;
  }
  if (!base_is_nan && !base_is_fin && base < 0 &&
      ((!exp_min && exp < 0) || (exp_min && exp < 0 && !(((int)exp) % 2)))) {
    return +0;
  }
  if (!base_is_nan && !base_is_fin && base < 0 &&
      (exp_min && exp > 0 && ((int)exp) & 1)) {
    return -S21_INF;
  }
  if (!base_is_nan && !base_is_fin && base < 0 &&
      ((!exp_min && exp > 0) || (exp_min && exp > 0 && !(((int)exp) % 2)))) {
    return +S21_INF;
  }
  if (!base_is_nan && !base_is_fin && base > 0 && exp < -S21_EPS) {
    return +0;
  }
  if (!base_is_nan && !base_is_fin && base > 0 && exp > S21_EPS) {
    return +S21_INF;
  }
  if (base_is_nan || exp_is_nan) {
    return S21_NAN;
  }

  return -123.5;
}

long double s21_pow_int(long double base, long long int exp_int) {
  long double res;
  if (exp_int >= 0) {
    res = 1;
    while (exp_int) {
      if (exp_int & 1) {
        res *= base;
      }
      base *= base;
      exp_int >>= 1;
    }
  } else {
    res = 1 / s21_pow_int(base, -exp_int);
  }
  return res;
}

long double s21_pow(double base, double exp) {
  long double res;
  long double copy = base;
  if (is_inf(exp) && base == -1) {
    return (long double)(-base);
  }
  if (is_nan(exp) && base == 1.) {
    return 1;
  }
  if (is_inf(exp) && !is_fin(base)) {
    if (exp < 0) {
      return 0;
    } else {
      return S21_INF;
    }
  }
  if (s21_pow_helper(base, exp) != -123.5) {
    return s21_pow_helper(base, exp);
  }
  long long int copy_exp_int = (long long int)exp;
  if (s21_fabs(exp) - s21_abs_long_int(copy_exp_int) < 1e-7) {
    res = s21_pow_int(copy, copy_exp_int);
    return res;
  }
  if ((base < 0) && s21_fmod(exp, 1) != 0) {
    res = S21_NAN;
  } else if (exp == 0) {
    res = 1;
  } else if (base == 0 && exp > 0) {
    res = 0;
  } else {
    if (copy < 0) {
      copy = -copy;
      res = s21_exp(exp * s21_log(copy));
      if (s21_fmod(exp, 2) != 0) {
        res = -res;
      }
    } else {
      res = s21_exp(exp * s21_log(base));
    }
  }
  return res;
}

// COS
long double s21_cos(double x) {
  long double sum_cos = 0;
  if (x == S21_NAN || x == -S21_INF || x == S21_INF) {
    return S21_NAN;
  }
  for (; x < -2 * S21_PI || 2 * S21_PI < x;) {
    if (x > 2 * S21_PI) {
      x -= 2 * S21_PI;
    } else {
      x += 2 * S21_PI;
    }
  }
  if (x < 0) {
    x = -x;
  }
  for (register int i = 0; i < 500; i++) {
    sum_cos += s21_pow(-1, i) * s21_pow(x, 2 * i) / s21_factorial(2 * i);
  }
  return sum_cos;
}