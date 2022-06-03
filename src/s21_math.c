#include "s21_math.h"

int s21_abs(int num) {
    return num < 0 ? num * -1 : num;
}

long double s21_acos(double x) {
    return s21_P_2 - s21_asin(x);
}

long double s21_asin(double x) {
    return s21_atan(x / s21_sqrt(1 - (x * x)));
}

long double s21_atan(double x) {
    int neg = 0, per = 0, sp = 0;
    x < 0.l ? x = -x, neg = 1.l : 0;
    x > 1.l ? x = (1.l / x), per = 1 : 0;
    while (x > s21_P / 12.l) {
        sp++;
        x = ((x * s21_sqrt(3)) - 1.l) * (1.l / (x + s21_sqrt(3)));
    }
    x = x * ((0.55913709l / (1.4087812l + x * x)) + 0.60310579l - 0.05160454l * (x * x));
    while (sp > 0) {
        x += (s21_P_6);
        sp--;
    }
    x = per ? (s21_P_2) - x : x;
    x = neg ? -x : x;
    return x;
}

long double s21_ceil(double num) {
    int res = num;
    return res == num ? res : res < num ? (long double)res + 1.0 : (long double)res;
}

long double s21_cos(double x) {
    int z = s21_fabs(x) / s21_P, i = x < 0 ? -1 : 1;
    z > 1 ? x -= (s21_P * z) * i : 0;
    double y = 1;
    double s = -1;
    for (int i = 2; i <= 100; i += 2) {
        y += s * (my_pow(x, i) / s21_fact(i));
        s *= -1;
    }
    return z == 1 ? y : z % 2 == 0 ? y : -y;
}

long double s21_sin(double x) {
    int z = s21_fabs(x) / s21_P, i = x < 0 ? -1 : 1;
    z > 1 ? x -= (s21_P * z) * i : 0;
    double y = x;
    double s = -1;
    for (int i = 3; i <= 100; i += 2) {
        y += s * (my_pow(x, i) / s21_fact(i));
        s *= -1;
    }
    return z == 1 ? y : z % 2 == 0 ? y : -y;
}

long double s21_tan(double x) {
     return (s21_sin(x)/s21_cos(x));
}

long double my_exp(double x, long double s, long double n, long double a) {
    a = a * x / n;
    return s21_fabs(a) <= ACC ? s : my_exp(x, s + a, n + 1, a);
}

long double s21_exp(double x) {
    long double s = 1;
    long double n = 1;
    long double a = 1;
    return x != x ? S21_NAN : my_exp(x, s, n, a);
}

long double s21_fabs(double num) {
    return num < 0 ? num * -1 : num;
}

long double s21_floor(double x) {
    int y = x;
    return y == x ? y : y < x ? (long double)y : (long double)y - 1;
}

void mfp(double x, double y, double *buf) {
    while (*buf <= x) *buf += y;
    *buf -= y;
}

long double my_fmod(double x, double y, double buf) {
    if (x < 0) { buf = s21_fabs(buf);
        y = s21_fabs(y);
        x = s21_fabs(x);
        mfp(x, y, &buf);
        buf = (x - buf) * -1;
    } else {  buf = s21_fabs(buf);
        y = s21_fabs(y);
        x = s21_fabs(x);
        mfp(x, y, &buf);
        buf = x - buf;
    }
    return buf;
}

long double s21_fmod(double x, double y) {
    long double res = 0.0;
    if (x != x || y != y || y == 0.0) { res = S21_NAN;
    } else if (y > x && x > 0) { res = x;
    } else if (y == x) { res = 0.0;
    } else { double buf = y;
        res = my_fmod(x, y, buf);
    }
    return res;
}

long double s21_fact(long double x) {
    return x > 1.l ? x * s21_fact(x - 1.l) : 1.l;
}

long double my_pow(long double base, long int exp) {
    return exp > 1 ? base * my_pow(base, exp - 1) : base;
}

long double s21_log(double x) {
    long double gr = 0, temp = 0, result = 0, z = 2, m = x < 0 ? 1 : x == 0 ? 2 : 0;
    while (((x >= 10) || (x < 1 && x > 0)) && m == 0)
        x < 1 && x > 0 ? (x *= 10.l, gr -= 1) : (x *= 0.1, gr += 1);
    x /= 10.l;
    for (x--, result = x, temp = x; m == 0 && \
        s21_fabs(result) > ACC; result *= -x * (z - 1) / z, z++, temp += result) {}
    return m == 1 ? S21_NAN : m == 2 ? -S21_INF : temp + (gr + 1) * S21_LN10;
}

long double s21_pow(double base, double exp) {
    long double res = exp == 0 ? 1.l : exp == 1 ? base : base < 0 && exp != (int)exp ? S21_NAN : 0;
    int r = base < 0 && (int)exp % 2 != 0 ? -1 : 1;
    int g = base < 0 ? -1 : 1, t = exp < 0 ? -1 : 1;
    exp == (int)exp && exp > 0 ? res = my_pow(base, exp) : 0;
    if (res == 0) {
        res = r * s21_exp((exp * t) * s21_log(base * g));
        t < 0 ? res = 1.l / res : 0;
    }
    return res;
}

long double s21_sqrt(double x) {
    long double y1 = 0, n = 2.l, y2 = x > 0 ? 1.l : S21_NAN;
    for (; s21_fabs(y2 - y1) >= ACC && x > 0; y1 = y2, y2 = y1 + (x / my_pow(y1, n - 1.l) - y1) / n) {}
    return x == 0 ? x : y2;
}
