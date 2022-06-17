
import 'dart:math';

/// @see https://pub.dev/documentation/extended_math/latest/extended_math/hypot.html
double hypot(double x, double y) {
double first = x.abs();
double second = y.abs();

if (y > x) {
first = y.abs();
second = x.abs();
}
if (first == 0.0) {
return second;
}

final t = second / first;
return first * sqrt(1 + t * t);
}

/// @see https://pub.dev/documentation/grizzly_distuv/latest/math/log1p.html
double log1p(double x) {
if (x.isInfinite && !x.isNegative) {
return x;
}
final double u = 1 + x;
final double d = u - 1;

if (d == 0) {
return x;
}
return log(u) * x / d;
}

double copySign(double m, double s) {
  if (s == 0 || s.isNaN || m.sign == s.sign) {
    return m;
  }
  return -m;
}

bool isSameSign(double x, double y) {
  return copySign(x, y) == x;
}

double iEEERemainder(double m, double s) {
  if (s == 0.0) { return double.nan; }
  return m - (s * (m / s).roundToDouble());
}

double toRadian(double x) => x * pi / 180;

double toDegree(double x) => x* 180 / pi;

num cbrt(num x) => pow(x, 1/3);
