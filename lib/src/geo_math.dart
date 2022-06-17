
import 'dart:math';

import 'package:geographiclib/src/pair.dart';

import 'math.dart';

class GeoMath {


  static const int digits = 53;

  static double sq(double x) => x * x;

  static double atanh(double x) {
    double y = x.abs();
    y = log1p(2 * y / (1 - y)) / 2;
    return x > 0 ? y : (x < 0 ? -y : x);
  }

  static void norm(Pair p, double sinx, double cosx) {
    double r = hypot(sinx, cosx) as double;
    p.first = sinx / r;
    p.second = cosx / r;
  }

  static void sum(Pair p, double u, double v) {
    double s = u + v;
    double up = s - v;
    double vpp = s - up;
    up -= u;
    vpp -= v;
    double t = s != 0 ? 0.0 - (up + vpp) : s;
    p.first = s;
    p.second = t;
  }

  static double polyval(int N, List<double> p, int s, double x ) {
    double y = N < 0 ? 0 : p[s++];
    while (--N >= 0) {
      y = y * x + p[s++];
    }
    return y;
  }

  static double angRound(double x) {
    final double z = 1 / 16.0;
    double y = x.abs();
    y = y < z ? z - (z - y) : y;
    return copySign(y, x);
  }

  static double angNormalize(double x) {
    double y = iEEERemainder(x, 360.0);
    return y.abs() == 180 ? copySign(180.0, x) : y;
  }

  static double latFix(double x) {
    return x.abs() > 90 ? double.nan : x;
  }

  static void angDiff(Pair p, double x, double y) {
    sum(p, iEEERemainder(-x, 360.0), iEEERemainder(y, 360.0));
    sum(p, iEEERemainder(p.first, 360.0), p.second);
    if (p.first == 0 || p.first.abs() == 180) {
      p.first = copySign(p.first, p.second == 0 ? y - x : -p.second);
    }
  }

  static void sincosd(Pair p, double x) {
    double r; int q;
    r = x % 360.0;
    q = (r / 90).round();
    r -= 90 * q;
    r = toRadian(r);
    double s = sin(r), c = cos(r);
    double sinx = double.nan, cosx = double.nan;
    switch (q & 3) {
      case 0:
        sinx = s;
        cosx = c;
        break;
      case 1:
        sinx = c;
        cosx = -s;
        break;
      case 2:
        sinx = -s;
        cosx = -c;
        break;
      default:
        sinx = -c;
        cosx = s;
        break;
    }
    if (sinx == 0) sinx = copySign(sinx, x);
    p.first = sinx;
    p.second = 0.0 + cosx;
  }

  static void sincosde(Pair p, double x, double t) {
    double r; int q;
    q = ( x / 90 ).round();
    r = x - 90 * q;
    r = toRadian(angRound(r + t));
    double s = sin(r), c = cos(r);
    double sinx = double.nan, cosx = double.nan;
    switch (q & 3) {
      case 0:
        sinx = s;
        cosx = c;
        break;
      case 1:
        sinx = c;
        cosx = -s;
        break;
      case 2:
        sinx = -s;
        cosx = -c;
        break;
      default:
        sinx = -c;
        cosx = s;
        break;
    }

    if (sinx == 0) sinx = copySign(sinx, x);
    p.first = sinx;
    p.second = 0.0 + cosx;
  }

  static double atan2d(double y, double x) {
    int q = 0;
    if (y.abs() > x.abs()) {
      double t;
      t = x;
      x = y;
      y = t;
      q = 2;
    }
    if (x < 0) {
      x = -x;
      ++q;
    }

    double ang = toDegree(atan2(y, x));
    switch (q) {
      case 1:
        ang = copySign(180.0, y) - ang;
        break;
      case 2:
        ang = 90 - ang;
        break;
      case 3:
        ang = -90 + ang;
        break;
      default:
        break;
    }
    return ang;
  }
}
