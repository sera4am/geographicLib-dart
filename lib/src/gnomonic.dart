

import 'dart:core';
import 'dart:math';

import 'package:geographiclib/src/pair.dart';
import 'geo_math.dart';
import 'geodesic.dart';
import 'geodesic_data.dart';
import 'geodesic_mask.dart';
import 'math.dart';

class Gnomonic {

  //  private static final double eps_ = 0.01 * Math.sqrt(Math.ulp(1.0));
  static const double eps_ = 1.4901161193847657E-10;
  static const int numit_ = 10;
  late Geodesic _earth;
  late double _a, _f;

  Gnomonic(Geodesic earth) {
    _earth = earth;
    _a = _earth.EquatorialRadius();
    _f = _earth.Flattening();
  }

  GnomonicData forward(double lat0, double lon0, double lat, double lon)
  {
    GeodesicData inv =
    _earth.Inverse(lat0, lon0, lat, lon,
        GeodesicMask.azimuth | GeodesicMask.geodesicscale |
        GeodesicMask.reducedlength);
    GnomonicData fwd = GnomonicData(lat0, lon0, lat, lon, double.nan, double.nan, inv.azi2, inv.M12);

    if (inv.M12 > 0) {
      double rho = inv.m12 / inv.M12;
      Pair p = Pair(double.nan, double.nan);
      GeoMath.sincosd(p, inv.azi1);
      fwd.x = rho * p.first;
      fwd.y = rho * p.second;
    }

    return fwd;
  }

  GnomonicData? reverse(double lat0, double lon0, double x, double y) {
    GnomonicData rev =
    GnomonicData(lat0, lon0, double.nan, double.nan, x, y, double.nan, double.nan);

    double azi0 = GeoMath.atan2d(x, y);
    double rho = hypot(x, y);
    double s = _a * atan(rho / _a);
    bool little = rho <= _a;

    if (!little)
      rho = 1 / rho;

    GeodesicLine line =
    _earth.Line(lat0, lon0, azi0, GeodesicMask.latitude
    | GeodesicMask.longitude | GeodesicMask.azimuth
    | GeodesicMask.distanceIn | GeodesicMask.reducedlength
    | GeodesicMask.geodesicscale);

    int count = numit_, trip = 0;
    late GeodesicData pos;

    while (count-- > 0) {
      pos =
          line.Position(s, GeodesicMask.longitude | GeodesicMask.latitude
          | GeodesicMask.azimuth | GeodesicMask.distanceIn
          | GeodesicMask.reducedlength
          | GeodesicMask.geodesicscale);

      if (trip > 0) {
        break;
      }

      double ds =
      little ? ((pos.m12 / pos.M12) - rho) * pos.M12 * pos.M12
          : (rho - (pos.M12 / pos.m12)) * pos.m12 * pos.m12;
      s -= ds;

      if (ds.abs() <= eps_ * _a) {
        trip++;
      }
    }

    if (trip == 0) {
      return rev;
    }

    rev.lat = pos.lat2;
    rev.lon = pos.lon2;
    rev.azi = pos.azi2;
    rev.rk = pos.M12;

    return rev;
  }

  /**
   * @return <i>a</i> the equatorial radius of the ellipsoid (meters).  This is
   *   the value inherited from the Geodesic object used in the constructor.
   **********************************************************************/
  public double EquatorialRadius() { return _a; }

  /**
   * @return <i>f</i> the  flattening of the ellipsoid.  This is
   *   the value inherited from the Geodesic object used in the constructor.
   **********************************************************************/
  public double Flattening() { return _f; }
}
