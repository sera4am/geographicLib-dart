
class GeodesicMask {

  static const int capNone = 0;
  static const int capC1 = 1<<0;
  static const int capC1p = 1<<1;
  static const int capC2 = 1<<2;
  static const int capC3 = 1<<3;
  static const int capC4 = 1<<4;
  static const int capAll = 0x1F;
  static const int capMask = capAll;
  static const int outAll = 0x7F80;
  static const int outMask = 0xFF80;

  static const int none = 0;
  static const int latitude = 1<<7 | capNone;
  static const int longitude = 1<<8 | capC3;
  static const int azimuth = 1<<9 | capNone;
  static const int distance = 1<<10 | capC1;
  static const int standard = latitude | longitude | azimuth | distance;

  static const int distanceIn = 1<<11 | capC1 | capC1p;
  static const int reducedlength = 1<<12 | capC1 | capC2;
  static const int geodesicscale = 1<<13 | capC1 | capC2;

  static const int area = 1<<14 | capC4;
  static const int all = outAll | capAll;
  static const int longUnroll = 1<<15;


}
