
class GeographicErr implements Exception {

  final String msg;

  GeographicErr(this.msg);

  @override
  String toString() => "GeographicLib exception : $msg";
}
