.onLoad = function(libname, pkgname) {
  .jpackage(pkgname, lib.loc = libname)
  .jaddClassPath(system.file("java",package="CBSr"))
}
