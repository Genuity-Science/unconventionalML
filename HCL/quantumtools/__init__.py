# quantumtools package, creates functions that are reused
# should add this package to Python path for it to work properly, and be
# importable from anywhere
#
# high to low level:
# gates > qtools > unitary > accessories
# if lower calls higher, it should do so within the function only,
# else you get error from circular references
