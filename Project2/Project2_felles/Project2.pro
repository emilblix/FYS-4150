TEMPLATE = app
CONFIG += console
CONFIG -= qt

win32 {
    INCLUDEPATH = C:\Dropbox\Emil\FYS4150\armadillo\include
    PRE_TARGETDEPS += "C:\Dropbox\Emil\FYS4150\armadillo\examples\lib_win64\lapack_win64_MT.lib"
    PRE_TARGETDEPS += "C:\Dropbox\Emil\FYS4150\armadillo\examples\lib_win64\blas_win64_MT.lib"
}

unix {
    LIBS += -larmadillo
}

release {
    DEFINES += ARMA_NO_DEBUG
}

SOURCES += main.cpp \
    Jacobi_rotate.cpp \
    offdiag.cpp

HEADERS += \
    Jacobi_rotate.h \
    offdiag.h

