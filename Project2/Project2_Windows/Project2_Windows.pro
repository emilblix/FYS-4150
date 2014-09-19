TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
win32 {
    INCLUDEPATH = C:\Dropbox\Emil\FYS4150\armadillo\include
}
unix {

}
SOURCES += \
    main.cpp \
    Jacobi_rotate.cpp \
    offdiag.cpp

HEADERS += \
    Jacobi_rotate.h \
    offdiag.h

