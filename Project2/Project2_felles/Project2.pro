TEMPLATE = app
CONFIG += console
CONFIG -= qt

win32 {
    INCLUDEPATH = C:\Dropbox\Emil\FYS4150\armadillo\include
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

