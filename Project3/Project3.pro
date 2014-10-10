TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    runge_kutta_4.cpp \
    vec3.cpp \
    solarsystem.cpp \
    celestialbody.cpp

HEADERS += \
    runge_kutta_4.h \
    vec3.h \
    solarsystem.h \
    celestialbody.h

