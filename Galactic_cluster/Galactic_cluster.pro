TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    celestialbody.cpp \
    vec3.cpp \
    solarsystem.cpp \
    rk_4.cpp \
    verlet.cpp

HEADERS += \
    celestialbody.h \
    vec3.h \
    solarsystem.h \
    rk_4.h \
    verlet.h

