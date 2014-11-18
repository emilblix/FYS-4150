TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    celestialbody.cpp \
    vec3.cpp \
    rk_4.cpp \
    verlet.cpp \
    cluster.cpp \
    rk4_adaptive.cpp

HEADERS += \
    celestialbody.h \
    vec3.h \
    rk_4.h \
    verlet.h \
    cluster.h \
    rk4_adaptive.h

