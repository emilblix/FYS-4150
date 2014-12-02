#TEMPLATE = app
#CONFIG += console
#CONFIG -= app_bundle
#CONFIG -= qt

#SOURCES += main.cpp \
#    openglwindow.cpp \
#    trianglewindow.cpp

#HEADERS += \
#    openglwindow.h \
#    trianglewindow.h

include(openglwindow.pri)

SOURCES += \
    main.cpp \
    openglwindow.cpp \
    trianglewindow.cpp

target.path = $$[QT_INSTALL_EXAMPLES]/gui/openglwindow
INSTALLS += target

HEADERS += \
    openglwindow.h \
    trianglewindow.h
