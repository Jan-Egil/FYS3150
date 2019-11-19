TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -fopenmp
LIBS += -larmadillo

QMAKE_CXXFLAGS += -O3


SOURCES += \
        main.cpp
