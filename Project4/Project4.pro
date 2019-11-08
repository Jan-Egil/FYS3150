TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -fopenmp
QMAKE_CXXFLAGS += -fopenmp -O3


SOURCES += \
        main.cpp
