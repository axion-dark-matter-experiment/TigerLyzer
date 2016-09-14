QT += core
QT -= gui

TARGET = NouveauAnalysis
CONFIG += console
CONFIG -= app_bundle

CONFIG += c++11

LIBS += -L/usr/local/lib -L/usr/lib -lboost_iostreams -lboost_system -lboost_filesystem

QMAKE_CXXFLAGS+= -fopenmp
#QMAKE_CXXFLAGS+= -O3
QMAKE_LFLAGS +=  -fopenmp

TEMPLATE = app

SOURCES += main.cpp \
    flatfileinterface.cpp \
    spectrum.cpp \
    spectrumfilter.cpp \
    plotter.cpp

HEADERS += \
    flatfileinterface.h \
    spectrum.h \
    spectrumfilter.h \
    plotter.h

