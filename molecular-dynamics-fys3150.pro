TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += optimize_full
TARGET = mol_dyn

QMAKE_CXX += -g
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE *= -O3

SOURCES += main.cpp \
    atom.cpp \
    system.cpp \
    velocityverlet.cpp \
    math/vec3.cpp \
    io.cpp \
    lennardjones.cpp \
    statisticssampler.cpp \
    unitconverter.cpp

HEADERS += \
    atom.h \
    system.h \
    velocityverlet.h \
    math/vec3.h \
    math/random.h \
    io.h \
    lennardjones.h \
    statisticssampler.h \
    unitconverter.h

