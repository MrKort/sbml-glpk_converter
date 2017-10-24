QT += core
QT -= gui

CONFIG += c++11

TARGET = source_code
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp
LIBS += -L/usr/lib64 -lsbml
LIBS += -lglpk
