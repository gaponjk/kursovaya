QT = core

CONFIG += c++17 cmdline

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
        coefficients.cpp \
        gauss_seidel.cpp \
        grid_level.cpp \
        main.cpp \
        prolongation.cpp \
        residual.cpp \
        restriction.cpp \
        test_problems.cpp \
        vcycle.cpp

TRANSLATIONS += \
    kursovaya_ru_BY.ts
CONFIG += lrelease
CONFIG += embed_translations

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    coefficients.h \
    gauss_seidel.h \
    grid_level.h \
    prolongation.h \
    residual.h \
    restriction.h \
    test_problems.h \
    vcycle.h
