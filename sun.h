#pragma once
#include "comunity.h"

class Sun {
public:
    // Условие стабильной работы: либо косинус обоих углов больше 0, либо меньше!
    // Особенность СК
    Sun(double _alfa, double _gamma)
        : alfa(_alfa)
        , gamma(_gamma) {}

    // Направление на Солнце
    const double alfa;
    const double gamma;
    const Vec direction = dir_to_vec(alfa, gamma);
};
