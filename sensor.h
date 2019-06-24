#pragma once
#include "comunity.h"

class Sensor {
public:
    // Предельный рабочий угол по большой оси в радианах
    const double alfa_range = ALFA;
    // Геометрический параметр L/h
    const double param_k_a = tan(alfa_range);
    // Геометрический параметр L/2h
    const double param_k_g = param_k_a / 2;
    // Предельный рабочий угол по малой оси в радианах
    const double gamma_range = atan(param_k_g);

    // Предельный угол попадания света в датчик
    // большая ось
    const double beta_range_a = atan(2 * param_k_a);
    // малая ось
    const double beta_range_g = atan(3 * param_k_g);

    // Высота и радиус орбиты
    const double orbit_rad = ORBIT_RAD;
    const double orbit_h = ORBIT_RAD - EARTH_RAD;

    /* Сигналы на фотодиодах
       Q21   Q22   Q23
       Q11   Q12   Q13
    */
    double Q11 = 0, Q21 = 0, Q12 = 0, Q22 = 0, Q13 = 0, Q23 = 0;
};
