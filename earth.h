#pragma once
#include "comunity.h"

double evaluate_alfa(double border_alfa, double border_gamma) {
    if(IS_EQL(cos(border_alfa), 1) && IS_EQL(cos(border_gamma), 1))
        throw invalid_argument("Некорректное направление на край Земли (нулевое направление, вырожденный случай)");
    if(!(cos(border_alfa) > 0) || !(cos(border_gamma) > 0))
        throw invalid_argument("Некорректное направление на край Земли (не в переднем полупространстве датчика)");

    double tan_beta = sqrt(SQR(tan(border_alfa)) + SQR(tan(border_gamma)));
    double beta2 = atan(tan_beta) + asin(EARTH_RAD / ORBIT_RAD);

    if(beta2 < M_PI_2 && beta2 > 0) {
        return atan(tan(border_alfa) * tan(beta2) / tan_beta);
    }
    if(beta2 < M_PI && beta2 > M_PI_2) {
        return M_PI - atan(tan(border_alfa) * tan(M_PI - beta2) / tan_beta);
    }

    throw invalid_argument("Ой, плохо, плохо!");
}

double evaluate_gamma(double border_alfa, double border_gamma) {
    return evaluate_alfa(border_gamma, border_alfa);
}

class Earth {
public:
    // Условие стабильной работы: косинус обоих углов больше 0!
    // Кроме того, запрещен случай, когда оба косинуса равны 1
    // Особенность СК
    Earth(double _border_alfa, double _border_gamma)
        : border_alfa(_border_alfa)
        , border_gamma(_border_gamma)
        , alfa(evaluate_alfa(_border_alfa, _border_gamma))
        , gamma(evaluate_gamma(_border_alfa, _border_gamma))
        , position(Norm(dir_to_vec(alfa, gamma), ORBIT_RAD))
    {}

    // Выделяем направление на край (ближнюю точку) Земли
    const double border_alfa;
    const double border_gamma;

    // Направление на центр Земли
    const double alfa;
    const double gamma;

    // Координаты центра Земли
    const Vec position;
};
