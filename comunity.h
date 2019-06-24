#pragma once
#include <cmath>
#include <iostream>

using namespace std;

#define SQR(x) ((x) * (x))
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define TO_RAD(x) ((x) * M_PI / 180)
#define TO_GRAD(x) ((x) * 180 / M_PI)
#define IS_EQL(x, y) (((x) - (y)) < 1e-30 && ((x) - (y)) > -1e-30) // Сравнение double
#define OUT(x) #x << " = " << (x)


// ***************************
// Константы задачи:
#define ALFA TO_RAD(46) // Предельный рабочий угол по большой оси в радианах
#define SUN_CONSTANT 1373 // Вт/м^2
#define EARTH_RAD 6.4e6 // Метры
#define ORBIT_RAD 6.8e6 // Метры
#define MID_ALBEDO 0.6 // Среднее альбедо Земли. Безразмерная

#define ANGLE_PARTITION 1000
#define DIOD_PARTITION 6
// ***************************


struct Direction {
    double alfa;
    double gamma;
};

class Vec {
public:
    Vec(double X, double Y, double Z) : x_(X), y_(Y), z_(Z) {
        if(IS_EQL(X, 0))
            x_ = 0;
        if(IS_EQL(Y, 0))
            y_ = 0;
        if(IS_EQL(Z, 0))
            z_ = 0;
        r_ = sqrt(SQR(X) + SQR(Y) + SQR(Z));
    }

    // Нормировка на R
    void Norm(double R) {
        if(R < 0)
            throw invalid_argument("Нормировочная константа < 0!");
        if(IS_EQL(r_, 0) && !IS_EQL(R, 0)) {
            throw invalid_argument("Попытка задать ненулевую длину нулевому вектору!");
        } else if(!IS_EQL(r_, 0)) {
            x_ = R * x_ / r_;
            y_ = R * y_ / r_;
            z_ = R * z_ / r_;
            r_ = R;
        }
    }

    Vec operator-() {
        return Vec(-x_, -y_, -z_);
    }

    double x() const {return x_;}
    double y() const {return y_;}
    double z() const {return z_;}
    double r() const {return r_;}

private:
    double x_ = 0;
    double y_ = 0;
    double z_ = 0;
    double r_ = 0;
};

Vec Norm(const Vec& v, double R) {
    Vec copy = v;
    copy.Norm(R);
    return copy;
}

ostream& operator<<(ostream& os, const Vec& v) {
    return os << "(" << v.x() << ", " << v.y() << ", " << v.z() << ") || r = " << v.r();
}

double scalar(const Vec& lhs, const Vec& rhs) {
    return lhs.x() * rhs.x() + lhs.y() * rhs.y() + lhs.z() * rhs.z();
}

// Условие стабильной работы: либо косинус обоих углов больше 0, либо меньше!
// Особенность СК
Vec dir_to_vec(double alfa, double gamma) {
    //cout << "dir_to_vec(" << TO_GRAD(alfa) << ", " << TO_GRAD(gamma) << ") started:   ";

    if(cos(alfa) * cos(gamma) < 0) {
        throw invalid_argument("alfa и gamma указывают в разные полупространства!");
    }
    if(IS_EQL(cos(alfa) * cos(gamma), 0)) {
        throw invalid_argument("alfa и gamma должны указывать либо в верхнее, либо в нижнее полупространства!");
    }

    double z = (cos(alfa) > 0 && cos(gamma) > 0) ? 1 : -1;
    //cout << "z = " << z << ", vec = " << Norm(Vec(tan(alfa) * z, tan(gamma) * z, z), 1) << endl;

    return Norm(Vec(tan(alfa) * z, tan(gamma) * z, z), 1);
}
