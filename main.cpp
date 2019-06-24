#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <set>

#include "comunity.h"
#include "sensor.h"
#include "earth.h"
#include "sun.h"

using namespace std;

/*
Расматривается датчик на квадрантном фотодиоде 3х2
Рабочий диапазон датчика по углам: 55x92
*/


Direction evaluate_direction(Sensor sensor) {
    double min = MIN(sensor.Q11, MIN(sensor.Q12, MIN(sensor.Q13, MIN(sensor.Q21, MIN(sensor.Q22, sensor.Q23)))));
    double alfa, gamma;

    sensor.Q21 -= min; sensor.Q22 -= min; sensor.Q23 -= min;
    sensor.Q11 -= min; sensor.Q12 -= min; sensor.Q13 -= min;

    if((sensor.Q11 + sensor.Q21) >= (sensor.Q13 + sensor.Q23)) {
        alfa = atan(sensor.param_k_a * (sensor.Q11 + sensor.Q21)
                    / (sensor.Q11 + sensor.Q12 + sensor.Q21 + sensor.Q22));
        gamma = atan(sensor.param_k_g * (sensor.Q11 + sensor.Q12 - sensor.Q21 - sensor.Q22)
                     / (sensor.Q11 + sensor.Q12 + sensor.Q21 + sensor.Q22));
    } else {
        alfa = atan(-sensor.param_k_a * (sensor.Q13 + sensor.Q23)
                    / (sensor.Q12 + sensor.Q13 + sensor.Q22 + sensor.Q23));
        gamma = atan(sensor.param_k_g * (sensor.Q12 + sensor.Q13 - sensor.Q22 - sensor.Q23)
                     / (sensor.Q12 + sensor.Q13 + sensor.Q22 + sensor.Q23));
    }

    return {alfa, gamma};
}


int main()
{
    // Задаем углы на Солнце и на край Земли
    Sun sun(TO_RAD(-23.5), TO_RAD(0));
    Earth earth(TO_RAD(23.5), TO_RAD(0));
    Sensor sensor;

    cout << "Солнце:" << endl
         << "alfa = " << TO_GRAD(sun.alfa) << endl
         << "gamma = " << TO_GRAD(sun.gamma) << endl
         << "вектор на Солнце: " << sun.direction << endl << endl;

    cout << "Земля:" << endl
         << "направление на край:" << endl
         << "border_alfa = " << TO_GRAD(earth.border_alfa) << endl
         << "border_gamma = " << TO_GRAD(earth.border_gamma) << endl
         << "направление на центр" << endl
         << "alfa = " << TO_GRAD(earth.alfa) << endl
         << "gamma = " << TO_GRAD(earth.gamma) << endl
         << "вектор в центр Земли" << earth.position << endl << endl;

    cout << "alfa_range = " << TO_GRAD(sensor.alfa_range) << endl
         << "gamma_range = " << TO_GRAD(sensor.gamma_range) << endl << endl;

    /*
    Небольшой ликбез:

    Уравнение поверхности Земли
    (x - earth.position.x())^2 + (y - earth.position.y())^2 + (z - earth.position.z())^2 = earth.rad^2

    Направление, заданное двумя углами alfa и gamma, преобразуем к виду двух пересекающихся плоскостей или зададим параметрически:

    |z = t
    |x = z * tan(alfa) = t * tan(alfa)
    |y = z * tan(gamma) = t * tan(gamma)
    |t > 0 Нижняя полуплоскость нас не интересует

    Подставляем в Ур-е пов-ти Земли:
    (t * tan(alfa) - earth.position.x())^2 + (t * tan(gamma) - earth.position.y())^2 + (t - earth.position.z())^2 = earth.rad^2 (*)

    Если есть решения t > 0, то в направлении (alfa, gamma) видна Земля. t дает точку пересечения.
    */


    /*
        Описание режимов работы:
        mod == 1: По каждому направлению выводится dI
        mod == 2: Иллюстрация по углам
        mod == 3:
        mod == 4: Вывод таблицы интенсивностей на разных участках диодов в стандартный поток вывода
    */
    set<int> mod = {4};

    double image[3 * DIOD_PARTITION][2 * DIOD_PARTITION] = {{0}};

    for(int a = 0; a < ANGLE_PARTITION; a++) {
        for(int g = 0; g < ANGLE_PARTITION; g++) {

            double current_alfa = -sensor.beta_range_a + (a + 0.5) * (2 * sensor.beta_range_a / ANGLE_PARTITION);
            double current_gamma = -sensor.beta_range_g + (g + 0.5) * (2 * sensor.beta_range_g / ANGLE_PARTITION);
            Vec current_direction = Norm(Vec(tan(current_alfa), tan(current_gamma), 1), 1);

            Vec intersection(0, 0, 0);
            Vec n(0, 0, 0);
            double cos_theta;
            double cos_phy;
            double cos_theta_sensor = scalar(current_direction, Vec(0, 0, 1)) / current_direction.r();
            double current_intensity_from_earth = 0;

            // Ищем дискриминант (*)
            double a = SQR(tan(current_alfa)) + SQR(tan(current_gamma)) + 1;
            double b = -2 * (earth.position.x() * tan(current_alfa) + earth.position.y() * tan(current_gamma) + earth.position.z());
            double c = SQR(ORBIT_RAD) - SQR(EARTH_RAD);

            double D = SQR(b) - 4 * a * c;

            if(D > 0) {
                double t = (-b - sqrt(D)) / (2 * a);
                if(t > 0) {
                    // Точка пересечения луча направления (current_alfa, current_gamma) и поверхности Земли
                    intersection = Vec(t * tan(current_alfa), t * tan(current_gamma), t);

                    // Нормаль (n.x(), n.y(), n.z()) к поверхности Земли в найденной точке (intersection.x(), intersection.y(), intersection.z())
                    n = Norm(Vec(intersection.x() - earth.position.x()
                                 , intersection.y() - earth.position.y()
                                 , intersection.z() - earth.position.z()), 1);

                    // Косинус падения солнечных лучей на поверхность Земли в точке intersection
                    cos_theta = scalar(sun.direction, n) / (sun.direction.r() * n.r());

                    // Косинус угла между нормалью к поверхности и направлением на датчик
                    cos_phy = scalar(-current_direction, n) / (current_direction.r() * n.r());

                    double d_omega = (2 * sensor.beta_range_a / ANGLE_PARTITION) * (2 * sensor.beta_range_g / ANGLE_PARTITION);
                    current_intensity_from_earth = SUN_CONSTANT * MID_ALBEDO * cos_theta * d_omega / M_PI;
                    if(current_intensity_from_earth < 0)
                        current_intensity_from_earth = 0;
                }
            }


            if(mod.count(1))
                cout << "a = " << TO_GRAD(current_alfa) << "; g = " << TO_GRAD(current_gamma) << "dI = " << current_intensity_from_earth << endl;


            // Вычислили интенсивность излучения с заданного направления, теперь наращиваем сигнал на приемнике
            if(!IS_EQL(current_intensity_from_earth, 0)) {

                double shift_x_ = -tan(current_alfa) / (sensor.param_k_a); // В единицах L [-1.5, 1.5]
                double left_x_ = max(-0.5 + shift_x_, -1.5);
                double right_x_ = min(0.5 + shift_x_, 1.5);
                int left_x = 3 * DIOD_PARTITION / 2 + static_cast<int>(round(left_x_ * DIOD_PARTITION)); // В дискретах
                int right_x = 3 * DIOD_PARTITION / 2 + static_cast<int>(round(right_x_ * DIOD_PARTITION)) - 1;

                double shift_y_ = -tan(current_gamma) / (2 * sensor.param_k_g); // В единицах L [-1.5, 1.5]
                double left_y_ = max(-0.5 + shift_y_, -1.);
                double right_y_ = min(0.5 + shift_y_, 1.);
                int left_y = DIOD_PARTITION + static_cast<int>(round(left_y_ * DIOD_PARTITION)); // В дискретах
                int right_y = DIOD_PARTITION + static_cast<int>(round(right_y_ * DIOD_PARTITION)) - 1;


                for(int current_x = left_x; current_x <= right_x; current_x++) {
                    for(int current_y = left_y; current_y <= right_y; current_y++) {
                        image[current_x][current_y] += current_intensity_from_earth * cos_theta_sensor;
                    }
                }
            }

            // Иллюстрация по углам
            if(mod.count(2)) {
                if(!IS_EQL(intersection.r(), 0))
                    cout << "+ ";
                else
                    cout << "- ";
            }

        }

        // Иллюстрация по углам
        if(mod.count(2))
            cout << endl;
    }


    // Вывод таблицы интенсивностей на разных участках диодов в стандартный поток вывода
    cout << "Распределение освещенности на фотодиодах:" << endl;
    if(mod.count(4)) {
        cout << endl;
        for(int iy = 2 * DIOD_PARTITION - 1; iy >=0; iy--) {
            cout << " ";
            for(int ix = 0; ix < 3 * DIOD_PARTITION; ix++) {
                cout << setw(5) << setprecision(4) << image[ix][iy] << "  ";
            }
            cout << endl;
        }
    }


    cout << setprecision(6) << endl
         << "Реальное направление на Солнце (плюс-минус результат, который выдаёт датчик без засветки от Земли):" << endl
         << OUT(TO_GRAD(sun.alfa)) << endl
         << OUT(TO_GRAD(sun.gamma)) << endl << endl;


    double cos_sun_direction = scalar(sun.direction, Vec(0, 0, 1));

    // В единицах L
    double min_x_board = max(-1.5, -0.5 - tan(sun.alfa) / (sensor.param_k_a));
    double max_x_board = min(1.5, 0.5 - tan(sun.alfa) / (sensor.param_k_a));
    double min_y_board = max(-1., -0.5 - tan(sun.gamma) / (2 * sensor.param_k_g));
    double max_y_board = min(1., 0.5 - tan(sun.gamma) / (2 * sensor.param_k_g));

    cout << OUT(min_x_board) << "  " << OUT(max_x_board) << endl
         << OUT(min_y_board) << "  " << OUT(max_y_board) << endl;

    sensor.Q21 = SUN_CONSTANT * cos_sun_direction
            * (min(1., max_y_board) - max(0., min_y_board))
            * (min(-0.5, max_x_board) - max(-1.5, min_x_board));
    if(sensor.Q21 < 0)
        sensor.Q21 = 0;
    sensor.Q22 = SUN_CONSTANT * cos_sun_direction
            * (min(1., max_y_board) - max(0., min_y_board))
            * (min(0.5, max_x_board) - max(-0.5, min_x_board));
    if(sensor.Q22 < 0)
        sensor.Q22 = 0;
    sensor.Q23 = SUN_CONSTANT * cos_sun_direction
            * (min(1., max_y_board) - max(0., min_y_board))
            * (min(1.5, max_x_board) - max(0.5, min_x_board));
    if(sensor.Q23 < 0)
        sensor.Q23 = 0;
    sensor.Q11 = SUN_CONSTANT * cos_sun_direction
            * (min(0., max_y_board) - max(-1., min_y_board))
            * (min(-0.5, max_x_board) - max(-1.5, min_x_board));
    if(sensor.Q11 < 0)
        sensor.Q11 = 0;
    sensor.Q12 = SUN_CONSTANT * cos_sun_direction
            * (min(0., max_y_board) - max(-1., min_y_board))
            * (min(0.5, max_x_board) - max(-0.5, min_x_board));
    if(sensor.Q12 < 0)
        sensor.Q12 = 0;
    sensor.Q13 = SUN_CONSTANT * cos_sun_direction
            * (min(0., max_y_board) - max(-1., min_y_board))
            * (min(1.5, max_x_board) - max(0.5, min_x_board));
    if(sensor.Q13 < 0)
        sensor.Q13 = 0;

    cout << OUT(sensor.Q21) << "  " << OUT(sensor.Q22) << "  " << OUT(sensor.Q23) << endl
         << OUT(sensor.Q11) << "  " << OUT(sensor.Q12) << "  " << OUT(sensor.Q13) << endl << endl;

    Direction dir = evaluate_direction(sensor);

    cout << "Расчетное направление на Солнце (без засветки от Земли):" << endl
         << "alfa == " << TO_GRAD(dir.alfa) << endl
         << "gamma == " << TO_GRAD(dir.gamma) << endl << endl;




    double Q11_E = 0, Q12_E = 0, Q13_E = 0, Q21_E = 0, Q22_E = 0, Q23_E = 0;
    for(int ix = 0; ix < DIOD_PARTITION; ix++)
        for(int iy = 0; iy < DIOD_PARTITION; iy++) {
            Q11_E += image[ix][iy];
            Q12_E += image[DIOD_PARTITION + ix][iy];
            Q13_E += image[2 * DIOD_PARTITION + ix][iy];
            Q21_E += image[ix][DIOD_PARTITION + iy];
            Q22_E += image[DIOD_PARTITION + ix][DIOD_PARTITION + iy];
            Q23_E += image[2 * DIOD_PARTITION + ix][DIOD_PARTITION + iy];
        }
    Q11_E /= SQR(DIOD_PARTITION);
    Q12_E /= SQR(DIOD_PARTITION);
    Q13_E /= SQR(DIOD_PARTITION);
    Q21_E /= SQR(DIOD_PARTITION);
    Q22_E /= SQR(DIOD_PARTITION);
    Q23_E /= SQR(DIOD_PARTITION);

    cout << OUT(Q11_E) << "  " << OUT(Q12_E) << "  " << OUT(Q13_E) << endl
         << OUT(Q21_E) << "  " << OUT(Q22_E) << "  " << OUT(Q23_E) << endl << endl;

    sensor.Q11 += Q11_E;
    sensor.Q12 += Q12_E;
    sensor.Q13 += Q13_E;
    sensor.Q21 += Q21_E;
    sensor.Q22 += Q22_E;
    sensor.Q23 += Q23_E;

    dir = evaluate_direction(sensor);

    cout << "Расчетное направление на Солнце (c засветкой от Земли):" << endl
         << "alfa == " << TO_GRAD(dir.alfa) << endl
         << "gamma == " << TO_GRAD(dir.gamma) << endl << endl;

    return 0;
}
