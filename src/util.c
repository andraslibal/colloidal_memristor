#include <math.h>
#include "util.h"

#define PI 3.14159265358979323846

int get_sign(double number){
    if (signbit(number)) {
        return -1;
    } else {
        return 1;
    }
}

double to_first_quarter(double angle_rad){
    if (angle_rad < 0)
        return 2*PI + angle_rad;
    if (angle_rad >= 0)
        return angle_rad;
    if (angle_rad > PI/2)
        return PI - angle_rad;
    if (angle_rad > PI)
        return angle_rad - PI;
    if (angle_rad > 3*PI/2)
        return 2*PI - angle_rad;
    if (angle_rad > 2*PI)
        return angle_rad - 2*PI;
}