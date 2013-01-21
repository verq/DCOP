#include <cmath>


double severity_of_function_changes = 0.5;
double severity_of_constraint_changes = 20;


double dynamic_p1(int problem, double t) {
    switch(problem) {
        case 4: return sin(severity_of_function_changes * M_PI * t + M_PI / 2);
    }
    return 0;
}

double dynamic_p2(int problem, double t) {
    switch(problem) {
        case 4: return 1;
    }
    return 0;
}

double dynamic_q1(int problem, double t) {
    switch(problem) {
        case 4: return 0;
    }
    return 0;  
}

double dynamic_q2(int problem, double t) {
    switch(problem) {
        case 4: return 0;
    }
    return 0;
}

double dynamic_r1(int problem, double t) {
    switch(problem) {
        case 4: return 1;
    }
    return 0;  
}

double dynamic_r2(int problem, double t) {
    switch(problem) {
        case 4: return 1;
    }
    return 0;
}

double dynamic_s1(int problem, double t) {
    switch(problem) {
        case 4: return 0;
    }
    return 0;
}

double dynamic_s2(int problem, double t) {
    switch(problem) {
        case 4: return t * 4.0/severity_of_constraint_changes;
    }
	return 0;
}

double constraint_1(double y1, double y2) {
	return -2 * pow(y1, 4) + 8 * pow(y1, 3) - 8 * pow(y1, 2) + y2 - 2;
}

double constraint_2(double y1, double y2) {
	return -4 * pow(y1, 4) + 32 * pow(y1, 3) - 88 * pow(y1, 2) + 96 * y1 + y2 - 36;
}