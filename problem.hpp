#ifndef problem_hpp
#define problem_hpp

double constraint_1();
double constraint_2();

double dynamic_p1(int problem, double t);
double dynamic_p2(int problem, double t);

double dynamic_q1(int problem, double t);
double dynamic_q2(int problem, double t);

double dynamic_r1(int problem, double t);
double dynamic_r2(int problem, double t);

double dynamic_s1(int problem, double t);
double dynamic_s2(int problem, double t);

double constraint_1(double y1, double y2);
double constraint_2(double y1, double y2);

#endif