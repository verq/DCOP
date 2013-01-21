#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include "problem.hpp"

int number_of_problem = 4;
const int number_of_possibly_feasible = 5;
const int number_of_feasible = 7;
int constraint_x1 = 3; //TODO
int constraint_x2 = 4; //TODO
double crossover_rate = 0.3;
double mutation_rate = 0.3;
int number_of_generations = 101;

struct individual {
    double x1;
    double x2;
	double fitness;
	int position;
};

struct population {
    individual possibly_feasible[number_of_possibly_feasible];
	individual feasible[number_of_feasible];
};

population individuals;

void initialize_possibly_feasible();
void initizalie_feasible();

void search();
void change_generation();

void repair(individual ind);
individual selection_feasible(int ranking);
individual selection_possibly_feasible(int ranking);
individual crossover(individual parent1, individual parent2);
individual mutation(individual parent);

bool is_feasible(individual ind);
double evaluation(individual ind);
std::pair<double, double> constraints(individual ind, double t);

double randomize_small_value();

void write_population();
void write_individual(individual ind);


int main() {
    srand(time(NULL));
	initialize_possibly_feasible();
	initizalie_feasible();
	for (int generation; generation < number_of_generations; generation++) {
		for (int i = 0; i < number_of_possibly_feasible; i++) {
			search();
		}
		if (generation % 100 == 0) {
			for (int i = 0; i < number_of_feasible; i++) {
				change_generation();
			}
		}
	}
	write_population();
    return 0;
}

void initialize_possibly_feasible() {
    for (int i = 0; i < number_of_possibly_feasible; i++) {
		individuals.possibly_feasible[i].x1 = rand() % constraint_x1;
		individuals.possibly_feasible[i].x2 = rand() % constraint_x2;
		individuals.possibly_feasible[i].fitness = evaluation(individuals.possibly_feasible[i]);
	}    
}

void initizalie_feasible() {
    for (int i = 0; i < number_of_feasible; i++) {
		individual ind;
		ind.x1 = rand() % constraint_x1;
		ind.x2 = rand() % constraint_x2;
		ind.fitness = evaluation(ind);
		while (!is_feasible(ind)) {
			ind.x1 = rand() % constraint_x1;
			ind.x2 = rand() % constraint_x2;
			ind.fitness = evaluation(ind);
			write_individual(ind);
		}
		individuals.feasible[i] = ind;
	}
}

void search() {
	double p1 = randomize_small_value();
	double p2 = randomize_small_value();
	individual child;
	if (p1 < crossover_rate) {
		individual parent1 = selection_possibly_feasible(number_of_possibly_feasible - 1);
		individual parent2 = selection_possibly_feasible(number_of_possibly_feasible - 2);
		child = crossover(parent1, parent2);
	}
	if (p2 < mutation_rate) {
		individual parent = selection_possibly_feasible(number_of_possibly_feasible - 1);
		child = mutation(parent);
	}
	if (p1 >= crossover_rate && p2 >= mutation_rate) {
		child = selection_possibly_feasible(number_of_possibly_feasible - 1);
	}
	child.fitness = evaluation(child);
	repair(child);
	individual worst = selection_possibly_feasible(0);
	worst = child;
}

void change_generation() {
	for (int i = 0; i < number_of_feasible; i++) {
		if (randomize_small_value() < crossover_rate) {
			individual parent1 = selection_feasible(number_of_feasible - 1);
			individual parent2 = selection_feasible(number_of_feasible - 2);
			individual child = crossover(parent1, parent2);
			if (is_feasible(child)) {
				child.fitness = evaluation(child);
				parent1.fitness = evaluation(parent1);
				if (child.fitness > parent1.fitness) {
					parent1 = child;
					parent1.fitness = child.fitness;
				}
			}
		}
		if (randomize_small_value() < mutation_rate) {
			individual parent = selection_feasible(number_of_feasible - 1);
			individual child = mutation(parent);
			if (is_feasible(child)) {
				child.fitness = evaluation(child);
				parent.fitness = evaluation(parent);
					if (child.fitness > parent.fitness) {
					parent = child;
					parent.fitness = child.fitness;
				}
			}
		}
	}
}

void repair(individual ind) {
	//TODO
}
individual selection_feasible(int ranking) {
	//TODO
}
individual selection_possibly_feasible(int ranking) {
	//TODO
}
individual crossover(individual parent1, individual parent2) {
	//TODO
}
individual mutation(individual parent) {
	//TODO
}

double evaluation(individual ind) {
	double t = ((int) time(NULL)) % 100;
	return dynamic_p1(number_of_problem, t) * (ind.x1 + dynamic_q1(number_of_problem, t)) +
			dynamic_p2(number_of_problem, t) * (ind.x2 + dynamic_q2(number_of_problem, t));
} 

bool is_feasible(individual ind) {
	double t = ((int) time(NULL)) % 100;
	std::pair<double, double> pair_of_constraints = constraints(ind, t);
	return (pair_of_constraints.first <= 0 && pair_of_constraints.second <= 0);
}


std::pair<double, double> constraints(individual ind, double t) {	
	double y1 = dynamic_r1(number_of_problem, t) * (ind.x1 + dynamic_s1(number_of_problem, t));
	double y2 = dynamic_r2(number_of_problem, t) * (ind.x2 + dynamic_s2(number_of_problem, t));
	printf("t %lf\n", t);
	printf("r1 %lf \t r2 %lf\n", dynamic_r1(number_of_problem, t), dynamic_r1(number_of_problem, t));
	printf("s1 %lf \t s2 %lf\n", dynamic_s1(number_of_problem, t), dynamic_s2(number_of_problem, t));
	printf("y1 %lf \t y2 %lf \t g_1 %lf \t g_2 %lf \n", y1, y2, constraint_1(y1, y2), constraint_2(y1, y2));
	switch (number_of_problem) {
		case 4 : return std::make_pair(constraint_1(y1, y2), constraint_2(y1, y2));
	}
	return std::make_pair(0, 0);
}

double randomize_small_value() {
	return 1.0 / (rand() % 100);
}

void write_population() {
	printf("Population \n\tPossibly feasible: \n");
	for (int i = 0; i < number_of_possibly_feasible; i++) {
		printf("\t%lf %lf fitness %lf\n", individuals.possibly_feasible[i].x1, individuals.possibly_feasible[i].x2, individuals.possibly_feasible[i].fitness);
	}
	printf("Feasible\n");
	for (int i = 0; i < number_of_feasible; i++) {
		printf("\t%lf %lf fitness %lf\n", individuals.feasible[i].x1, individuals.feasible[i].x2, individuals.feasible[i].fitness);	
	}
	printf("End of population \n \n");
}

void write_individual(individual ind) {
	printf("%lf %lf fitness %lf\n", ind.x1, ind.x2, ind.fitness);
}
