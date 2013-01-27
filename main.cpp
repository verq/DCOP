#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <algorithm>
#include "problem.hpp"

//from paper
int constraint_x1 = 3;
int constraint_x2 = 4;
double crossover_rate = 0.1;
double mutation_rate = 0.15;
double replacement_rate = 0.0;              //in paper: table VI, p.5 (Baldwinian evolution)
int number_of_changes = 10;                 //in paper: table VI, p.5
int change_frequency = 1000;                //in paper: table VI, p.5
const int number_of_possibly_feasible = 20; //in paper: table VI, p.5
const int number_of_feasible = 5;           //in paper: table VI, p.5

//not from paper
int number_of_problem = 4;
int number_of_generations = 1000;

int function_evaluations_counter = 0;       // number of function evaluations so far

struct individual {
    double x1;
    double x2;
	double fitness;
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

individual repair(individual ind);
individual* selection_feasible(int ranking);
individual* selection_possibly_feasible(int ranking);
individual arithmetic_crossover(individual parent1, individual parent2);
individual uniform_crossover(individual parent1, individual parent2);
individual mutation(individual parent);

bool is_feasible(individual ind);
double evaluation(individual ind);
std::pair<double, double> constraints(individual ind, double t);
bool compare(individual i,individual j);

double randomize_small_value(int);
double count_time();

void write_population();
void write_individual(individual ind);


int main() {
    srand(time(NULL));
	initialize_possibly_feasible();
	initizalie_feasible();
	function_evaluations_counter = 0;
	// function should change every <change_frequency> function evaluations
	// there should be <number_of_changes> changes of function
	// hence:
	while(function_evaluations_counter < number_of_changes * change_frequency)
    {
		for (int i = 0; i < number_of_possibly_feasible; i++) {
			search();
		}
		if (function_evaluations_counter % 100 == 0) {
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
		individuals.possibly_feasible[i].x1 = randomize_small_value(constraint_x1);
		individuals.possibly_feasible[i].x2 = randomize_small_value(constraint_x2);
		individuals.possibly_feasible[i].fitness = evaluation(individuals.possibly_feasible[i]);
	}
}

void initizalie_feasible() {
    for (int i = 0; i < number_of_feasible; i++) {
		individual ind;
		ind.x1 = randomize_small_value(constraint_x1);
		ind.x2 = randomize_small_value(constraint_x2);
		ind.fitness = evaluation(ind);
		while (!is_feasible(ind)) {
			ind.x1 = randomize_small_value(constraint_x1);
			ind.x2 = randomize_small_value(constraint_x2);
			ind.fitness = evaluation(ind);
		}
		individuals.feasible[i] = ind;
	}
}

void search() {
	double p1 = randomize_small_value(1);
	double p2 = randomize_small_value(1);
	individual child;
	if (p1 < crossover_rate) {
	    //TODO: nonlinear ranking + should choose one of the worst, not the worst
		individual parent1 = *selection_possibly_feasible(number_of_possibly_feasible - 1);
		individual parent2 = *selection_possibly_feasible(number_of_possibly_feasible - 2);  //TODO: nonlinear ranking
		child = uniform_crossover(parent1, parent2);
        child.fitness = evaluation(child);
        child = repair(child);
        individual *worst = selection_possibly_feasible(0);  //TODO: replace ONE OF the worst, not the worst + nonlinear ranking
        *worst = child;      // does it work?
	}
	if (p2 < mutation_rate) {
		individual parent = *selection_possibly_feasible(number_of_possibly_feasible - 1);   //TODO: nonlinear ranking
		child = mutation(parent);
        child.fitness = evaluation(child);
        ichild = repair(child);
        individual *worst = selection_possibly_feasible(0);  //TODO: replace ONE OF the worst, not the worst + nonlinear ranking
        *worst = child;      // does it work?
	}
	if (p1 >= crossover_rate && p2 >= mutation_rate) {
		child = *selection_possibly_feasible(number_of_possibly_feasible - 1);   //TODO: nonlinear ranking
		//TODO: check if child hasn't been evaluated since the last generation
		child.fitness = evaluation(child);
        child = repair(child);
        individual *worst = selection_possibly_feasible(0);  //TODO: replace ONE OF the worst, not the worst + nonlinear ranking
        *worst = child; // does it work?
	}
}

void change_generation() {
	//for (int i = 0; i < number_of_feasible; i++) {    this loop is already in main function
    if (randomize_small_value(1) < crossover_rate) {
        individual parent1 = *selection_feasible(number_of_feasible - 1);    //TODO: nonlinear ranking
        individual parent2 = *selection_feasible(number_of_feasible - 2);    //TODO: nonlinear ranking
        individual child = uniform_crossover(parent1, parent2);
        if (is_feasible(child)) {
            child.fitness = evaluation(child);
            //from paper: "evaluate r and x, the better of the two parents"
            parent1.fitness = evaluation(parent1);  //don't know if it's necessary. maybe fitness value is up to date
            parent2.fitness = evaluation(parent2);
            individual *better_parent;
            if(parent1.fitness > parent2.fitness) {
                better_parent = &parent1;
            }
            else {
                better_parent = &parent2;
            }
            //parent1.fitness = evaluation(parent1);
            if (child.fitness > (*better_parent).fitness) {
                *better_parent = child;
                //better_parent.fitness = child.fitness;   not necessary
            }
        }
    }
    if (randomize_small_value(1) < mutation_rate) {
        individual parent = *selection_feasible(number_of_feasible - 1);
        individual child = mutation(parent);
        if (is_feasible(child)) {
            child.fitness = evaluation(child);
            parent.fitness = evaluation(parent);
                if (child.fitness > parent.fitness) {
                parent = child;
                //parent.fitness = child.fitness;   //not necessary
            }
        }
    }
}

individual repair(individual ind) {
	int random_feasible_index = rand() % number_of_feasible;
	individual repaired;
	int iter = 0; //number_of_trials = 100 in paper
	do {
		double p = randomize_small_value(1);
		repaired.x1 = p * ind.x1 + (1 - p) * individuals.feasible[random_feasible_index].x1; //linear interpolation
		repaired.x2 = p * ind.x2 + (1 - p) * individuals.feasible[random_feasible_index].x2; //linear interpolation
		iter++;
	} while (!is_feasible(repaired) && iter != 100);
	if (iter == 100) {
		repaired = individuals.feasible[random_feasible_index];
	}
	repaired.fitness = evaluation(repaired);
	if (repaired.fitness > individuals.feasible[random_feasible_index].fitness) {
		individuals.feasible[random_feasible_index] = repaired;
	}
	ind.fitness = repaired.fitness;
	return ind;
}

individual* selection_feasible(int ranking) { //linear selection unfortunetly
	std::sort(individuals.feasible, individuals.feasible + number_of_feasible, compare);
	return &individuals.feasible[ranking];
}
individual* selection_possibly_feasible(int ranking) { //linear selection unfortunetly
	std::sort(individuals.feasible, individuals.possibly_feasible + number_of_possibly_feasible, compare);
	return &individuals.possibly_feasible[ranking];
}
individual arithmetic_crossover(individual parent1, individual parent2) {
	individual child;
	child.x1 = parent1.x1 + parent2.x1;
	child.x2 = parent1.x2 + parent2.x2;
	child.fitness = evaluation(child);
	return child;
}

individual uniform_crossover(individual parent1, individual parent2) { //in paper
	individual child;

	if(rand() % 2 == 0)
        child.x1 = parent1.x1;
    else
        child.x1 = parent2.x1;

    if(rand() % 2 == 0)
        child.x2 = parent1.x2;
    else
        child.x2 = parent2.x2;

	return child;
}

individual mutation(individual parent) {
	individual child;
	child.x1 = parent.x1;
	child.x2 = parent.x2;
	//int random_gen = rand() % 2;

	int small_value1 = randomize_small_value(1);
	int small_value2 = randomize_small_value(1);

	if(rand() % 2 == 0)
        child.x1 += small_value1;
    else
        child.x1 -= small_value1;

    if(rand() % 2 == 0)
        child.x2 += small_value2;
    else
        child.x2 -= small_value2;

	//if (random_gen == 0) child.x1 = randomize_small_value() * constraint_x1;
	//else child.x2 = randomize_small_value() * constraint_x2;

	return child;
}

double evaluation(individual ind) {
    function_evaluations_counter++;     // should be increased with each function evaluation
	int t = function_evaluations_counter / change_frequency;    // t is increased every <change_frequency> function evaluations
	return dynamic_p1(number_of_problem, t) * (ind.x1 + dynamic_q1(number_of_problem, t)) +
			dynamic_p2(number_of_problem, t) * (ind.x2 + dynamic_q2(number_of_problem, t));
}

bool is_feasible(individual ind) {
	//double t = count_time();
	int t = function_evaluations_counter / change_frequency;
	std::pair<double, double> pair_of_constraints = constraints(ind, t);
	return (pair_of_constraints.first <= 0 && pair_of_constraints.second <= 0);
}


std::pair<double, double> constraints(individual ind, double t) {
	double y1 = dynamic_r1(number_of_problem, t) * (ind.x1 + dynamic_s1(number_of_problem, t));
	double y2 = dynamic_r2(number_of_problem, t) * (ind.x2 + dynamic_s2(number_of_problem, t));
	switch (number_of_problem) { //some of the problems don't have contraints, some have only one
		case 4 : return std::make_pair(constraint_1(y1, y2), constraint_2(y1, y2));
	}
	return std::make_pair(0, 0);
}

double randomize_small_value(int not_bigger_than) {
	//return 1.0 / (rand() % 100);    could be divided by zero
	return (rand() % (not_bigger_than * 100)) * 0.01;
}

//double count_time() {
//	return ((int) time(NULL)) % 100;
//}

bool compare(individual i,individual j) {
	return (i.fitness < j.fitness);
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
