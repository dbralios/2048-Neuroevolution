#ifndef GAME_H
#define GAME_H
#include "matrix.h"
#include "neural.h"
#include <vector>
#include <utility>
#include <random>
#include <iostream>

class game{
public:
	game();
	game(unsigned int n);
	friend int play(Net, unsigned int n);
	friend int play(Net, ostream &out, unsigned int n);
private:
	bool Move(int id);
	bool Up();
	bool Down();
	bool Left();
	bool Right();

	void random();
	void restart();
	bool isOver();
	int  getScore();
	void printGrid(ostream &out);
	matrix<int> getGrid();

	void availableCells(std::vector<std::pair<int, int>> &cells);
	bool transform(std::vector<int> &arr);
	bool isValid(std::vector<int> &arr);

	int score;
	const unsigned int N;
	matrix<int> grid;
};
#endif
