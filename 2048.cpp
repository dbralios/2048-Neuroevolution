#include "2048.h"

game::game() : game(4) { }

game::game(unsigned int n) : score(0), N(n)
{
	grid = matrix<int>(n, n);
	for (auto it = grid.begin(); it != grid.end(); it++){
		*it = 0;
	}
}

int game::getScore()
{
	return score;
}

void game::restart() 
{
	score = 0;
	for (auto it = grid.begin(); it != grid.end(); it++){
		*it = 0;
	}
}

void game::random()
{
	std::vector<std::pair<int, int>> cells;
	availableCells(cells);

	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> dist_int(0, cells.size() - 1), dist(1, 10);

	grid(cells[dist_int(mt)].first, cells[dist_int(mt)].second) = dist(mt) > 1 ? 2 : 2; // 4;
}

void game::availableCells(std::vector<std::pair<int, int>> &cells)
{
	for (unsigned i = 0; i < grid.getRows(); ++i){
		for (unsigned j = 0; j < grid.getColumns(); ++j){
			if (grid(i, j) == 0){
				cells.push_back({ i, j });
			}
		}
	}
}

matrix<int> game::getGrid()
{
	return grid;
}

bool game::isOver()
{
	std::vector<std::pair<int, int>> cells;
	availableCells(cells);
	if(cells.size() != 0)
		return false;

	std::vector<int> arr;
	bool valid = false;
	for (unsigned i = 0; i < grid.getColumns(); ++i){
		arr.clear();
		for (unsigned j = 0; j < grid.getRows(); ++j){
			arr.push_back(grid(j, i));
		}
		valid ? isValid(arr) : valid = isValid(arr);
		std::reverse(arr.begin(), arr.end());
		valid ? isValid(arr) : valid = isValid(arr);
	}
	for (unsigned i = 0; i < grid.getRows(); ++i){
		arr.clear();
		for (unsigned j = 0; j < grid.getColumns(); ++j){
			arr.push_back(grid(i, j));
		}
		valid ? isValid(arr) : valid = isValid(arr);
		std::reverse(arr.begin(), arr.end());
		valid ? isValid(arr) : valid = isValid(arr);
	}
	return !valid;
}

bool game::isValid(std::vector<int> &arr)
{
	int value;
	bool valid = false;
	for (unsigned i = 0; i < arr.size() - 1; ++i) {
		value = arr[i];
		for (unsigned j = i + 1; j < arr.size(); ++j) {
			if (arr[j] != 0) {
				if (value == arr[j] || value == 0){
					arr[i] += arr[j];
					arr[j] = 0;
					if (!value) --i;
					valid = true;
				}
				break;
			}
		}
	}
	return valid;
}

bool game::transform(std::vector<int> &arr)
{
	/*int value;
	bool valid = false;
	for (unsigned i = 0; i < arr.size() - 1; ++i) {
		value = arr[i];
		for (unsigned j = i + 1; j < arr.size(); ++j) {
			if (arr[j] != 0) {
				if (value == arr[j] || value == 0){
					arr[i] += arr[j];
					arr[j] = 0;
					if (value == 0) {
						--i;
					}
					else {
						score += arr[i];
					}
					valid = true;
				}
				break;
			}
		}
	}
	return valid;*/

	int val = -1, pos = -1;
	bool valid = false;
	std::vector<int> A;
	for (unsigned i = 0; i < arr.size(); ++i) {		
		if (val == arr[i]) {
			valid = true;
			A.push_back(2 * val);
			score += 2 * val;
			val = -1;
		}
		else if (arr[i] != 0 && val != arr[i]) {
			if (pos + 1 < i) valid = true;
			if (val != -1) A.push_back(val);
			val = arr[i]; pos = i;
		}
	}
	if (val != -1) A.push_back(val);

	for (unsigned i = 0; i < arr.size(); ++i) {
		if (i >= A.size()) {
			arr[i] = 0;
		}
		else {
			arr[i] = A[i];
		}
	}
	return valid;
}

bool game::Move(int id)
{
	switch (id){
	case 0:
		return Up();
	case 1:
		return Down();
	case 2:
		return Left();
	case 3:
		return Right();
	default:
		return false;
	}
}

bool game::Up()
{
	std::vector<int> arr;
	bool valid = false;
	for (unsigned i = 0; i < grid.getColumns(); ++i){
		arr.clear();
		for (unsigned j = 0; j < grid.getRows(); ++j){
			arr.push_back(grid(j, i));
		}
		valid ?	transform(arr) : valid = transform(arr);
		for (unsigned j = 0; j < grid.getRows(); ++j){
			grid(j, i) = arr[j];
		}
	}
	return valid;
}

bool game::Down()
{
	std::vector<int> arr;
	bool valid = false;
	for (unsigned i = 0; i < grid.getColumns(); ++i){
		arr.clear();
		for (unsigned j = 0; j < grid.getRows(); ++j){
			arr.push_back(grid(j, i));
		}
		std::reverse(arr.begin(), arr.end());
		valid ? transform(arr) : valid = transform(arr);
		std::reverse(arr.begin(), arr.end());
		for (unsigned j = 0; j < grid.getRows(); ++j){
			grid(j, i) = arr[j];
		}
	}
	return valid;
}

bool game::Left()
{
	std::vector<int> arr;
	bool valid = false;
	for (unsigned i = 0; i < grid.getRows(); ++i){
		arr.clear();
		for (unsigned j = 0; j < grid.getColumns(); ++j){
			arr.push_back(grid(i, j));
		}
		valid ? transform(arr) : valid = transform(arr);
		for (unsigned j = 0; j < grid.getColumns(); ++j){
			grid(i, j) = arr[j];
		}
	}
	return valid;
}

bool game::Right()
{
	std::vector<int> arr;
	bool valid = false;
	for (unsigned i = 0; i < grid.getRows(); ++i){
		arr.clear();
		for (unsigned j = 0; j < grid.getColumns(); ++j){
			arr.push_back(grid(i, j));
		}
		std::reverse(arr.begin(), arr.end());
		valid ? transform(arr) : valid = transform(arr);
		std::reverse(arr.begin(), arr.end());
		for (unsigned j = 0; j < grid.getColumns(); ++j){
			grid(i, j) = arr[j];
		}
	}
	return valid;
}

void game::printGrid(ostream &out)
{
	for (unsigned i = 0; i < grid.getRows(); ++i) {
		for (unsigned j = 0; j < grid.getColumns(); ++j) {
			out << grid(i, j) << ' ';
		}
		out << std::endl;
	}
}

int play(Net network, unsigned int n)
{
	game gInst(n);
	const unsigned int N = gInst.N;

	gInst.random();
	gInst.random();

	Matrix input(N*N, 1), output(4, 1);
	matrix<int>  grid;
	while (true){  // !gInst.isOver()
		// Get input matrix
		grid = gInst.getGrid();
		int maxTile = 0;
		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
				if (grid(i, j) != 0){
					input(i * N + j, 0) = log2(grid(i, j));
					if (maxTile < input(i * N + j, 0))
						maxTile = (int)input(i * N + j, 0);
				}
				else {
					input(i * N + j, 0) = 0;
				}
			}
		}

		for (unsigned i = 0; i < input.getColumns(); i++){
			input(i, 0) /= 16;//maxTile;
		}

		// Feed forward
		network.feedForward(input, output);

		std::vector<std::pair<double, int>> decision;
		for (unsigned i = 0; i < output.getRows(); i++){
			decision.push_back({ output(i, 0), i });
		}
		// Sort output
		std::sort(decision.begin(), decision.end());
		// Moves
		if (gInst.Move(decision[3].second)){
			gInst.random();
		}
		else{
			break;
		}
	}
	return gInst.getScore();
}

int play(Net network, ostream &out, unsigned int n)
{
	std::vector<std::string> M;
	M.push_back("UP");
	M.push_back("DOWN");
	M.push_back("LEFT");
	M.push_back("RIGHT");

	game gInst(n);
	const unsigned int N = gInst.N;
	gInst.random();
	gInst.random();

	out << "\tGAME\t" << std::endl;
	gInst.printGrid(out);

	Matrix input(N*N, 1), output(4, 1);
	matrix<int>  grid;
	while (true){  // !gInst.isOver()
		// Get input matrix
		grid = gInst.getGrid();
		int maxTile = 0;
		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
				if (grid(i, j) != 0){
					input(i * N + j, 0) = log2(grid(i, j));
					if (maxTile < input(i * N + j, 0))
						maxTile = (int)input(i * N + j, 0);
				}
				else {
					input(i * N + j, 0) = 0;
				}
			}
		}

		for (unsigned i = 0; i < input.getColumns(); i++){
			input(i, 0) /= 16;//maxTile;
		}

		// Feed forward
		network.feedForward(input, output);

		std::vector<std::pair<double, int>> decision;
		for (unsigned i = 0; i < output.getRows(); i++){
			decision.push_back({ output(i, 0), i });
		}
		// Sort output
		std::sort(decision.begin(), decision.end());
		// Moves
		out << M[decision[3].second].c_str() << " with confidence: " << decision[3].first << std::endl;
		out << M[decision[2].second].c_str() << " " << decision[2].first << " " 
			<< M[decision[1].second].c_str() << " " << decision[1].first << " " 
			<< M[decision[0].second].c_str() << " " << decision[0].first << std::endl;
		if (gInst.Move(decision[3].second)){
			gInst.random();
		}
		else{
			break;
		}		
		
		out << std::endl;
		out << "Score: " << gInst.getScore() << std::endl;
		gInst.printGrid(out);		
		//Sleep(2000);		
	}
	out << "Score: " << gInst.getScore() << std::endl;
	out << std::endl;
	return gInst.getScore();
}