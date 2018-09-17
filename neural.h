#ifndef NET_H
#define NET_H
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <random>
#include "matrix.h"
#include <Windows.h>
#include <algorithm>
#include <fstream>
#include <string>

using namespace std;
typedef matrix<double> Matrix;

class Net
{
public:
	Net(const vector<unsigned> &topology);
	Net(const vector<unsigned> &topology, double range);
	Net(const vector<unsigned> &topology, char * fileName);

	void feedForward(const Matrix &inputVals, Matrix &resultVals);
	void feedForward(const Matrix &inputVals, Matrix &resultVals, vector<Matrix> &A);
	void backProp(const Matrix &inputVals, const Matrix &targetVals, vector<Matrix> &alpha, vector<Matrix> &delta);

	double computeCost(const Matrix &inputVals, const Matrix &targetVals);
	double computeClassifCost(const Matrix &inputVals, const Matrix &targetVals);
	void computeGradients(const Matrix &inputVals, const Matrix &targetVals, vector<Matrix> &grads);
	void checkGradients(const Matrix &inputVals, const Matrix &targetVals, vector<Matrix> &grads);
	void gradientDescent(const Matrix &inputVals, const Matrix &targetVals);
	static double sigmoid(double x);

	void save(char * fileName);
	void load(char * fileName);

	void setWeights(Net& firstParent, Net& secondParent);
	void setWeightsMean(Net& firstParent, Net& secondParent);
	void mutate();
	void mutate(int prop, double sigma);
private:
	void randWeights(Matrix &weightVals);
	void randWeights(Matrix &weightVals, double epsilon);
	void sigmoid(Matrix &mat);

	vector<unsigned> m_layers; 
	vector<Matrix> m_weights;
	vector<Matrix> m_sigma;
};

struct creature
{
	double fitness;
	double avgScore;
	Net network;
};

struct idCreature
{
	int id;
	int numChildren;
	double fitness;
	Net network;
};
#endif





