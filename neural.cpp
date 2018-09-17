#include "neural.h"

void Net::save(char * fileName)
{
	std::ofstream out(fileName);

	for (unsigned l = 0; l < m_layers.size() - 1; ++l) {
		for (unsigned i = 0; i < m_weights[l].getRows(); ++i) {
			for (unsigned j = 0; j < m_weights[l].getColumns(); ++j) {
				out << m_weights[l](i, j) << ' ';
			}
			out << std::endl;
		}
		out << endl;
	}

	out.close();
}

void Net::load(char * fileName)
{
	std::ifstream in(fileName);
	std::streambuf *cinbuf = std::cin.rdbuf(); //save old buf
	std::cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!

	for (unsigned l = 0; l < m_layers.size() - 1; ++l) {
		for (unsigned i = 0; i < m_weights[l].getRows(); ++i) {
			for (unsigned j = 0; j < m_weights[l].getColumns(); ++j) {
				cin >> m_weights[l](i, j);
			}
		}
	}

	std::cin.rdbuf(cinbuf);   //reset to standard input again
}

double Net::sigmoid(double x)
{
	return 1.0 / (1.0 + exp(-x));
}

void Net::sigmoid(Matrix &mat)
{
	for (auto it = mat.begin(); it != mat.end(); ++it) {
		(*it) = Net::sigmoid(*it);
	}
}

void Net::feedForward(const Matrix &inputVals, Matrix &resultVals)
{
	Matrix alpha, zeta;
	// input vals to a_1 
	alpha = inputVals;

	// for all hidden layers
	for (unsigned i = 1; i < m_layers.size(); ++i) {
		//	add a_i-1_0
		zeta.resizeNoCopy(alpha.getRows() + 1, 1);
		zeta(0, 0) = 1.0;
		for (unsigned j = 1; j < zeta.getRows(); ++j) {
			zeta(j, 0) = alpha(j - 1, 0);
		}
		alpha = m_weights[i - 1] * zeta;
		sigmoid(alpha);
	}
	// return a_i
	resultVals = alpha;
}

void Net::feedForward(const Matrix &inputVals, Matrix &resultVals, vector<Matrix> &A)
{
	Matrix alpha, zeta;
	// input vals to a_1 
	alpha = inputVals;

	// for all hidden layers
	for (unsigned i = 1; i < m_layers.size(); ++i) {
		//	add a_i-1_0
		Matrix temp(alpha.getRows() + 1, 1);
		temp(0, 0) = 1.0;
		auto temp_it = ++(temp.begin());
		for (auto alpha_it = alpha.begin(); alpha_it != alpha.end(); ++alpha_it, ++temp_it) {
			(*temp_it) = (*alpha_it);
		}
		alpha = temp;
		A.push_back(alpha);
		//	z_i = theta_i-1 * a_i-1
		zeta = m_weights[i - 1] * alpha;
		//	a_i = sigmoid(z_i);
		alpha = zeta;
		sigmoid(alpha);
	}
	// return a_i
	A.push_back(alpha);
	resultVals = alpha;
}

double Net::computeCost(const Matrix &inputVals, const Matrix &targetVals)
{
	Matrix resultVals;
	feedForward(inputVals, resultVals);

	unsigned numOutputs = resultVals.getRows();
	double sum = 0.0;
	for (unsigned n = 0; n < numOutputs; ++n) {
		sum += 0.5 * (targetVals(n, 0) - resultVals(n, 0)) * (targetVals(n, 0) - resultVals(n, 0));
	}
	return sum;
}

double Net::computeClassifCost(const Matrix &inputVals, const Matrix &targetVals)
{
	Matrix resultVals;
	feedForward(inputVals, resultVals);

	unsigned numOutputs = resultVals.getRows();
	double sum = 0.0;
	for (unsigned n = 0; n < numOutputs; ++n) {
		sum += targetVals(n, 0) * log(resultVals(n, 0)) + (1 - targetVals(n, 0)) * log(1 - resultVals(n, 0));
	}
	return sum;
}

void Net::backProp(const Matrix &inputVals, const Matrix &targetVals, vector<Matrix> &alpha, vector<Matrix> &delta)
{
	Matrix resultVals;
	feedForward(inputVals, resultVals, alpha);

	unsigned numLayers = m_layers.size();
	// Last Layer
	// d_j = A[L]_j - targetVals_j
	delta[numLayers - 1] = (resultVals - targetVals);

	// Rest Layers
	for (unsigned l = numLayers - 2; l > 0; --l) {
		// d[l] = m_weights[l]^T*d[l+1]
		delta[l] = m_weights[l].trans() *delta[l + 1];
		// d[l] .*= alpha[l - 1].*(1-alpha[l - 1])
		for (unsigned i = 0; i < delta[l].getRows(); ++i) {
			delta[l](i, 0) = alpha[l - 1](i, 0) * (1 - alpha[l - 1](i, 0));
		}
	}
}

void Net::computeGradients(const Matrix &inputVals, const Matrix &targetVals, vector<Matrix> &grads)
{
	vector<Matrix> alpha, delta;
	for (unsigned i = 0; i < m_layers.size(); ++i) {
		delta.push_back({ { 0 } });
	}
	// compute alphas deltas (back propagation)
	backProp(inputVals, targetVals, alpha, delta);

	// grards[l].(i, j) = alpha[l]_j*delta[l+1]_i
	for (unsigned l = 0; l < m_layers.size() - 1; ++l) {
		grads.push_back(Matrix(m_layers[l + 1], m_layers[l]));
		for (unsigned i = 0; i < m_layers[l + 1]; ++i) {
			for (unsigned j = 0; j < m_layers[l]; ++j) {
				grads[l](i, j) = alpha[l](j, 0) * delta[l + 1](i, 0);
			}
		}
	}
}

void Net::checkGradients(const Matrix &inputVals, const Matrix &targetVals, vector<Matrix> &grads)
{
	double epsilon = 10e-4;

	double costPlusEpsilon, costMinusEpsilon;
	double originalWieght;

	// for all grads
	for (unsigned l = 0; l < m_layers.size() - 1; ++l) {
		grads.push_back(Matrix(m_layers[l + 1], m_layers[l]));
		for (unsigned i = 0; i < m_layers[l + 1]; ++i) {
			for (unsigned j = 0; j < m_layers[l]; ++j) {
				originalWieght = m_weights[l](i, j);
				m_weights[l](i, j) = originalWieght + epsilon;
				costPlusEpsilon = computeCost(inputVals, targetVals);

				m_weights[l](i, j) = originalWieght - epsilon;
				costMinusEpsilon = computeCost(inputVals, targetVals);
				m_weights[l](i, j) = originalWieght;

				grads[l](i, j) = (costPlusEpsilon - costMinusEpsilon) / (2.0 * epsilon);
			}
		}
	}
}

void Net::gradientDescent(const Matrix &inputVals, const Matrix &targetVals)
{
	vector<Matrix> grads, gradsCheck;
	//computeGradients(inputVals, targetVals, grads);
	checkGradients(inputVals, targetVals, gradsCheck);

	/*for (unsigned l = 0; l < m_layers.size() - 1; ++l) {
	for (unsigned i = 0; i < m_layers[l + 1].size(); ++i) {
	for (unsigned j = 0; j < m_layers[l].size(); ++j) {
	cout << "Grad: " << grads[l](i, j) << " Check: " << gradsCheck[l](i, j) << endl;
	cout << "Error: " << grads[l](i, j) - gradsCheck[l](i, j) << endl;	}	}	}*/

	double alpha = 2.0;
	for (unsigned l = 0; l < m_layers.size() - 1; ++l) {
		for (unsigned i = 0; i < m_layers[l + 1]; ++i) {
			for (unsigned j = 0; j < m_layers[l]; ++j) {
				m_weights[l](i, j) -= alpha * gradsCheck[l](i, j);
			}
		}
	}
}

void Net::randWeights(Matrix &weightVals)
{
	double epsilon = 7.12;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(-epsilon, epsilon);

	for (auto it = weightVals.begin(); it != weightVals.end(); ++it) {
		(*it) = dist(mt);
	}
}

void Net::randWeights(Matrix &weightVals, double epsilon)
{
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(-epsilon, epsilon);

	for (auto it = weightVals.begin(); it != weightVals.end(); ++it) {
		(*it) = dist(mt);
	}
}

Net::Net(const vector<unsigned> &topology)
{
	unsigned numLayers = topology.size();
	m_layers = topology;
	for (unsigned layerNum = 0; layerNum < numLayers; ++layerNum) {
		unsigned numOutputs = layerNum == topology.size() - 1 ? 0 : topology[layerNum + 1];

		if (numOutputs > 0)	{
			// W matrix (Layer i+1 size) x (Layer i size + 1)		
			m_weights.push_back(Matrix(numOutputs, topology[layerNum] + 1));
			// Randomize weights
			randWeights(m_weights.back());
		}
	}
}

Net::Net(const vector<unsigned> &topology, double range)
{
	unsigned numLayers = topology.size();
	m_layers = topology;
	for (unsigned layerNum = 0; layerNum < numLayers; ++layerNum) {
		unsigned numOutputs = layerNum == topology.size() - 1 ? 0 : topology[layerNum + 1];

		if (numOutputs > 0)	{
			// W matrix (Layer i+1 size) x (Layer i size + 1)		
			m_weights.push_back(Matrix(numOutputs, topology[layerNum] + 1));
			// Randomize weights
			randWeights(m_weights.back(), range);
		}
	}
}

Net::Net(const vector<unsigned> &topology, char * fileName)
{
	unsigned numLayers = topology.size();
	m_layers = topology;
	for (unsigned layerNum = 0; layerNum < numLayers; ++layerNum) {
		unsigned numOutputs = layerNum == topology.size() - 1 ? 0 : topology[layerNum + 1];

		if (numOutputs > 0)	{
			// W matrix (Layer i+1 size) x (Layer i size + 1)		
			m_weights.push_back(Matrix(numOutputs, topology[layerNum] + 1));
		}
	}

	load(fileName);
}

void Net::setWeightsMean(Net& firstParent, Net& secondParent)
{
	for (unsigned l = 0; l < m_layers.size() - 1; ++l) {
		for (unsigned i = 0; i < m_layers[l + 1] - 1; ++i) {
			for (unsigned j = 0; j < m_layers[l]; ++j) {
				m_weights[l](i, j) = (firstParent.m_weights[l](i, j) + secondParent.m_weights[l](i, j)) * 0.5;
			}
		}
	}
}

void Net::setWeights(Net& firstParent, Net& secondParent)
{
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> dist(0, 1);
	for (unsigned l = 0; l < m_layers.size() - 1; ++l) {
		for (unsigned i = 0; i < m_layers[l + 1] - 1; ++i) {
			for (unsigned j = 0; j < m_layers[l]; ++j) {
				m_weights[l](i, j) = dist(mt) ? firstParent.m_weights[l](i, j) : secondParent.m_weights[l](i, j);
			}
		}
	}
}

void Net::mutate()
{
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> dist_int(0, 100);
	std::uniform_real_distribution<double> dist_real(-2.0, 2.0);
	for (unsigned l = 0; l < m_layers.size() - 1; ++l) {
		for (unsigned i = 0; i < m_layers[l + 1] - 1; ++i) {
			for (unsigned j = 0; j < m_layers[l]; ++j) {
				dist_int(mt) < 2 ? m_weights[l](i, j) += dist_real(mt) : 0;
			}
		}
	}
}

void Net::mutate(int prop, double sigma)
{
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> dist_int(0, prop);
	std::uniform_real_distribution<double> dist_real(-sigma, sigma);
	for (unsigned l = 0; l < m_layers.size() - 1; ++l) {
		for (unsigned i = 0; i < m_layers[l + 1] - 1; ++i) {
			for (unsigned j = 0; j < m_layers[l]; ++j) {
				dist_int(mt) < 2 ? m_weights[l](i, j) += dist_real(mt) : 0;
			}
		}
	}
}


