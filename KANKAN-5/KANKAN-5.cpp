//Concept: Andrew Polar and Mike Poluektov
//Developer Andrew Polar

// License
// If the end user somehow manages to make billions of US dollars using this code,
// and happens to meet the developer begging for change outside a McDonald's,
// he or she is under no obligation to buy the developer a sandwich.

// Symmetry Clause
// Likewise, if the developer becomes rich and famous by publishing this code,
// and meets an unfortunate end user who went bankrupt using it,
// the developer is also under no obligation to buy the end user a sandwich.

//Publications:
//https://www.sciencedirect.com/science/article/abs/pii/S0016003220301149
//https://www.sciencedirect.com/science/article/abs/pii/S0952197620303742
//https://arxiv.org/abs/2305.08194

#include <iostream>
#include "Helper.h"
#include "Urysohn.h"
#include "Layer.h"

///////////// Determinat dataset
std::vector<std::vector<double>> GenerateInput(int nRecords, int nFeatures, double min, double max) {
	std::vector<std::vector<double>> x(nRecords);
	for (int i = 0; i < nRecords; ++i) {
		x[i] = std::vector<double>(nFeatures);
		for (int j = 0; j < nFeatures; ++j) {
			x[i][j] = static_cast<double>((rand() % 10000) / 10000.0);
			x[i][j] *= (max - min);
			x[i][j] += min;
		}
	}
	return x;
}

double determinant(const std::vector<std::vector<double>>& matrix) {
	int n = (int)matrix.size();
	if (n == 1) {
		return matrix[0][0];
	}
	if (n == 2) {
		return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
	}
	double det = 0.0;
	for (int col = 0; col < n; ++col) {
		std::vector<std::vector<double>> subMatrix(n - 1, std::vector<double>(n - 1));
		for (int i = 1; i < n; ++i) {
			int subCol = 0;
			for (int j = 0; j < n; ++j) {
				if (j == col) continue;
				subMatrix[i - 1][subCol++] = matrix[i][j];
			}
		}
		det += (col % 2 == 0 ? 1 : -1) * matrix[0][col] * determinant(subMatrix);
	}
	return det;
}

double ComputeDeterminant(const std::vector<double>& input, int N) {
	std::vector<std::vector<double>> matrix(N, std::vector<double>(N, 0.0));
	int cnt = 0;
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			matrix[i][j] = input[cnt++];
		}
	}
	return determinant(matrix);
}

std::vector<double> ComputeDeterminantTarget(const std::vector<std::vector<double>>& x, int nMatrixSize) {
	int nRecords = (int)x.size();
	std::vector<double> target(nRecords);
	int counter = 0;
	while (true) {
		target[counter] = ComputeDeterminant(x[counter], nMatrixSize);
		if (++counter >= nRecords) break;
	}
	return target;
}
///////// End determinant data

///////// Areas of faces of tetrahedron
double Area(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) {
	double a1 = (y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1);
	double a2 = (x2 - x1) * (z3 - z1) - (z2 - z1) * (x3 - x1);
	double a3 = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
	double A = 0.5 * sqrt(a1 * a1 + a2 * a2 + a3 * a3);
	return A;
}

std::vector<std::vector<double>> MakeRandomMatrix(int rows, int cols, double min, double max) {
	std::vector<std::vector<double>> matrix(rows);
	for (int i = 0; i < rows; ++i) {
		matrix[i] = std::vector<double>(cols);
		for (int j = 0; j < cols; ++j) {
			matrix[i][j] = static_cast<double>((rand() % 1000) / 1000.0) * (max - min) + min;
		}
	}
	return matrix;
}

std::vector<std::vector<double>> ComputeTargetMatrix(const std::vector<std::vector<double>>& X) {
	int rows = (int)X.size();
	std::vector<std::vector<double>> matrix(rows);
	for (int i = 0; i < rows; ++i) {
		matrix[i] = std::vector<double>(4);
		matrix[i][0] = Area(X[i][0], X[i][1], X[i][2], X[i][3], X[i][4], X[i][5], X[i][6], X[i][7], X[i][8]);
		matrix[i][1] = Area(X[i][0], X[i][1], X[i][2], X[i][3], X[i][4], X[i][5], X[i][9], X[i][10], X[i][11]);
		matrix[i][2] = Area(X[i][0], X[i][1], X[i][2], X[i][6], X[i][7], X[i][8], X[i][9], X[i][10], X[i][11]);
		matrix[i][3] = Area(X[i][3], X[i][4], X[i][5], X[i][6], X[i][7], X[i][8], X[i][9], X[i][10], X[i][11]);
	}
	return matrix;
}
//////////// End tetrahedron

///////// Medians
double Median1(double x1, double y1, double x2, double y2, double x3, double y3) {
	double t1 = x1 - (x2 + x3) / 2.0;
	double t2 = y1 - (y2 + y3) / 2.0;
	t1 *= t1;
	t2 *= t2;
	return sqrt(t1 + t2);
}

double Median2(double x1, double y1, double x2, double y2, double x3, double y3)
{
	double t1 = x2 - (x1 + x3) / 2.0;
	double t2 = y2 - (y1 + y3) / 2.0;
	t1 *= t1;
	t2 *= t2;
	return sqrt(t1 + t2);
}

double Median3(double x1, double y1, double x2, double y2, double x3, double y3)
{
	double t1 = x3 - (x2 + x1) / 2.0;
	double t2 = y3 - (y2 + y1) / 2.0;
	t1 *= t1;
	t2 *= t2;
	return sqrt(t1 + t2);
}

std::vector<std::vector<double>> GenerateInputsMedians(int nRecords, int nFeatures, double min, double max) {
	std::vector<std::vector<double>> x(nRecords);
	for (int i = 0; i < nRecords; ++i) {
		x[i] = std::vector<double>(nFeatures);
		for (int j = 0; j < nFeatures; ++j) {
			x[i][j] = static_cast<double>((rand() % 10000) / 10000.0);
			x[i][j] *= (max - min);
			x[i][j] += min;
		}
	}
	return x;
}

std::vector<std::vector<double>> ComputeTargetsMedians(const std::vector<std::vector<double>>& x) {
	int nRecords = (int)x.size();
	std::vector<std::vector<double>> y(nRecords);
	for (int i = 0; i < nRecords; ++i) {
		y[i] = std::vector<double>(3);
		for (int j = 0; j < 3; ++j) {
			y[i][0] = Median1(x[i][0], x[i][1], x[i][2], x[i][3], x[i][4], x[i][5]);
			y[i][1] = Median2(x[i][0], x[i][1], x[i][2], x[i][3], x[i][4], x[i][5]);
			y[i][2] = Median3(x[i][0], x[i][1], x[i][2], x[i][3], x[i][4], x[i][5]);
		}
	}
	return y;
}
///////// End medians

///////// Random triangles
std::vector<std::vector<double>> MakeRandomMatrixForTriangles(int rows, int cols, double min, double max) {
	std::vector<std::vector<double>> matrix(rows);
	for (int i = 0; i < rows; ++i) {
		matrix[i] = std::vector<double>(cols);
		for (int j = 0; j < cols; ++j) {
			matrix[i][j] = static_cast<double>((rand() % 1000) / 1000.0) * (max - min) + min;
		}
	}
	return matrix;
}
double AreaOfTriangle(double x1, double y1, double x2, double y2, double x3, double y3) {
	double A = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
	return A;
}
std::vector<std::vector<double>> ComputeAreasOfTriangles(std::vector<std::vector<double>>& matrix) {
	int N = (int)matrix.size();
	std::vector<std::vector<double>> u(N);
	for (int i = 0; i < N; ++i) {
		u[i] = std::vector<double>(1);
		u[i][0] = AreaOfTriangle(matrix[i][0], matrix[i][1], matrix[i][2], matrix[i][3], matrix[i][4], matrix[i][5]);
	}
	return u;
}
///////// End of random triangles

//We don't need 4 layers here, it is only a demo how to make 4 layers
void AreasOfTriangles() {
	int nFeatures = 6;
	int nTargets = 1;
	int nTrainingRecords = 10000;
	int nValidationRecords = 2000;
	auto features_training = MakeRandomMatrixForTriangles(nTrainingRecords, nFeatures, 0.0, 1.0);
	auto features_validation = MakeRandomMatrixForTriangles(nValidationRecords, nFeatures, 0.0, 1.0);
	auto targets_training = ComputeAreasOfTriangles(features_training);
	auto targets_validation = ComputeAreasOfTriangles(features_validation);

	//data is ready, we start training
	clock_t start_application = clock();
	clock_t current_time = clock();

	std::vector<double> argmin;
	std::vector<double> argmax;
	Helper::FindMinMaxMatrix2(argmin, argmax, features_training);

	double targetMin = DBL_MAX;
	double targetMax = -DBL_MAX;
	for (int i = 0; i < nTrainingRecords; ++i) {
		for (int j = 0; j < nTargets; ++j) {
			if (targets_training[i][j] < targetMin) targetMin = targets_training[i][j];
			if (targets_training[i][j] > targetMax) targetMax = targets_training[i][j];
		}
	}

	int nU0 = 50;
	int nU1 = 8;
	int nU2 = 4;
	int nU3 = nTargets;
	double alpha = 0.05;

	auto layer0 = std::make_unique<Layer>(nU0, argmin, argmax, targetMin, targetMax, 2);
	auto layer1 = std::make_unique<Layer>(nU1, nU0, targetMin, targetMax, 12);
	auto layer2 = std::make_unique<Layer>(nU2, nU1, targetMin, targetMax, 12);
	auto layer3 = std::make_unique<Layer>(nU3, nU2, targetMin, targetMax, 22);

	std::vector<double> models0(nU0);
	std::vector<double> models1(nU1);
	std::vector<double> models2(nU2);
	std::vector<double> models3(nU3);

	std::vector<double> deltas3(nU3);
	std::vector<double> deltas2(nU2);
	std::vector<double> deltas1(nU1);
	std::vector<double> deltas0(nU0);

	std::vector<std::vector<double>> derivatives1(nU1, std::vector<double>(nU0, 0.0));
	std::vector<std::vector<double>> derivatives2(nU2, std::vector<double>(nU1, 0.0));
	std::vector<std::vector<double>> derivatives3(nU3, std::vector<double>(nU2, 0.0));

	auto actual0 = std::make_unique<double[]>(nValidationRecords);
	auto computed0 = std::make_unique<double[]>(nValidationRecords);

	printf("Targets are areas of random triangles, %d training records\n", nTrainingRecords);
	for (int epoch = 0; epoch < 128; ++epoch) {
		for (int i = 0; i < nTrainingRecords; ++i) {
			layer0->Input2Output(features_training[i], models0);
			layer1->Input2Output(models0, models1, derivatives1);
			layer2->Input2Output(models1, models2, derivatives2);
			layer3->Input2Output(models2, models3, derivatives3);

			for (int j = 0; j < nTargets; ++j) {
				deltas3[j] = (targets_training[i][j] - models3[j]) * alpha;
			}

			layer3->ComputeDeltas(derivatives3, deltas3, deltas2);
			layer2->ComputeDeltas(derivatives2, deltas2, deltas1);
			layer1->ComputeDeltas(derivatives1, deltas1, deltas0);

			layer3->Update(models2, deltas3);
			layer2->Update(models1, deltas2);
			layer1->Update(models0, deltas1);
			layer0->Update(features_training[i], deltas0);
		}

		double error = 0.0;
		for (int i = 0; i < nValidationRecords; ++i) {
			layer0->Input2Output(features_validation[i], models0);
			layer1->Input2Output(models0, models1);
			layer2->Input2Output(models1, models2);
			layer3->Input2Output(models2, models3);

			for (int j = 0; j < nTargets; ++j) {
				double err = targets_validation[i][j] - models3[j];
				error += err * err;
			}

			actual0[i] = targets_validation[i][0];
			computed0[i] = models3[0];
		}

		//pearsons for correlated targets
		double p1 = Helper::Pearson(computed0, actual0, nValidationRecords);

		//mean error
		error /= nTargets;
		error /= nValidationRecords;
		error = sqrt(error);
		current_time = clock();
		printf("Epoch %d, RMSE %f, Pearson: %f, time %2.3f\n", epoch, error, p1,
			(double)(current_time - start_application) / CLOCKS_PER_SEC);

		if (p1 > 0.985) break;
	}
	printf("\n");
}

void Medians() {
	int nTrainingRecords = 10000;
	int nValidationRecords = 2000;
	int nFeatures = 6;
	int nTargets = 3;
	double min = 0.0;
	double max = 1.0;
	auto features_training = GenerateInputsMedians(nTrainingRecords, nFeatures, min, max);
	auto features_validation = GenerateInputsMedians(nValidationRecords, nFeatures, min, max);
	auto targets_training = ComputeTargetsMedians(features_training);
	auto targets_validation = ComputeTargetsMedians(features_validation);

	//data is ready, we start training
	clock_t start_application = clock();
	clock_t current_time = clock();

	std::vector<double> argmin;
	std::vector<double> argmax;
	Helper::FindMinMaxMatrix2(argmin, argmax, features_training);

	double targetMin = DBL_MAX;
	double targetMax = -DBL_MAX;
	for (int i = 0; i < nTrainingRecords; ++i) {
		for (int j = 0; j < nTargets; ++j) {
			if (targets_training[i][j] < targetMin) targetMin = targets_training[i][j];
			if (targets_training[i][j] > targetMax) targetMax = targets_training[i][j];
		}
	}

	int nU0 = 20;
	int nU1 = 10;
	int nU2 = 4;
	int nU3 = nTargets;
	double alpha = 0.05;

	auto layer0 = std::make_unique<Layer>(nU0, argmin, argmax, targetMin, targetMax, 2);
	auto layer1 = std::make_unique<Layer>(nU1, nU0, targetMin, targetMax, 12);
	auto layer2 = std::make_unique<Layer>(nU2, nU1, targetMin, targetMax, 12);
	auto layer3 = std::make_unique<Layer>(nU3, nU2, targetMin, targetMax, 22);

	std::vector<double> models0(nU0);
	std::vector<double> models1(nU1);
	std::vector<double> models2(nU2);
	std::vector<double> models3(nU3);

	std::vector<double> deltas3(nU3);
	std::vector<double> deltas2(nU2);
	std::vector<double> deltas1(nU1);
	std::vector<double> deltas0(nU0);

	std::vector<std::vector<double>> derivatives1(nU1, std::vector<double>(nU0, 0.0));
	std::vector<std::vector<double>> derivatives2(nU2, std::vector<double>(nU1, 0.0));
	std::vector<std::vector<double>> derivatives3(nU3, std::vector<double>(nU2, 0.0));

	auto actual0 = std::make_unique<double[]>(nValidationRecords);
	auto actual1 = std::make_unique<double[]>(nValidationRecords);
	auto actual2 = std::make_unique<double[]>(nValidationRecords);

	auto computed0 = std::make_unique<double[]>(nValidationRecords);
	auto computed1 = std::make_unique<double[]>(nValidationRecords);
	auto computed2 = std::make_unique<double[]>(nValidationRecords);

	printf("Targets are medians of random triangles, %d training records\n", nTrainingRecords);
	for (int epoch = 0; epoch < 128; ++epoch) {
		for (int i = 0; i < nTrainingRecords; ++i) {
			layer0->Input2Output(features_training[i], models0);
			layer1->Input2Output(models0, models1, derivatives1);
			layer2->Input2Output(models1, models2, derivatives2);
			layer3->Input2Output(models2, models3, derivatives3);

			for (int j = 0; j < nTargets; ++j) {
				deltas3[j] = (targets_training[i][j] - models3[j]) * alpha;
			}

			layer3->ComputeDeltas(derivatives3, deltas3, deltas2);
			layer2->ComputeDeltas(derivatives2, deltas2, deltas1);
			layer1->ComputeDeltas(derivatives1, deltas1, deltas0);

			layer3->Update(models2, deltas3);
			layer2->Update(models1, deltas2);
			layer1->Update(models0, deltas1);
			layer0->Update(features_training[i], deltas0);
		}

		double error = 0.0;
		for (int i = 0; i < nValidationRecords; ++i) {
			layer0->Input2Output(features_validation[i], models0);
			layer1->Input2Output(models0, models1);
			layer2->Input2Output(models1, models2);
			layer3->Input2Output(models2, models3);

			for (int j = 0; j < nTargets; ++j) {
				double err = targets_validation[i][j] - models3[j];
				error += err * err;
			}

			actual0[i] = targets_validation[i][0];
			actual1[i] = targets_validation[i][1];
			actual2[i] = targets_validation[i][2];

			computed0[i] = models3[0];
			computed1[i] = models3[1];
			computed2[i] = models3[2];
		}

		//pearsons for correlated targets
		double p1 = Helper::Pearson(computed0, actual0, nValidationRecords);
		double p2 = Helper::Pearson(computed1, actual1, nValidationRecords);
		double p3 = Helper::Pearson(computed2, actual2, nValidationRecords);

		//mean error
		error /= nTargets;
		error /= nValidationRecords;
		error = sqrt(error);
		current_time = clock();
		printf("Epoch %d, RMSE %f, Pearsons: %f, %f, %f, time %2.3f\n", epoch, error, p1, p2, p3,
			(double)(current_time - start_application) / CLOCKS_PER_SEC);

		if (p1 > 0.985 && p2 > 0.985 && p3 > 0.985) break;
	}
	printf("\n");
}

//Demo how to use Layers without KANKAN wrapper
void Det_4_4() {
	int nTrainingRecords = 100000;
	int nValidationRecords = 20000;
	int nMatrixSize = 4;
	int nFeatures = nMatrixSize * nMatrixSize;
	int nTargets = 1;
	double min = 0.0;
	double max = 1.0;
	auto features_training = GenerateInput(nTrainingRecords, nFeatures, min, max);
	auto features_validation = GenerateInput(nValidationRecords, nFeatures, min, max);
	auto targets_training = ComputeDeterminantTarget(features_training, nMatrixSize);
	auto targets_validation = ComputeDeterminantTarget(features_validation, nMatrixSize);

	clock_t start_application = clock();
	clock_t current_time = clock();

	//find limits
	std::vector<double> argmin;
	std::vector<double> argmax;
	double targetMin;
	double targetMax;
	Helper::FindMinMax2(argmin, argmax, targetMin, targetMax, features_training, targets_training);

	//normalize targets, it is not necessary, but sometimes converges faster
	for (int i = 0; i < nTrainingRecords; ++i) {
		targets_training[i] = (targets_training[i] - targetMin) / (targetMax - targetMin);
	}
	for (int i = 0; i < nValidationRecords; ++i) {
		targets_validation[i] = (targets_validation[i] - targetMin) / (targetMax - targetMin);
	}

	//configuration
	int nU0 = 64;
	int nU1 = 1;
	double alpha = 0.4;

	//instantiation of layers
	auto layer0 = std::make_unique<Layer>(nU0, argmin, argmax, targetMin, targetMax, 3);
	auto layer1 = std::make_unique<Layer>(nU1, nU0, targetMin, targetMax, 30);

	//auxiliary data buffers for a quick moving data between methods
	std::vector<double> models0(nU0);
	std::vector<double> models1(nU1);

	std::vector<double> deltas1(nU1);
	std::vector<double> deltas0(nU0);

	std::vector<std::vector<double>> derivatives1(nU1, std::vector<double>(nU0, 0.0));
	auto actual_validation = std::make_unique<double[]>(nValidationRecords);

	//training
	printf("Targets are determinants of random 4 * 4 matrices, %d training records\n", nTrainingRecords);
	for (int epoch = 0; epoch < 128; ++epoch) {
		for (int i = 0; i < nTrainingRecords; ++i) {
			//forward feeding by two layers
			layer0->Input2Output(features_training[i], models0);
			layer1->Input2Output(models0, models1, derivatives1);

			//computing residual error
			for (int j = 0; j < nTargets; ++j) {
				deltas1[j] = (targets_training[i] - models1[j]) * alpha;
			}

			//back propagation
			layer1->ComputeDeltas(derivatives1, deltas1, deltas0);

			//updating of two layers
			layer1->Update(models0, deltas1);
			layer0->Update(features_training[i], deltas0);
		}

		//validation at the end of each epoch
		double error = 0.0;
		for (int i = 0; i < nValidationRecords; ++i) {
			layer0->Input2Output(features_validation[i], models0);
			layer1->Input2Output(models0, models1);
			actual_validation[i] = models1[0];
			error += (targets_validation[i] - models1[0]) * (targets_validation[i] - models1[0]);
		}
		double pearson = Helper::Pearson(targets_validation, actual_validation, nValidationRecords);
		error /= nValidationRecords;
		error = sqrt(error);
		current_time = clock();
		printf("Epoch %d, current relative error %f, pearson %f, time %2.3f\n", epoch, error, pearson, (double)(current_time - start_application) / CLOCKS_PER_SEC);

		if (pearson > 0.97) break;
	}
	printf("\n");
}

//Here I show how to use Layers directly without KANKAN wrapper
void Tetrahedron() {
	const int nTrainingRecords = 500000;
	const int nValidationRecords = 50000;
	const int nFeatures = 12;
	const int nTargets = 4;
	const double min = 0.0;
	const double max = 1.0;
	auto features_training = MakeRandomMatrix(nTrainingRecords, nFeatures, min, max);
	auto features_validation = MakeRandomMatrix(nValidationRecords, nFeatures, min, max);
	auto targets_training = ComputeTargetMatrix(features_training);
	auto targets_validation = ComputeTargetMatrix(features_validation);

	//data is ready, we start training
	clock_t start_application = clock();
	clock_t current_time = clock();

	std::vector<double> argmin;
	std::vector<double> argmax;
	Helper::FindMinMaxMatrix2(argmin, argmax, features_training);

	double targetMin = DBL_MAX;
	double targetMax = -DBL_MAX;
	for (int i = 0; i < nTrainingRecords; ++i) {
		for (int j = 0; j < nTargets; ++j) {
			if (targets_training[i][j] < targetMin) targetMin = targets_training[i][j];
			if (targets_training[i][j] > targetMax) targetMax = targets_training[i][j];
		}
	}

	int nU0 = 50;
	int nU1 = 10;
	int nU2 = nTargets;
	double alpha = 0.1;

	auto layer0 = std::make_unique<Layer>(nU0, argmin, argmax, targetMin, targetMax, 2);
	auto layer1 = std::make_unique<Layer>(nU1, nU0, targetMin, targetMax, 12);
	auto layer2 = std::make_unique<Layer>(nU2, nU1, targetMin, targetMax, 22);

	std::vector<double> models0(nU0);
	std::vector<double> models1(nU1);
	std::vector<double> models2(nU2);

	std::vector<double> deltas2(nU2);
	std::vector<double> deltas1(nU1);
	std::vector<double> deltas0(nU0);

	std::vector<std::vector<double>> derivatives1(nU1, std::vector<double>(nU0, 0.0));
	std::vector<std::vector<double>> derivatives2(nU2, std::vector<double>(nU1, 0.0));

	auto actual0 = std::make_unique<double[]>(nValidationRecords);
	auto actual1 = std::make_unique<double[]>(nValidationRecords);
	auto actual2 = std::make_unique<double[]>(nValidationRecords);
	auto actual3 = std::make_unique<double[]>(nValidationRecords);

	auto computed0 = std::make_unique<double[]>(nValidationRecords);
	auto computed1 = std::make_unique<double[]>(nValidationRecords);
	auto computed2 = std::make_unique<double[]>(nValidationRecords);
	auto computed3 = std::make_unique<double[]>(nValidationRecords);

	printf("Targets are areas of faces of random tetrahedrons, %d\n", nTrainingRecords);
	for (int epoch = 0; epoch < 128; ++epoch) {
		for (int i = 0; i < nTrainingRecords; ++i) {
			layer0->Input2Output(features_training[i], models0);
			layer1->Input2Output(models0, models1, derivatives1);
			layer2->Input2Output(models1, models2, derivatives2);

			for (int j = 0; j < nTargets; ++j) {
				deltas2[j] = (targets_training[i][j] - models2[j]) * alpha;
			}

			layer2->ComputeDeltas(derivatives2, deltas2, deltas1);
			layer1->ComputeDeltas(derivatives1, deltas1, deltas0);

			layer2->Update(models1, deltas2);
			layer1->Update(models0, deltas1);
			layer0->Update(features_training[i], deltas0);
		}

		double error = 0.0;
		for (int i = 0; i < nValidationRecords; ++i) {
			layer0->Input2Output(features_validation[i], models0);
			layer1->Input2Output(models0, models1);
			layer2->Input2Output(models1, models2);

			for (int j = 0; j < nTargets; ++j) {
				double err = targets_validation[i][j] - models2[j];
				error += err * err;
			}

			actual0[i] = targets_validation[i][0];
			actual1[i] = targets_validation[i][1];
			actual2[i] = targets_validation[i][2];
			actual3[i] = targets_validation[i][3];

			computed0[i] = models2[0];
			computed1[i] = models2[1];
			computed2[i] = models2[2];
			computed3[i] = models2[3];
		}
		double p1 = Helper::Pearson(computed0, actual0, nValidationRecords);
		double p2 = Helper::Pearson(computed1, actual1, nValidationRecords);
		double p3 = Helper::Pearson(computed2, actual2, nValidationRecords);
		double p4 = Helper::Pearson(computed3, actual3, nValidationRecords);

		error /= nTargets;
		error /= nValidationRecords;
		error = sqrt(error);
		current_time = clock();
		printf("Epoch %d, RMSE %f, Pearsons: %f %f %f %f, time %2.3f\n", epoch, error, p1, p2, p3, p4,
			(double)(current_time - start_application) / CLOCKS_PER_SEC);

		if (p1 > 0.975 && p2 > 0.975 && p3 > 0.975 && p4 > 0.975) break;
	}
	printf("\n");
}

int main() {
	srand((unsigned int)time(NULL));

	//Areas of random triangles.
	AreasOfTriangles();

	//Related targets, the medians of random triangles.
	Medians();

	//This simple unit test, features are random matrices of 4 by 4, targets are their determinants.
	Det_4_4();

	//Related targets, the areas of the faces of tetrahedron given by random vertices.
	Tetrahedron();
}

