#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <vector>
#include <iostream>

using namespace std;
using Eigen::MatrixXf;

MatrixXf    ConvTrans_Vector2DEigen(const vector<vector<float>> &vector2D);
MatrixXf    Conv_Vector2DEigen(const vector<vector<float>> &vector2D);
MatrixXf    Conv_Vector9fEigen3f(const vector<float> &vector9f);
void        printEigen2D(const Eigen::MatrixXf &eigenMat);

int main () {
        std::vector<std::vector<float>> vector2D{
                                                {1,2,3,4,2,3,4},
                                                {3,4,5,6,4,5,6},
                                                {9,10,11,12,10,11,12}
                                                };
        std::vector<float> vector9f = {1,2,3,4,5,6,7,8,9};

    // test ConvTrans_Vector2DEigen
        MatrixXf eigen2D;
        eigen2D = ConvTrans_Vector2DEigen(vector2D);
                printEigen2D(eigen2D);

    // test Conv_Vector2DEigen
        eigen2D = Conv_Vector2DEigen(vector2D);
                printEigen2D(eigen2D);

    // test Conv_Vector9fEigen3f
        MatrixXf eigen3f;
        eigen3f = Conv_Vector9fEigen3f(vector9f);
                printEigen2D(eigen3f);

    return 0;
};

void printEigen2D(const Eigen::MatrixXf &eigen2D) {
    // print out the result
    cout << "---------------------------" << endl;
    for (int i{0}; i < eigen2D.rows(); ++i) {
        for (int j{0}; j < eigen2D.cols(); ++j) {
            cout << eigen2D(i, j) << " ";
        }
        cout << endl;
    }
}

MatrixXf ConvTrans_Vector2DEigen(const vector<vector<float>> &vector2D) {
    // vector2D[i][j] >> eigen2D(j, i)
        // convert and transpose in the same time

    MatrixXf eigen2D = Eigen::MatrixXf::Constant(vector2D[0].size(), vector2D.size(), 0.0f);

    for (int i{0}; i < vector2D.size(); ++i) {
        eigen2D.col(i) = Eigen::VectorXf::Map(&vector2D[i][0], vector2D[i].size());
    }
    
    return eigen2D;
}

MatrixXf Conv_Vector2DEigen(const vector<vector<float>> &vector2D) {
    // vector2D[i][j] >> eigen2D(i, j)
        // convert

    MatrixXf eigen2D = Eigen::MatrixXf::Constant(vector2D.size(), vector2D[0].size(), 0.0f);

    for (int i{0}; i < vector2D.size(); ++i) {
        eigen2D.row(i) = Eigen::VectorXf::Map(&vector2D[i][0], vector2D[i].size());
    }
    
    return eigen2D;
}

MatrixXf Conv_Vector9fEigen3f(const vector<float> &vector9f) {
    // vector2D[i] >> eigen2D(i, j)
        // convert 9 element vector 1D into 3x3 eigen, row-major
            // eigen3f  = Eigen::Map<Eigen::Matrix<float, 3, 3>> (&vector9f);
            // eigen3f  = Eigen::MatrixXf::Map(&vector9f, 3, 3);

    MatrixXf eigen3f    = Eigen::MatrixXf::Constant(3, 3, 0.0f);
    eigen3f.row(0)      = Eigen::VectorXf::Map(&vector9f[0], 3);
    eigen3f.row(1)      = Eigen::VectorXf::Map(&vector9f[3], 3);
    eigen3f.row(2)      = Eigen::VectorXf::Map(&vector9f[6], 3);

    return eigen3f;
}