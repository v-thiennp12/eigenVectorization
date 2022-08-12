#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <vector>
#include <iostream>

using namespace std;

int main () {

    //https://itecnote.com/tecnote/c-how-to-convert-an-stdvector-to-a-matrix-in-eigen/
        //version 1
        std::vector<float> test_vector1 = { 2,1,3 };
        float* test_array = test_vector1.data();
        Eigen::MatrixXf test1 = Eigen::Map<Eigen::Matrix<float, 3, 1>>(test_array);

        //version 2
        std::vector<float> test_vector2 = { 2,1,3 };
        Eigen::MatrixXf test2 = Eigen::Map<Eigen::Matrix<float, 3, 1>>(test_vector2.data());

        //version 3
        std::vector<float> test_vector3 = { 2,1,3 };
        Eigen::Map<Eigen::Matrix<float, 3, 1> > dangerousVec (test_vector3.data());

        // //version 4 not work
        // std::vector<float> test_vector4 = { 2,1,3 };
        // Eigen::MatrixXf test4 = Eigen::Map<Eigen::Matrix<float, 3, 1> >(test_vector4);


    //https://stackoverflow.com/questions/26094379/typecasting-eigenvectorxd-to-stdvector
    //init a first vector
    std::vector<float> v1;
    v1.push_back(0.5);
    v1.push_back(1.5);
    v1.push_back(2.5);
    v1.push_back(3.5);

        //from v1 to an eignen vector
        float* ptr_data = &v1[0];
        Eigen::VectorXf v2 = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(v1.data(), v1.size());

        Eigen::VectorXf v4 = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(ptr_data, v1.size());

        //from the eigen vector to the std vector
        std::vector<float> v3(&v2[0], v2.data()+v2.cols()*v2.rows());

        //to check
        for(int i = 0; i < v1.size() ; i++){
            std::cout << std::to_string(v1[i]) << " | " << std::to_string(v2[i]) << " | " << std::to_string(v3[i]) << " | " << std::to_string(v4[i]) << " | " << std::endl;
        }

    //2D vector to eigen
    std::vector<vector<float>>  vec2D{
                                        {1,2,3,4,2,3,4},
                                        {3,4,5,6,4,5,6},
                                        {9,10,11,12,10,11,12}
                                     };

    for (auto vec:vec2D) {
        for (auto val:vec) {
            cout << val << " ";
        }
        cout << endl;
    }

    // float                       *ptr_vec2D = &vec2D[0][0];
    // std::vector<vector<float>>     *ptr_vec2D = &vec2D;
    // float* ptr_vec2D                    = &vec2D.data();
    // Eigen::MatrixXf eig_vec2D;
    // // eig_vec2D = Eigen::Map<Eigen::Matrix<float, 3, 7, Eigen::RowMajor>>(ptr_vec2D);
    // eig_vec2D = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(vec2D.data());
    // // eig_vec2D = Eigen::Map<Eigen::Matrix<float, 3, 7, Eigen::RowMajor>>(vec2D);
    //     // does not work with 2D vector

    // cout << "---------------------------" << endl;
    // for (int i{0}; i < eig_vec2D.rows(); ++i) {
    //     for (int j{0}; j < eig_vec2D.cols(); ++j) {
    //         cout << eig_vec2D(i, j) << " ";
    //     }
    //     cout << endl;
    // }

    // ------------------------------------------------------
    Eigen::MatrixXf eig_vec2D_2(vec2D[0].size(), vec2D.size());
    for (int i{0}; i < vec2D.size(); ++i) {
        eig_vec2D_2.col(i) = Eigen::VectorXf::Map(&vec2D[i][0], vec2D[i].size());

        // float* ptr_data_2D = &vec2D[i][0];
        // eig_vec2D_2.row(i) = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned> (ptr_data_2D, vec2D[i].size());  
    }
 
    cout << "---------------------------" << endl;
    for (int i{0}; i < eig_vec2D_2.rows(); ++i) {
        for (int j{0}; j < eig_vec2D_2.cols(); ++j) {
            cout << eig_vec2D_2(i, j) << " ";
        }
        cout << endl;
    }

    return 0;
};