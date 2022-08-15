#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <vector>
#include <iostream>

using namespace std;
using Eigen::MatrixXf;
using Eigen::VectorXf;

MatrixXf    ConvTrans_Vector2DEigen(const vector<vector<float>> &vector2D);
MatrixXf    Conv_Vector2DEigen(const vector<vector<float>> &vector2D);
MatrixXf    Conv_Vector9fEigen3f(const vector<float> &vector9f);
MatrixXf    Conv_Vector1DEigenVec(const vector<float> &vector1D);
void        printEigen2D(const Eigen::MatrixXf &eigenMat);

struct camIFs_strct {
    vector<float> matrixD       = {1,2,3,4};
    vector<float> matrixK       = {1,2,3,4};
    vector<float> matrixR       = {1,2,3,4,5,6,7,8,9};
    vector<float> vectT         = {1,2,3};
};

int main () {
    //     std::vector<std::vector<float>> vector2D{
    //                                             {1,2,3,4,2,3,4},
    //                                             {3,4,5,6,4,5,6},
    //                                             {9,10,11,12,10,11,12}
    //                                             };
    //     std::vector<float> vector9f = {1,2,3,4,5,6,7,8,9};

    // // test ConvTrans_Vector2DEigen
    //     MatrixXf eigen2D;
    //     eigen2D = ConvTrans_Vector2DEigen(vector2D);
    //             printEigen2D(eigen2D);

    // // test Conv_Vector2DEigen
    //     eigen2D = Conv_Vector2DEigen(vector2D);
    //             printEigen2D(eigen2D);

    // // test Conv_Vector9fEigen3f
    //     MatrixXf eigen3f;
    //     eigen3f = Conv_Vector9fEigen3f(vector9f);
    //             // printEigen2D(eigen3f);

    //------------------------
    const int playGroundVerticesSize  = 5;
    vector<vector<float>>       playGroundVertices{
                                        {1,2,3},
                                        {3,4,5},
                                        {5,6,7},
                                        {7,8,9},
                                        {9,10,11}
                                     };
    int i=0;    

    vector<camIFs_strct> camIFs;
    camIFs.push_back(camIFs_strct());
    camIFs.push_back(camIFs_strct());
    vector<float> ImageOffset   = {0.1, 0.2};
    float PLAYGROUND_TEXTURE_HEIGHT = 5;
    float PLAYGROUND_TEXTURE_WIDTH  = 2.5;
    float RAW_IMG_HEIGHT            = 600;
    float RAW_IMG_WIDTH             = 800;
    
    vector<vector<float>> uvaOuts;
    // Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>

    Eigen::MatrixXf 				eig_playGroundVertices;
    Eigen::Matrix3f 				eig_matrixR;
    Eigen::MatrixXf					eig_vectT;    

        eig_matrixR 					= Conv_Vector9fEigen3f(camIFs[i].matrixR);
		eig_vectT 						= Conv_Vector1DEigenVec(camIFs[i].vectT).replicate(1, playGroundVerticesSize);
        eig_playGroundVertices          = ConvTrans_Vector2DEigen(playGroundVertices);
        printEigen2D(eig_matrixR);
        printEigen2D(eig_vectT);
        printEigen2D(eig_playGroundVertices);        
    
    //vecPointTransformed
    Eigen::Matrix<float,3,playGroundVerticesSize,Eigen::RowMajor>  	eig_vecWorldPosRotate;
    Eigen::Matrix<float,3,playGroundVerticesSize,Eigen::RowMajor> 	eig_vecPointTransformed;
			
			eig_vecWorldPosRotate.row(0) 	= eig_playGroundVertices.row(2);
			eig_vecWorldPosRotate.row(1) 	= eig_playGroundVertices.row(0);
			eig_vecWorldPosRotate.row(2) 	= eig_playGroundVertices.row(1);			
			eig_vecPointTransformed 		= eig_matrixR*eig_vecWorldPosRotate + eig_vectT;
            printEigen2D(eig_vecWorldPosRotate);
            printEigen2D(eig_vecPointTransformed);
    
    //lenVec2D, lenVec3D
    Eigen::Array<float, 1, playGroundVerticesSize> 					eig_array_lenVec2D, eig_array_lenVec3D;
    Eigen::Array<float, 3, playGroundVerticesSize> 					eig_array_vecPointTransformed_p2;			
			eig_array_vecPointTransformed_p2		= eig_vecPointTransformed.array().square(); // or .pow(2);			
					                                // Block of size (p,q), starting at (i,j)	: matrix.block(i,j,p,q);			
			eig_array_lenVec2D 						= eig_array_vecPointTransformed_p2.topRows(2).colwise().sum().sqrt();
			eig_array_lenVec3D 						= eig_array_vecPointTransformed_p2.topRows(3).colwise().sum().sqrt();
			
    Eigen::Matrix<float,1,playGroundVerticesSize> 	eig_lenVec2D, eig_lenVec3D;
			eig_lenVec2D 							= eig_array_lenVec2D.matrix();
			eig_lenVec3D 							= eig_array_lenVec3D.matrix();
            printEigen2D(eig_lenVec2D);
            printEigen2D(eig_lenVec3D);

    // //xd, yd
			//masking : if (vecPointTransformed[2] > ImageOffset[i])
			//masking : if (lenVec2D > 1e-6f)
			Eigen::MatrixXf					eig_mask_vecPointTransformed;
			Eigen::MatrixXf					eig_mask_lenVec2D;
				
			eig_mask_vecPointTransformed	= (eig_vecPointTransformed.array().row(2) > ImageOffset[i]).cast<float>();
			eig_mask_lenVec2D				= (eig_lenVec2D.array() > 1e-6f).cast<float>();
                printEigen2D(eig_mask_vecPointTransformed);
                printEigen2D(eig_mask_lenVec2D);		
			//
			Eigen::Array<float,1,playGroundVerticesSize> 	eig_ratio;
			Eigen::Array<float,1,playGroundVerticesSize> 	eig_theta, eig_theta_p2, eig_theta_p4, eig_theta_p6, eig_theta_p8;
			Eigen::Array<float,1,playGroundVerticesSize> 	eig_cdist;
			Eigen::Array<float,1,playGroundVerticesSize> 	eig_xd, eig_yd;
			
			//eig_ratio = (eig_lenVec2D.array() / eig_lenVec3D.array()).matrix();
			eig_ratio 		= eig_lenVec2D.array()*eig_lenVec3D.cwiseInverse().array();
			eig_theta 		= eig_ratio.asin();
			eig_theta_p2 	= eig_theta.square();
			eig_theta_p4 	= eig_theta_p2.square();
			eig_theta_p6 	= eig_theta_p2.cube();
			eig_theta_p8 	= eig_theta_p4.square();
			
			eig_cdist 		= eig_theta*(1.0 + 
										camIFs[i].matrixD[0] * eig_theta_p2 + 
										camIFs[i].matrixD[1] * eig_theta_p4 +
                                        camIFs[i].matrixD[2] * eig_theta_p6 + 
                                        camIFs[i].matrixD[3] * eig_theta_p8);
                printEigen2D(eig_ratio);
                printEigen2D(eig_theta);
                printEigen2D(eig_theta_p2);
                printEigen2D(eig_theta_p4);
                printEigen2D(eig_theta_p6);
                printEigen2D(eig_theta_p8);
                printEigen2D(eig_cdist);
	
			eig_xd = (eig_vecPointTransformed.row(0).array()*eig_lenVec2D.cwiseInverse().array())*eig_cdist;
			eig_yd = (eig_vecPointTransformed.row(1).array()*eig_lenVec2D.cwiseInverse().array())*eig_cdist;

                printEigen2D(eig_xd);
                printEigen2D(eig_yd);

    // vecPoint2D
            Eigen::Array<float,2,playGroundVerticesSize> 	eig_vecPoint2D;	
			eig_vecPoint2D.row(0) = (camIFs[i].matrixK[0] * eig_xd * eig_mask_lenVec2D.array()) + camIFs[i].matrixK[1];
			eig_vecPoint2D.row(1) = (camIFs[i].matrixK[2] * eig_yd * eig_mask_lenVec2D.array()) + camIFs[i].matrixK[3];

            //masking : if ((vecPoint2D[1] >= 0.0 && vecPoint2D[1] <= PLAYGROUND_TEXTURE_HEIGHT) && (vecPoint2D[0] >= 0.0 && vecPoint2D[0] <= PLAYGROUND_TEXTURE_WIDTH))			
            Eigen::MatrixXf					eig_mask_vecPoint2D, eig_maskNOT_vecPoint2D;
			eig_mask_vecPoint2D     = (((eig_vecPoint2D.row(1) >= 0.0) && (eig_vecPoint2D.row(1) <= PLAYGROUND_TEXTURE_HEIGHT)) && ((eig_vecPoint2D.row(0) >= 0.0) && (eig_vecPoint2D.row(0) <= PLAYGROUND_TEXTURE_WIDTH))).cast<float>();
			eig_maskNOT_vecPoint2D  = (!(((eig_vecPoint2D.row(1) >= 0.0) && (eig_vecPoint2D.row(1) <= PLAYGROUND_TEXTURE_HEIGHT)) && ((eig_vecPoint2D.row(0) >= 0.0) && (eig_vecPoint2D.row(0) <= PLAYGROUND_TEXTURE_WIDTH)))).cast<float>();

			//masking : if (vecPointTransformed[2] > ImageOffset[i] + 0.1)
			Eigen::MatrixXf					eig_mask_vecPointTransformed_0p1;				
			eig_mask_vecPointTransformed_0p1	= (eig_vecPointTransformed.array().row(2) > (ImageOffset[i] + 0.1)).cast<float>();
                    printEigen2D(eig_vecPoint2D);
                    cout << "eig_mask_vecPoint2D" << endl;
                    printEigen2D(eig_mask_vecPoint2D);
                    cout << "eig_maskNOT_vecPoint2D" << endl;
                    printEigen2D(eig_maskNOT_vecPoint2D);
                    cout << "eig_mask_vecPointTransformed_0p1" << endl;
                    printEigen2D(eig_mask_vecPointTransformed_0p1);

        //UVA
			Eigen::Array<float,playGroundVerticesSize, 3, Eigen::RowMajor>	eig_UVA;
            Eigen::ArrayXf					eig_UVA_default(3,1);
            // eig_UVA_default <<  -1.0,
            //                     -1.0,
            //                      0.0;       //eig_UVA = eig_UVA_default.replicate(1, playGroundVerticesSize);
			
			eig_UVA.col(0) 					= (eig_vecPoint2D.row(0) / RAW_IMG_WIDTH)*eig_mask_vecPoint2D.array() + (-1.0*eig_maskNOT_vecPoint2D.array());
			eig_UVA.col(1) 					= (eig_vecPoint2D.row(1) / RAW_IMG_HEIGHT)*eig_mask_vecPoint2D.array() + (-1.0*eig_maskNOT_vecPoint2D.array());
			eig_UVA.col(2) 					= (eig_mask_vecPointTransformed_0p1.row(0).array() + 1.0)*eig_mask_vecPoint2D.array();
			
            printEigen2D(eig_UVA);

		for (int j = 0; j < playGroundVerticesSize; ++j ){
            // std::vector<float> v3 (eig_UVA.row(j).data(),  eig_UVA.row(j).data() + eig_UVA.row(j).size() );
            std::vector<float> v3 (&eig_UVA.row(j)(0), &eig_UVA.row(j)(0) + eig_UVA.row(j).cols()*eig_UVA.row(j).rows());
            // std::vector<float> v3 = {eig_UVA.row(j)(0), eig_UVA.row(j)(1), eig_UVA.row(j)(2)};
            uvaOuts.push_back(v3);
		}

        int j = 0;
        cout << "mem adr " << &eig_UVA.row(j)(0) << " mem adr " <<  &eig_UVA.row(j)(1) <<  " mem adr " << &eig_UVA.row(j)(2) << endl;
        cout << "mem adr " << &eig_UVA.row(j)(0) << " mem adr " << &eig_UVA.row(j)(0) + eig_UVA.row(j).cols()*eig_UVA.row(j).rows() << endl;

        // check uvaOuts
        cout << "uvaOuts" << endl;
		for (int j = 0; j < uvaOuts.size(); j++ ){
            for (int k = 0; k < uvaOuts[j].size(); k++ ){
                cout << uvaOuts[j][k] << " ";
            }
            cout << endl;
		}

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

MatrixXf Conv_Vector1DEigenVec(const vector<float> &vector1D) {
    // vector2D[i] >> eigen2D(i, j)
        // convert 9 element vector 1D into 3x3 eigen, row-major
            // eigen3f  = Eigen::Map<Eigen::Matrix<float, 3, 3>> (&vector9f);
            // eigen3f  = Eigen::MatrixXf::Map(&vector9f, 3, 3);

    // MatrixXf eigenVec;
    // eigen1D = Eigen::VectorXf::Map(&vector1D[0], vector1D.size());
    // eigenVec = Eigen::Map<Eigen::VectorXf>(&vector1D[0], vector1D.size());

        // float* ptr_data = &v1[0];
        // Eigen::VectorXf v2 = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(v1.data(), v1.size());

        // Eigen::VectorXf v4 = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(ptr_data, v1.size());
    Eigen::MatrixXf eigenVec    = Eigen::VectorXf::Map(&vector1D[0], vector1D.size());

    return eigenVec;
}