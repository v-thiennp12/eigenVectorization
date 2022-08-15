#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "glm.hpp" // glm::vec3

using namespace std;
using Eigen::MatrixXf;
using Eigen::VectorXf;

MatrixXf    Conv_Vertices2Eigen(const vector<glm::vec3> &vertices);
MatrixXf    Conv_Vector9x1fEigen3x3f(const vector<float> &vector9f);
MatrixXf    Conv_Vector1DEigenVec(const vector<float> &vector1D);
void        printEigen2D(const Eigen::MatrixXf &eigenMat);

//history : redo the GLM processing
struct camIFs_strct {
    vector<float> matrixD       = {1,2,3,4};
    vector<float> matrixK       = {1,2,3,4};
    vector<float> matrixR       = {1,2,3,4,5,6,7,8,9};
    vector<float> vectT         = {1,2,3};
};

int main () {
    //---------------------------------------------
    const int playGroundVerticesSize  = 5;
    vector<glm::vec3>       playGroundVertices{
                                    glm::vec3(1,2,3),
                                    glm::vec3(3,4,5),
                                    glm::vec3(5,6,7),
                                    glm::vec3(7,8,9),
                                    glm::vec3(9,10,11)
                                };
    int i   =   0;
    vector<camIFs_strct> camIFs;
    camIFs.push_back(camIFs_strct());
    camIFs.push_back(camIFs_strct());
    vector<float> ImageOffset       = {0.1, 0.2};
    float PLAYGROUND_TEXTURE_HEIGHT = 5;
    float PLAYGROUND_TEXTURE_WIDTH  = 2.5;
    float RAW_IMG_HEIGHT            = 600;
    float RAW_IMG_WIDTH             = 800;
    
    vector<glm::vec3> uvaOuts;

    //---------------------------------------------

    //** playGroundVertices[i]

    Eigen::MatrixXf 				eig_playGroundVertices;
    Eigen::Matrix3f 				eig_matrixR;
    Eigen::MatrixXf					eig_vectT;    

        eig_matrixR 				= Conv_Vector9x1fEigen3x3f(camIFs[i].matrixR);
		eig_vectT 					= Conv_Vector1DEigenVec(camIFs[i].vectT).replicate(1, playGroundVerticesSize);
        eig_playGroundVertices      = Conv_Vertices2Eigen(playGroundVertices);
            printEigen2D(eig_matrixR);
            cout << "vect_T" << endl;
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
    Eigen::Array<float, 1, playGroundVerticesSize> 					eig_lenVec2D, eig_lenVec3D;
    Eigen::Array<float, 3, playGroundVerticesSize> 					eig_vecPointTransformed_p2;			
			eig_vecPointTransformed_p2          = eig_vecPointTransformed.array().square(); // or .pow(2);			
					                                // Block of size (p,q), starting at (i,j)	: matrix.block(i,j,p,q);			
			eig_lenVec2D 						= eig_vecPointTransformed_p2.topRows(2).colwise().sum().sqrt();
			eig_lenVec3D 						= eig_vecPointTransformed_p2.topRows(3).colwise().sum().sqrt();
			
    // //xd, yd
			//masking : if (vecPointTransformed[2] > ImageOffset[i])
			//masking : if (lenVec2D > 1e-6f)
			Eigen::Array<float, 1, playGroundVerticesSize>			eig_mask_vecPointTransformed;
			Eigen::Array<float, 1, playGroundVerticesSize>			eig_maskNOT_vecPointTransformed;
			Eigen::Array<float, 1, playGroundVerticesSize>			eig_mask_lenVec2D;
				
			eig_mask_vecPointTransformed	= (eig_vecPointTransformed.array().row(2) > ImageOffset[i]).cast<float>();
			eig_maskNOT_vecPointTransformed	= (!(eig_vecPointTransformed.array().row(2) > ImageOffset[i])).cast<float>();
			eig_mask_lenVec2D				= (eig_lenVec2D > 1e-6f).cast<float>();
                cout << "eig_mask_vecPointTransformed" << endl;
                printEigen2D(eig_mask_vecPointTransformed);
                cout << "eig_maskNOT_vecPointTransformed" << endl;
                printEigen2D(eig_maskNOT_vecPointTransformed);
                printEigen2D(eig_mask_lenVec2D);		
			//
			Eigen::Array<float,1,playGroundVerticesSize> 	eig_ratio;
			Eigen::Array<float,1,playGroundVerticesSize> 	eig_theta, eig_theta_p2, eig_theta_p4, eig_theta_p6, eig_theta_p8;
			Eigen::Array<float,1,playGroundVerticesSize> 	eig_cdist;
			Eigen::Array<float,1,playGroundVerticesSize> 	eig_xd, eig_yd;
			
			//eig_ratio = (eig_lenVec2D.array() / eig_lenVec3D.array()).matrix();
			eig_ratio 		= eig_lenVec2D*eig_lenVec3D.inverse();
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
	
			eig_xd = (eig_vecPointTransformed.row(0).array()*eig_lenVec2D.inverse())*eig_cdist;
			eig_yd = (eig_vecPointTransformed.row(1).array()*eig_lenVec2D.inverse())*eig_cdist;

                printEigen2D(eig_xd);
                printEigen2D(eig_yd);

    // vecPoint2D
            Eigen::Array<float,2,playGroundVerticesSize> 	eig_vecPoint2D;	
			eig_vecPoint2D.row(0) = (camIFs[i].matrixK[0] * eig_xd * eig_mask_lenVec2D) + camIFs[i].matrixK[1];
			eig_vecPoint2D.row(1) = (camIFs[i].matrixK[2] * eig_yd * eig_mask_lenVec2D) + camIFs[i].matrixK[3];

            //masking : if ((vecPoint2D[1] >= 0.0 && vecPoint2D[1] <= PLAYGROUND_TEXTURE_HEIGHT) && (vecPoint2D[0] >= 0.0 && vecPoint2D[0] <= PLAYGROUND_TEXTURE_WIDTH))			
            Eigen::Array<float, 1, playGroundVerticesSize>  eig_mask_vecPoint2D;
            Eigen::Array<float, 1, playGroundVerticesSize>  eig_maskNOT_vecPoint2D;
			eig_mask_vecPoint2D     = (((eig_vecPoint2D.row(1) >= 0.0) && (eig_vecPoint2D.row(1) <= PLAYGROUND_TEXTURE_HEIGHT)) && ((eig_vecPoint2D.row(0) >= 0.0) && (eig_vecPoint2D.row(0) <= PLAYGROUND_TEXTURE_WIDTH))).cast<float>();
			eig_maskNOT_vecPoint2D  = (!(((eig_vecPoint2D.row(1) >= 0.0) && (eig_vecPoint2D.row(1) <= PLAYGROUND_TEXTURE_HEIGHT)) && ((eig_vecPoint2D.row(0) >= 0.0) && (eig_vecPoint2D.row(0) <= PLAYGROUND_TEXTURE_WIDTH)))).cast<float>();

			//masking : if (vecPointTransformed[2] > ImageOffset[i] + 0.1)
			Eigen::Array<float, 1, playGroundVerticesSize>	eig_mask_vecPointTransformed_0p1;				
			eig_mask_vecPointTransformed_0p1	= (eig_vecPointTransformed.array().row(2) > (ImageOffset[i] + 0.1)).cast<float>();
                    printEigen2D(eig_vecPoint2D);
                    cout << "eig_mask_vecPoint2D" << endl;
                    printEigen2D(eig_mask_vecPoint2D);
                    cout << "eig_maskNOT_vecPoint2D" << endl;
                    printEigen2D(eig_maskNOT_vecPoint2D);
                    cout << "eig_mask_vecPointTransformed_0p1" << endl;
                    printEigen2D(eig_mask_vecPointTransformed_0p1);

        //UVA
            // transposed to size x 3 matrice as vertices 
			// Eigen::Array<float,playGroundVerticesSize, 3, Eigen::RowMajor>	eig_UVA;
			// eig_UVA.col(0) 					= (eig_vecPoint2D.row(0) / RAW_IMG_WIDTH)*eig_mask_vecPoint2D + (-1.0*eig_maskNOT_vecPoint2D);
			// eig_UVA.col(1) 					= (eig_vecPoint2D.row(1) / RAW_IMG_HEIGHT)*eig_mask_vecPoint2D + (-1.0*eig_maskNOT_vecPoint2D);
			// eig_UVA.col(2) 					= (eig_mask_vecPointTransformed_0p1.row(0) + 1.0)*eig_mask_vecPoint2D;
			
            // eig_UVA.col(0) 					= eig_UVA.col(0)*eig_mask_vecPointTransformed.row(0) + (-1.0*eig_maskNOT_vecPointTransformed.row(0));
            // eig_UVA.col(1) 					= eig_UVA.col(1)*eig_mask_vecPointTransformed.row(0) + (-1.0*eig_maskNOT_vecPointTransformed.row(0));
            // eig_UVA.col(2) 					= eig_UVA.col(2)*eig_mask_vecPointTransformed.row(0);
            // printEigen2D(eig_UVA);

            // keep to 3xsize matrice
			Eigen::Array<float, 3, playGroundVerticesSize, Eigen::RowMajor>	eig_UVA;
			eig_UVA.row(0) 					= (eig_vecPoint2D.row(0) / RAW_IMG_WIDTH)*eig_mask_vecPoint2D + (-1.0*eig_maskNOT_vecPoint2D);
			eig_UVA.row(1) 					= (eig_vecPoint2D.row(1) / RAW_IMG_HEIGHT)*eig_mask_vecPoint2D + (-1.0*eig_maskNOT_vecPoint2D);
			eig_UVA.row(2) 					= (eig_mask_vecPointTransformed_0p1.row(0) + 1.0)*eig_mask_vecPoint2D;
			
            eig_UVA.row(0) 					= eig_UVA.row(0)*eig_mask_vecPointTransformed + (-1.0*eig_maskNOT_vecPointTransformed);
            eig_UVA.row(1) 					= eig_UVA.row(1)*eig_mask_vecPointTransformed + (-1.0*eig_maskNOT_vecPointTransformed);
            eig_UVA.row(2) 					= eig_UVA.row(2)*eig_mask_vecPointTransformed;
            printEigen2D(eig_UVA);

        // push_back to uvaOuts
            //for-loop
            // for (int k = 0; k < playGroundVerticesSize; ++k ){
            //     uvaOuts.push_back(glm::vec3(eig_UVA(k, 0), eig_UVA(k, 1), eig_UVA(k, 2)));
            // }

            for (int k = 0; k < playGroundVerticesSize; ++k ){
                uvaOuts.push_back(glm::vec3(eig_UVA(0, k), eig_UVA(1, k), eig_UVA(2, k)));
            }            
            //at once
            // https://forums.codeguru.com/showthread.php?519570-Add-multiple-values-to-a-vector
            // https://stackoverflow.com/questions/10516495/calling-a-function-on-every-element-of-a-c-vector
            // uvaOuts()
            // array<int, 5> tmp = {4, 57, 5786, 578, 8841}; 
            // vector<int> vec(tmp.begin(), tmp.end());


        int j = 0;
        cout << "mem adr " << &eig_UVA.row(j)(0) << " mem adr " << &eig_UVA.row(j)(1) <<  " mem adr " << &eig_UVA.row(j)(2) << endl;
        cout << "mem adr " << &eig_UVA.row(j)(0) << " mem adr " << &eig_UVA.row(j)(0) + 1 <<  " mem adr " << &eig_UVA.row(j)(0) + 2 << " mem adr " << &eig_UVA.row(j)(0) + 3 << endl;
        cout << "mem adr " << &eig_UVA.row(j)(0) << " mem adr " << &eig_UVA.row(j)(0) + eig_UVA.row(j).cols()*eig_UVA.row(j).rows() << endl;

        // check uvaOuts
        cout << "uvaOuts" << endl;
		for (int j = 0; j < uvaOuts.size(); j++ ){
            for (int k = 0; k < 3; k++ ){
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

MatrixXf Conv_Vertices2Eigen(const vector<glm::vec3> &vertices) {
    // vertices[i][j] >> eigen2D(j, i)
        // convert vector of vertices (vec3) into 3x_ eigen 2D RowMajor

    Eigen::Matrix<float, 3, Eigen::Dynamic, Eigen::RowMajor> eigen2D = Eigen::MatrixXf::Constant(3, vertices.size(), 0.0f);

    for (int i{0}; i < vertices.size(); ++i) {
        eigen2D(0,i)  = vertices[i][0];
        eigen2D(1,i)  = vertices[i][1];
        eigen2D(2,i)  = vertices[i][2];
    }
    return eigen2D;
}

MatrixXf Conv_Vector9x1fEigen3x3f(const vector<float> &vector9f) {
    // vector1D[k] >> eigen2D(i, j)
        // convert 9 element vector 1D into 3x3 eigen, row-major

    MatrixXf eigen3x3f          = Eigen::MatrixXf::Constant(3, 3, 0.0f);
        eigen3x3f.row(0)        = Eigen::VectorXf::Map(&vector9f[0], 3);
        eigen3x3f.row(1)        = Eigen::VectorXf::Map(&vector9f[3], 3);
        eigen3x3f.row(2)        = Eigen::VectorXf::Map(&vector9f[6], 3);

    return eigen3x3f;
}

MatrixXf Conv_Vector1DEigenVec(const vector<float> &vector1D) {
    Eigen::MatrixXf eigenVec    = Eigen::VectorXf::Map(&vector1D[0], vector1D.size());
    return eigenVec;
}