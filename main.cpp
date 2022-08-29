// #include <Eigen/Core>
// #include <Eigen/Geometry>
// #include <Eigen/Dense>

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/Dense"
#include <vector>
#include <iostream>
#include "glm.hpp"      // for vertices glm::vec3
#include <chrono>       // for timer
#include <ctime>        // for random by time
#include <stdlib.h>     // for random by time
#include <cmath>        // for abs()

// integrated by ThanhPhan and debugged
// 18 aug 2022
// optimize mask by .select() function
// prevent crash with 1/0 at lenvec2D/lenvec3D

using namespace std;
using Eigen::MatrixXf;
using Eigen::VectorXf;

MatrixXf    Conv_Vertices2Eigen(const vector<glm::vec3> &vertices);
MatrixXf    Conv_Vector9x1fEigen3x3f(const vector<float> &vector9f);
MatrixXf    Conv_Vector1DEigenVec(const vector<float> &vector1D);
void        printEigen2D(const Eigen::MatrixXf &eigenMat);

//history : redo the GLM processing
struct camIFs_strct {
    vector<float> matrixD       = {1,1,0,0};
    vector<float> matrixK       = {1,1,1,1};
    vector<float> matrixR       = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
    vector<float> vectT         = {1,2,3};
};

int main () {
    //---------------------------------------------
    const int addedVerticesSize       = 1000000;
    const int playGroundVerticesSize  = 5 + addedVerticesSize;
    vector<glm::vec3>       playGroundVertices{
                                    glm::vec3(100,200,300),
                                    glm::vec3(30,4,5),
                                    glm::vec3(5,6,7),
                                    glm::vec3(-7,8,-9),
                                    glm::vec3(90,10,110)
                                };

        // add random point into vertices
        float min = -10;
        float max = 10;
        float randX, randY, randZ;
        //
        for (int count{0}; count < addedVerticesSize; ++count) {
            // int r = minN + rand() % (maxN + 1 - minN)
            // srand((int)time(0)); //!not worked
            // srand(time(NULL)); !not worked
            randX =  min + rand()/(float)RAND_MAX*(max - min);
            randY =  min + rand()/(float)RAND_MAX*(max - min);
            randZ =  min + rand()/(float)RAND_MAX*(max - min);

            // cout << randX << " " << randY << " " << randZ << endl;
            playGroundVertices.push_back(glm::vec3(randX,randY,randZ));            
        }

    int i   =   0;
    vector<camIFs_strct> camIFs;
    camIFs.push_back(camIFs_strct());
    camIFs.push_back(camIFs_strct());
    vector<float> ImageOffset       = {0.1, 0.1};
    float PLAYGROUND_TEXTURE_HEIGHT = 10;
    float PLAYGROUND_TEXTURE_WIDTH  = 10;
    float RAW_IMG_HEIGHT            = 1;
    float RAW_IMG_WIDTH             = 1;
    
    vector<glm::vec3> uvaOuts;
    vector<glm::vec3> uvaOuts_forLoop;

            // // floating point error showcase
            // float a, b, c_tol;
            // a       = 1.0f;
            // b       = 5.999;
            // c_tol   = 1e-20; // tolerance

            // if (abs(a*a*a*a*a*b*a*a*a*a*a*a*a*a*a - b) > c_tol) {
            //     cout << "x*1.0 != x if x is float";
            // }

            // // .inverse() of zero, singularity showcase >> nan
            // Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> zeroMat;
            // zeroMat         = Eigen::MatrixXf::Constant(1, 10, 0.0).array();
            // zeroMat         = zeroMat*zeroMat.inverse() + 1.0;
            // for (int ij{0}; ij < 5; ++ij) {
            //     cout << zeroMat(0,ij);
            // }
            // cout << endl;

            // conditional selection with Eigen
            Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> zeroMat;
            zeroMat         = Eigen::MatrixXf::Zero(1, 10).array();
            zeroMat(0, 1)   = 5.0;

            Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> inverseMat;
            inverseMat      = zeroMat.inverse();

            Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> maskMat;
            // maskMat         = !(inverseMat.isFinite()); // is inf
            maskMat         = !(Eigen::isfinite(inverseMat)); // is inf
                    // isinf(), isnan(), ..

            Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> constMat;
            constMat         = Eigen::MatrixXf::Constant(1, 10, 0.0).array();

            Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> inverseMatFixed;
            //template argument
            // inverseMatFixed  = Eigen::Select<Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> , Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>,Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(maskMat, constMat, inverseMat);
            inverseMatFixed  = maskMat.select(constMat, inverseMat);

            cout << "zeroMat" << endl;
            for (int ij{0}; ij < 10; ++ij) {
                cout << zeroMat(0,ij) << " ";
            }
            cout << endl;

            cout << "inverseMat" << endl;
            for (int ij{0}; ij < 10; ++ij) {
                cout << inverseMat(0,ij) << " ";
            }
            cout << endl;

            cout << "maskMat" << endl;
            for (int ij{0}; ij < 10; ++ij) {
                cout << maskMat(0,ij) << " ";
            }
            cout << endl;           

            cout << "constMat" << endl;
            for (int ij{0}; ij < 10; ++ij) {
                cout << constMat(0,ij) << " ";
            }
            cout << endl;    

            cout << "inverseMatFixed" << endl;
            for (int ij{0}; ij < 10; ++ij) {
                cout << inverseMatFixed(0,ij) << " ";
            }
            cout << endl;                                         

            // fixing nan
            // >> use mask array then select()


    //---------------------------------------------
    auto start_timer        = std::chrono::high_resolution_clock::now();
    auto stop_timer         = std::chrono::high_resolution_clock::now();
    auto duration_time      = std::chrono::duration_cast<std::chrono::milliseconds>(stop_timer - start_timer); // std::chrono::duration<float>

    start_timer             = std::chrono::high_resolution_clock::now();

    //** playGroundVertices[i] new 29 aug 2022
    //** playGroundVertices iteration----------------vectorization---------------------------	
	 	Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_playGroundVertices;
        Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_matrixR;
            eig_matrixR 				    = Conv_Vector9x1fEigen3x3f(camIFs[i].matrixR).array();
            eig_playGroundVertices          = Conv_Vertices2Eigen(playGroundVertices).array();
		
        // // // vecPointTransformed
        Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_vecWorldPosRotate;
        Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_vecPointTransformed;
            eig_vecWorldPosRotate           = Eigen::MatrixXf::Constant(3, playGroundVerticesSize, 0.0).array();
            eig_vecWorldPosRotate.row(0) 	= eig_playGroundVertices.row(2);
            eig_vecWorldPosRotate.row(1) 	= eig_playGroundVertices.row(0);
            eig_vecWorldPosRotate.row(2) 	= eig_playGroundVertices.row(1);
            
            // rotate
            eig_vecPointTransformed 		= (eig_matrixR.matrix()*eig_vecWorldPosRotate.matrix()).array();
            
            // translation
            eig_vecPointTransformed.row(0)  += camIFs[i].vectT[0];
            eig_vecPointTransformed.row(1)  += camIFs[i].vectT[1];
            eig_vecPointTransformed.row(2)  += camIFs[i].vectT[2];

	    // // // lenVec2D, lenVec3D	
        Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_lenVec2D, eig_lenVec3D;
		Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_vecPointTransformed_p2;
            eig_vecPointTransformed_p2      = eig_vecPointTransformed.square();
            eig_lenVec2D 				    = eig_vecPointTransformed_p2.topRows(2).colwise().sum().sqrt();
            eig_lenVec3D 					= eig_vecPointTransformed_p2.topRows(3).colwise().sum().sqrt();

		// // // xd, yd
        // masking : if (vecPointTransformed[2] > ImageOffset[i])
        // masking : if (lenVec2D > 1e-6f)        
        Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_mask_vecPointTransformed;
        Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_maskNOT_vecPointTransformed;
        Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_mask_lenVec2D;
            eig_mask_vecPointTransformed	= (eig_vecPointTransformed.row(2) > ImageOffset[i]).cast<float>();
            eig_maskNOT_vecPointTransformed	= (!(eig_vecPointTransformed.row(2) > ImageOffset[i])).cast<float>();
            eig_mask_lenVec2D				= (eig_lenVec2D > 1e-6f).cast<float>();
	
        // // // ratio, theta
        Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_ratio;
        Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_theta, eig_theta_p2, eig_theta_p4, eig_theta_p6, eig_theta_p8;
        Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_cdist;
        Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_xd, eig_yd;
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
            // // // xd, yd
            eig_xd          = (eig_vecPointTransformed.row(0).array()*eig_lenVec2D.inverse())*eig_cdist;
            eig_yd          = (eig_vecPointTransformed.row(1).array()*eig_lenVec2D.inverse())*eig_cdist;

    //         // printEigen2D(eig_xd);
    //         // printEigen2D(eig_yd);

	// 	// // // vecPoint2D
        Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_vecPoint2D;
            eig_vecPoint2D        = Eigen::MatrixXf::Constant(2, playGroundVerticesSize, 0.0).array();
            eig_vecPoint2D.row(0)   = (camIFs[i].matrixK[0] * eig_xd * eig_mask_lenVec2D) + camIFs[i].matrixK[1];
            eig_vecPoint2D.row(1)   = (camIFs[i].matrixK[2] * eig_yd * eig_mask_lenVec2D) + camIFs[i].matrixK[3];


            // ***************
            //masking : if ((vecPoint2D[1] >= 0.0 && vecPoint2D[1] <= PLAYGROUND_TEXTURE_HEIGHT) && (vecPoint2D[0] >= 0.0 && vecPoint2D[0] <= PLAYGROUND_TEXTURE_WIDTH))
            // Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_mask_vecPoint2D, eig_maskNOT_vecPoint2D;
            // eig_mask_vecPoint2D     = (((eig_vecPoint2D.row(1) >= 0.0) && (eig_vecPoint2D.row(1) <= PLAYGROUND_TEXTURE_HEIGHT)) && ((eig_vecPoint2D.row(0) >= 0.0) && (eig_vecPoint2D.row(0) <= PLAYGROUND_TEXTURE_WIDTH))).cast<float>();
            // eig_maskNOT_vecPoint2D  = (!(((eig_vecPoint2D.row(1) >= 0.0) && (eig_vecPoint2D.row(1) <= PLAYGROUND_TEXTURE_HEIGHT)) && ((eig_vecPoint2D.row(0) >= 0.0) && (eig_vecPoint2D.row(0) <= PLAYGROUND_TEXTURE_WIDTH)))).cast<float>();

            // ***************
            Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_mask_vecPoint2D;
            eig_mask_vecPoint2D     = (((eig_vecPoint2D.row(1) >= 0.0) && (eig_vecPoint2D.row(1) <= PLAYGROUND_TEXTURE_HEIGHT)) && ((eig_vecPoint2D.row(0) >= 0.0) && (eig_vecPoint2D.row(0) <= PLAYGROUND_TEXTURE_WIDTH)));
            // ***************

            //masking : if (vecPointTransformed[2] > ImageOffset[i] + 0.1)

            // // **********************
            // Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_mask_vecPointTransformed_0p1;
            // eig_mask_vecPointTransformed_0p1 = (eig_vecPointTransformed.array().row(2) > (ImageOffset[i] + 0.1)).cast<float>();
            // // **********************
            
            Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> eig_mask_vecPointTransformed_offset0p1;
            eig_mask_vecPointTransformed_offset0p1 = (eig_vecPointTransformed.array().row(2) > (ImageOffset[i] + 0.1));

    //             // printEigen2D(eig_vecPoint2D.matrix());
    //             // printEigen2D(eig_mask_vecPoint2D.matrix());
    //             // printEigen2D(eig_maskNOT_vecPoint2D.matrix());
    //             // cout << "eig_mask_vecPointTransformed_0p1" << endl;
    //             // printEigen2D(eig_mask_vecPointTransformed_0p1.matrix());

	// 		//UVA
	// 			// keep to 3xsize matrice
                Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic> eig_rowOnes, eig_rowMinusOnes, eig_rowZeros;
                    eig_rowOnes         = Eigen::MatrixXf::Constant(1, playGroundVerticesSize, 1.0).array();
                    eig_rowMinusOnes    = Eigen::MatrixXf::Constant(1, playGroundVerticesSize, -1.0).array();
                    eig_rowZeros        = Eigen::MatrixXf::Constant(1, playGroundVerticesSize, 0.0).array();

                // ****************
                Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic> eig_UVA;
                eig_UVA         = Eigen::MatrixXf::Constant(3, playGroundVerticesSize, 0.0).array();

                eig_UVA.row(0) 	= eig_mask_vecPoint2D.select(eig_vecPoint2D.row(0) / RAW_IMG_WIDTH, eig_rowMinusOnes);
                eig_UVA.row(1) 	= eig_mask_vecPoint2D.select(eig_vecPoint2D.row(1) / RAW_IMG_HEIGHT, eig_rowMinusOnes);
                eig_UVA.row(2) 	= eig_mask_vecPoint2D.select(eig_mask_vecPointTransformed_offset0p1.select(eig_rowOnes, eig_rowZeros), eig_rowZeros);
                // ****************


				// eig_UVA.row(0) 	= (eig_vecPoint2D.row(0) / RAW_IMG_WIDTH)*eig_mask_vecPoint2D + (-1.0*eig_maskNOT_vecPoint2D);
				// eig_UVA.row(1) 	= (eig_vecPoint2D.row(1) / RAW_IMG_HEIGHT)*eig_mask_vecPoint2D + (-1.0*eig_maskNOT_vecPoint2D);
				// eig_UVA.row(2) 	= (eig_mask_vecPointTransformed_0p1*1.0)*eig_mask_vecPoint2D;
				
                // cout << "eig_UVA" << endl;
                // printEigen2D(eig_UVA.matrix());

				eig_UVA.row(0) 	= eig_UVA.row(0)*eig_mask_vecPointTransformed + (-1.0*eig_maskNOT_vecPointTransformed);
				eig_UVA.row(1) 	= eig_UVA.row(1)*eig_mask_vecPointTransformed + (-1.0*eig_maskNOT_vecPointTransformed);
				eig_UVA.row(2) 	= eig_UVA.row(2)*eig_mask_vecPointTransformed;

    //             // printEigen2D(eig_UVA.matrix());

			// push_back to uvaOuts
				for (int k = 0; k < playGroundVerticesSize; ++k ){
					uvaOuts.push_back(glm::vec3(eig_UVA(0, k), eig_UVA(1, k), eig_UVA(2, k)));
				}

    // new 29 aug 2022
	// 	//** playGroundVertices iteration----------------vectorization---------------------------
	// 	//-----------vectorization ^^^^^^^^^^^^^^^^^^^^^^^^^^
    // //** playGroundVertices iteration----------------vectorization---------------------------

    stop_timer 	        = std::chrono::high_resolution_clock::now();
    duration_time 	    = std::chrono::duration_cast<std::chrono::milliseconds>(stop_timer - start_timer);
    std::cout << "elapsed time vetorization : " << duration_time.count() << " [ms] " << std::endl;

    start_timer         = std::chrono::high_resolution_clock::now();
        //------------------------------------ for-loop

        for (int j = 0; j < playGroundVerticesSize; j++ ){
            glm::vec3 UVA = glm::vec3(-1.0, -1.0, 0.0);
            glm::vec3 vecVerticeWorld = playGroundVertices[j];
            std::vector<float> vecWorldPosRotate;
            vecWorldPosRotate.push_back(vecVerticeWorld[2]);
            vecWorldPosRotate.push_back(vecVerticeWorld[0]);
            vecWorldPosRotate.push_back(vecVerticeWorld[1]);

            std::vector<float> vecPointTransformed;
            vecPointTransformed.push_back(camIFs[i].matrixR[0] * vecWorldPosRotate[0] + camIFs[i].matrixR[1] * vecWorldPosRotate[1] + camIFs[i].matrixR[2] * vecWorldPosRotate[2] + camIFs[i].vectT[0]);
            vecPointTransformed.push_back(camIFs[i].matrixR[3] * vecWorldPosRotate[0] + camIFs[i].matrixR[4] * vecWorldPosRotate[1] + camIFs[i].matrixR[5] * vecWorldPosRotate[2] + camIFs[i].vectT[1]);
            vecPointTransformed.push_back(camIFs[i].matrixR[6] * vecWorldPosRotate[0] + camIFs[i].matrixR[7] * vecWorldPosRotate[1] + camIFs[i].matrixR[8] * vecWorldPosRotate[2] + camIFs[i].vectT[2]);

            float lenVec3D = sqrt(
                vecPointTransformed[0] * vecPointTransformed[0] + vecPointTransformed[1] * vecPointTransformed[1] +
                vecPointTransformed[2] * vecPointTransformed[2]);
            float lenVec2D = sqrt(
                vecPointTransformed[0] * vecPointTransformed[0] + vecPointTransformed[1] * vecPointTransformed[1]);

            if (vecPointTransformed[2] > ImageOffset[i]){
                float ratio = 0.0;
                float xd, yd;
                float theta = 0.0f;
                std::vector<float> vecPoint2D(2, 0.0);
                if (lenVec2D > 1e-6f){
                    ratio = lenVec2D / lenVec3D;
                    theta = asin(ratio);
                    float thetaSq = theta * theta;

                    float cdist = theta * (1.0 + camIFs[i].matrixD[0] * thetaSq + camIFs[i].matrixD[1] * thetaSq * thetaSq +
                                        camIFs[i].matrixD[2] * thetaSq * thetaSq * thetaSq +
                                        camIFs[i].matrixD[3] * thetaSq * thetaSq * thetaSq * thetaSq);

                    xd = vecPointTransformed[0] / lenVec2D * cdist;
                    yd = vecPointTransformed[1] / lenVec2D * cdist;

                    vecPoint2D[0] = camIFs[i].matrixK[0] * xd + camIFs[i].matrixK[1];
                    vecPoint2D[1] = camIFs[i].matrixK[2] * yd + camIFs[i].matrixK[3];
                }else{
                    vecPoint2D[0] = camIFs[i].matrixK[1];
                    vecPoint2D[1] = camIFs[i].matrixK[3];
                }

                if ((vecPoint2D[1] >= 0.0 && vecPoint2D[1] <= PLAYGROUND_TEXTURE_HEIGHT) && (vecPoint2D[0] >= 0.0 && vecPoint2D[0] <= PLAYGROUND_TEXTURE_WIDTH)){
                    if (vecPointTransformed[2] > ImageOffset[i] + 0.1){
                        UVA = glm::vec3(vecPoint2D[0] / RAW_IMG_WIDTH, vecPoint2D[1] / RAW_IMG_HEIGHT, 1.0);
                    }else{
                        UVA = glm::vec3(vecPoint2D[0] / RAW_IMG_WIDTH, vecPoint2D[1] / RAW_IMG_HEIGHT, 0.0);
                    }
                }else{
                    UVA = glm::vec3(-1.0, -1.0, 0.0);	
                }
            }
            uvaOuts_forLoop.push_back(UVA);
        }

    //     //------------------------------------ for-loop ^^^^^^^^^^^^^^^^^^^^^^^

    stop_timer 	        = std::chrono::high_resolution_clock::now();
    duration_time 	    = std::chrono::duration_cast<std::chrono::milliseconds>(stop_timer - start_timer);
    std::cout << "elapsed time for-loop : " << duration_time.count() << " [ms] " << std::endl;

    //     // check uvaOuts
    //     // cout << "uvaOuts" << endl;
	// 	// for (int j = 0; j < uvaOuts.size(); j++ ){
    //     //     for (int k = 0; k < 3; k++ ){
    //     //         cout << uvaOuts[j][k] << " ";
    //     //     }
    //     //     cout << endl;
	// 	// }

    //     // cout << "uvaOuts_forLoop" << endl;
	// 	// for (int j = 0; j < uvaOuts_forLoop.size(); j++ ){
    //     //     for (int k = 0; k < 3; k++ ){
    //     //         cout << uvaOuts_forLoop[j][k] << " ";
    //     //     }
    //     //     cout << endl;
	// 	// }

        cout << "uvaOuts - uvaOuts_forLoop :: comparing result " << endl;
        bool    different   = false;
        int     difCount    = 0;
        float   tol         = 1e-4;
		for (int j = 0; j < uvaOuts.size(); j++ ){
            for (int k = 0; k < 3; k++ ){
                if (abs(uvaOuts[j][k] - uvaOuts_forLoop[j][k]) > tol) {
                    different = true;
                }
                
            }

            if (different) {
                ++difCount;
                cout << "aie aie aie aie aie aie aie aie aie aie aie aie" << endl;
                cout << "j " << j << endl;
                cout << "uvaOuts         " << uvaOuts[j][0] << " " << uvaOuts[j][1] << " " <<uvaOuts[j][2] << " " <<endl;
                cout << "uvaOuts_forLoop " << uvaOuts_forLoop[j][0] << " " << uvaOuts_forLoop[j][1] << " " << uvaOuts_forLoop[j][2] << " " << endl;
            }

            different = false;
            // cout << endl;
		}
        if (difCount == 0) {
            cout << "PASSED!        no different at tolerance " << tol << endl;
        } else {
            cout << "NOT passed!!!  different at tolerance "    << tol << endl;
            cout << difCount << endl;
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
    int verticesSize = vertices.size();
    for (int i=0; i < verticesSize; ++i) {
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
    Eigen::MatrixXf eigenVec    = Eigen::VectorXf::Map(&vector1D[0], 3);
    return eigenVec;
}