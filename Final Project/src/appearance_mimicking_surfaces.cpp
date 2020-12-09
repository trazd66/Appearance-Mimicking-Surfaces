// Given a mesh, a viewpoint, and linear constraints:
//      lambdaMin -- combined with lambdaMax this will "squeeze" the mesh into the defined thickness by projecting the vertices
//      lambdaMax
//      mu -- different parts of the shape can be given different thickness by segmenting the vertices into groups
//      fixed vertices
//      weights -- scales the difference of vertex normals so that more visible vertices are given preserved "more" than less visible
//      optional depth ordering
// output a bas-relief deformed mesh whose appearance is preserved
// from the viewpoint
//
//
// Inputs:
//    V            #V by 3 list of the vertex positions of the model
//    F            #F by 3 list of triangle indices into V
//    view         3D vector of the coordinates of the viewpoint
//    lambdaMin    #Constraints by 2 list of vertex indices and the lambdaMinValue: [vertexIdx,lambdaMinValue]
//    lambdaMax    #Constraints by 2 list of [vertexIdx,lambdaMaxValue]
//    bf           #fixedVertices vector of indices into V
//    weights      #V length list of vertex weights
//    mu           #V length list of mu indices
//                      (separates mesh into independent regions each with their own thickness constraint)
//                      (default should be all 1's)
//  Outputs:
//    DV           #V by 3 list of the vertex positions of the deformed model

#include "appearance_mimicking_surfaces.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/repdiag.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/active_set.h>
#include <iostream>

/*
 * Given a mesh, a viewpoint and linear constraints, outputs a bas-relief deformed mesh whose appearance is preserved.
 */
void appearance_mimicking_surfaces(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const Eigen::Vector3d &view,
        const Eigen::MatrixXd &lambdaMin,
		const Eigen::MatrixXd &lambdaMax,
        const Eigen::VectorXi &bf,
        const Eigen::VectorXd &weights,
        const Eigen::VectorXd &mu,
        Eigen::MatrixXd &DV) {

	int num_vertices = V.rows();
	int mu_len = mu.size();
	DV.resize(num_vertices, 3);

	// TODO:  We need to the "Voronoi area" of each vertex.  Do we call igl::massmatrix() or igl::point_areas()?
	// I spent too much time looking into point_areas() I think massmatrix() is actually the right way to get the
	// Voronoi areas.

	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    M = M.cwiseSqrt();

	// D_A(i, j) = sqrt(A_i^0) if i = j else 0.
	Eigen::SparseMatrix<double> D_A(3 * num_vertices, 3 * num_vertices);
    igl::repdiag(M,3,D_A);

	// D_w(i, j) = w_i if i = j else 0.
	Eigen::SparseMatrix<double> D_w(3 * num_vertices, 3 * num_vertices);
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> w_matrix;
    w_matrix.diagonal() = weights;
    igl::repdiag((Eigen::SparseMatrix<double>)w_matrix,3,D_w);

	// D_V_tilde = Diagonal matrix containing row-wise stacked elements of V.
	Eigen::DiagonalMatrix<double,Eigen::Dynamic> D_V_tilde;
    Eigen::MatrixXd V_tilde(V.rows(),3);

    V_tilde = V.rowwise() - view.transpose();
    D_V_tilde.diagonal() = Eigen::Map<const Eigen::VectorXd,Eigen::RowMajor>(V_tilde.data(),V_tilde.size());

    //Lambda_0 is the per vertex lambda value (||v||) of the undeformed mesh
    Eigen::VectorXd Lambda_0(num_vertices);
    Lambda_0 = V_tilde.rowwise().norm();

	// Cotangent Laplace-Beltrami operator.
	Eigen::SparseMatrix<double> L(num_vertices, num_vertices);
	igl::cotmatrix(V, F, L);

	// L^0~ = L^0 (x) I_3, where (x) is the Kronecker product.
	Eigen::SparseMatrix<double> L_0_tilde(3 * num_vertices, 3 * num_vertices);

	// Three entries per non-zero coefficient in L.
	std::vector<Eigen::Triplet<double>>L_triplets;
	L_triplets.reserve(3 * L.nonZeros());

	for (int row = 0; row < L.outerSize(); row++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(L, row); it; ++it) {
			L_triplets.push_back(Eigen::Triplet<double>(3 * it.row(), 3 * it.col(), it.value()));
			L_triplets.push_back(Eigen::Triplet<double>(3 * it.row() + 1, 3 * it.col() + 1, it.value()));
			L_triplets.push_back(Eigen::Triplet<double>(3 * it.row() + 2, 3 * it.col() + 2, it.value()));
		}
	}
	L_0_tilde.setFromTriplets(L_triplets.begin(), L_triplets.end());

	// S = selector matrix = I_n (x) [1, 1, 1]^T.
	Eigen::SparseMatrix<double> S(3 * num_vertices, num_vertices);
	std::vector<Eigen::Triplet<int>> S_triplets;
	S_triplets.reserve(3 * num_vertices);
	for (int vertex = 0; vertex < num_vertices; vertex++) {
		for (int i = 0; i < 3; i++) {
			S_triplets.push_back(Eigen::Triplet<int>(3 * vertex + 1, vertex, 1));
		}
	}
	S.setFromTriplets(S_triplets.begin(), S_triplets.end());

    //intermediate step to make compiler happy
    Eigen::VectorXd S_Lambda_0,L_theta;
    S_Lambda_0 = S * Lambda_0;

    L_theta = S_Lambda_0.asDiagonal().inverse() * L_0_tilde * D_V_tilde * S_Lambda_0;


    //D_L_theta : 3n x 3n sparsematrix
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> D_L_theta;
    D_L_theta.diagonal() = L_theta;

    Eigen::MatrixXd F_top_left = D_A * D_w * (L_0_tilde * D_V_tilde - (Eigen::SparseMatrix<double>) D_L_theta) * S;
    std::cout << F_top_left.rows() << ", " << F_top_left.cols() << std::endl;
    Eigen::MatrixXd F(3 * num_vertices + mu_len, num_vertices + mu_len);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(n);

    // TODO:  Find an optimal solution for lambda.
	Eigen::VectorXd lambda(num_vertices);

	// Set the vertices of the bas-relief shape.
	for (int i = 0; i < num_vertices; i++) {
		DV.row(i) = lambda(i) * v_hat.row(i);
	}
}

//  if you are feeling fancy, this is how they implemented the depth constraints
//    zOrder       #Constraints x 8 [order,vIdx,vAIdx,vBIdx,vCIdx,u,v,w]
//                 Contains the depth order constraints
//                 for vertex <-> triangle point pairs
//                 order: 1 vertex in front of triangle; -1 behind triangle
//                 vIdx:  index of vertex
//                 vAIdx: index of first vertex (A) of triangle
//                 vBIdx: index of second vertex (B) of triangle
//                 vCIdx: index of third vertex (C) of triangle
//                 u: barycentric coordinate of point of triangle A
//                 v: barycentric coordinate of point of triangle B
//                 w: barycentric coordinate of point of triangle C
