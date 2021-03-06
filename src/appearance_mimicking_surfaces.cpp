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
#include <igl/ray_mesh_intersect.h>

/*
 * Shoot raycasts at each vertex to build up the inequality constraints matrix C_I.
 */
void depth_constraint(const int num_vertices,
					  const Eigen::Vector3d &view,
					  const Eigen::MatrixXd &V,
					  const Eigen::MatrixXi &F,
					  int &curr_C_I_row,
					  std::vector<Eigen::Triplet<double>> &C_I_triplets) {

	for (int i = 0; i < num_vertices; i++) {
   		std::vector<igl::Hit> hits;
   		Eigen::RowVectorXd direction = V.row(i);

   		if (igl::ray_mesh_intersect(view, direction, V, F, hits))
        {
        	// Each raycast hit corresponds to a single row of C_I.
            for (int hit_no = 0; hit_no < hits.size() - 1; hit_no++) {

            	igl::Hit front_hit = hits[hit_no];
                int front_face = front_hit.id;

                // Front face vertices.
                int u0 = F(front_face, 0);
                int u1 = F(front_face, 1);
                int u2 = F(front_face, 2);

                // Front face barycentric coordinates.
                double u1_bary = front_hit.u;
                double u2_bary = front_hit.v;
                double u0_bary = 1 - u1_bary - u2_bary;

                C_I_triplets.push_back(Eigen::Triplet<double>(curr_C_I_row, u0, u0_bary));
                C_I_triplets.push_back(Eigen::Triplet<double>(curr_C_I_row, u1, u1_bary));
                C_I_triplets.push_back(Eigen::Triplet<double>(curr_C_I_row, u2, u2_bary));

                igl::Hit back_hit = hits[hit_no + 1];
				int back_face = back_hit.id;

                // Back face vertices.
                int v0 = F(back_face, 0);
                int v1 = F(back_face, 1);
                int v2 = F(back_face, 2);

                // Back face barycentric coordinates.
                double v1_bary = back_hit.u;
                double v2_bary = back_hit.v;
                double v0_bary = 1 - v1_bary - v2_bary;

                C_I_triplets.push_back(Eigen::Triplet<double>(curr_C_I_row, v0, -v0_bary));
                C_I_triplets.push_back(Eigen::Triplet<double>(curr_C_I_row, v1, -v1_bary));
                C_I_triplets.push_back(Eigen::Triplet<double>(curr_C_I_row, v2, -v2_bary));

                curr_C_I_row++;
            }
        }
   	}
}

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
        const Eigen::VectorXi &mu,
        Eigen::MatrixXd &DV) {

    // In the paper, an alpha value of 10^-7 was used.
    double ALPHA = 0.0000001;

	int num_vertices = V.rows();
	int mu_len = mu.maxCoeff();
	DV.resize(num_vertices, 3);

	// Get the mass matrix with Voronoi areas, and square root each entry.
	Eigen::SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    M = M.cwiseSqrt();

	// D_A(i, j) = sqrt(M_i^0) if i = j else 0, where M_i is the Voronoi area of vertex i.
	Eigen::SparseMatrix<double> D_A(3 * num_vertices, 3 * num_vertices);
    igl::repdiag(M,3,D_A);

	// D_w(i, j) = w_i if i = j else 0.
	Eigen::SparseMatrix<double> D_w(3 * num_vertices, 3 * num_vertices);
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> w_matrix;
    w_matrix.diagonal() = weights;
    igl::repdiag((Eigen::SparseMatrix<double>)w_matrix,3,D_w);

	// D_V_tilde = Diagonal matrix containing row-wise stacked elements of V.
	Eigen::DiagonalMatrix<double,Eigen::Dynamic> D_V_tilde;
    Eigen::MatrixXd V_tilde(num_vertices, 3);

    V_tilde = V.rowwise() - view.transpose();
    D_V_tilde.diagonal() = Eigen::Map<const Eigen::VectorXd,Eigen::RowMajor>(V_tilde.data(),V_tilde.size());

    //Lambda_0 is the per vertex lambda value (||v||) of the undeformed mesh
    Eigen::VectorXd Lambda_0(num_vertices);
    Lambda_0 = V_tilde.rowwise().norm();
    for (int i = 0; i < num_vertices; i++)
    {
        V_tilde.row(i) /= Lambda_0[i];
    }

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
			for (int k = 0; k < 3; k++) {
				L_triplets.push_back(Eigen::Triplet<double>(3 * it.row() + k, 3 * it.col() + k, it.value()));
			}
		}
	}
	L_0_tilde.setFromTriplets(L_triplets.begin(), L_triplets.end());

	// S = selector matrix = I_n (x) [1, 1, 1]^T.
	Eigen::SparseMatrix<double> S(3 * num_vertices, num_vertices);
	std::vector<Eigen::Triplet<int>> S_triplets;
	S_triplets.reserve(3 * num_vertices);
	for (int vertex = 0; vertex < num_vertices; vertex++) {
		for (int i = 0; i < 3; i++) {
			S_triplets.push_back(Eigen::Triplet<int>(3 * vertex + i, vertex, 1));
		}
	}
	S.setFromTriplets(S_triplets.begin(), S_triplets.end());

    //intermediate step to make compiler happy
    Eigen::VectorXd S_Lambda_0,L_theta;
    S_Lambda_0 = S * Lambda_0;

    L_theta = S_Lambda_0.asDiagonal().inverse() * L_0_tilde * D_V_tilde * S_Lambda_0;

    //D_L_theta : 3n x 3n sparsematrix whose diagonal is L_theta.
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> D_L_theta;
    D_L_theta.diagonal() = L_theta;

    //size 3N X N
    Eigen::SparseMatrix<double> F_top_left = D_A * D_w * (L_0_tilde * D_V_tilde - (Eigen::SparseMatrix<double>) D_L_theta) * S;

    //the combined matrix F (name is used so I'm using A here)
	Eigen::SparseMatrix<double> A(3 * num_vertices + mu_len,num_vertices + mu_len);
    std::vector<Eigen::Triplet<double>>A_triplets;
	A_triplets.reserve(F.nonZeros() + mu_len);
    //setting the top left corner
    for (int k=0; k<F_top_left.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(F_top_left,k); it; ++it){
            A_triplets.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
        }
    }

    //setting the bottom right corner
    double sqrt_alpha = sqrt(ALPHA);
    for (int i = 0; i < mu_len; i++)
    {
        A_triplets.push_back(Eigen::Triplet<double>(3 * num_vertices + i, num_vertices + i, sqrt_alpha));
    }

    A.setFromTriplets(A_triplets.begin(),A_triplets.end());

    // size n + mu_len X n + mu_len
    Eigen::SparseMatrix<double> A_T_A = A.transpose() * A;

    Eigen::VectorXd f = Eigen::VectorXd::Zero(num_vertices + mu_len);

	Eigen::VectorXd lambda(num_vertices);
    Eigen::VectorXd x (num_vertices + mu_len);
    x.setZero();

    // Fix a subset of the vertices so that their lambda is the average of lambda_min and lambda_max.
    Eigen::VectorXd fixed_verticies_values(bf.size());
    for (int i = 0; i < bf.size(); i++)
    {
    	// Search lambdaMin and lambdaMax to find the lambda bounds for this vertex.
    	double lambdaMin_i = 0;
    	double lambdaMax_i = 0;
    	bool found_lambda_min = false;
    	bool found_lambda_max = false;
    	for (int lambda_row = 0; lambda_row < lambdaMin.rows(); lambda_row++) {
    		if (lambdaMin(lambda_row, 0) == bf[i]) {
    			lambdaMin_i = lambdaMin(lambda_row, 1);
    			found_lambda_min = true;
    		}
    		if (lambdaMax(lambda_row, 0) == bf[i]) {
    			lambdaMax_i = lambdaMax(lambda_row, 1);
    			found_lambda_max = true;
    		}
    		// If we have found both lambda bound values for this vertex, exit.
    		if (found_lambda_min && found_lambda_max) {
    			break;
    		}
    	}

    	// Set the fixed lambda value to the average of the min and max.
    	double avg_lambda_min_max = (lambdaMin_i + lambdaMax_i) / 2.0;
        fixed_verticies_values[i]  = avg_lambda_min_max;
    }

    std::vector<Eigen::Triplet<double>> C_I_triplets;

    // Shoot a raycast from the viewpoint in the direction of each vertex.
    int curr_C_I_row = 0;

    /* Note:  This function is meant to add in the depth constraints as described in the "Self-intersection avoidance" section
     * of the paper.  However it was giving us issues, so we decided not to include it.
     */
   	depth_constraint(num_vertices, view, V, F, curr_C_I_row, C_I_triplets);

    for (int i = 0; i < num_vertices; i++) {
        int vertex = lambdaMin(i, 0);
        double lambda_min_i = lambdaMin(i, 1);
        double lambda_max_i = lambdaMax(i, 1);

        // Constraint for lambda_i <= mu_g * lambda_i_max.
        C_I_triplets.push_back(Eigen::Triplet<double>(curr_C_I_row, vertex, 1));
        C_I_triplets.push_back(Eigen::Triplet<double>(curr_C_I_row, num_vertices + mu[i] - 1, -lambda_max_i));
        curr_C_I_row++;

        // // Constraint for mu_g * lambda_i_min <= lambda_i.
        C_I_triplets.push_back(Eigen::Triplet<double>(curr_C_I_row, vertex, -1));
        C_I_triplets.push_back(Eigen::Triplet<double>(curr_C_I_row, num_vertices + mu[i] - 1, lambda_min_i));
        curr_C_I_row++;
    }

   	Eigen::SparseMatrix<double> C_I(curr_C_I_row, num_vertices + mu_len);

   	C_I.setFromTriplets(C_I_triplets.begin(), C_I_triplets.end());

   	Eigen::VectorXd d = Eigen::VectorXd::Zero(curr_C_I_row);

    Eigen::SparseMatrix<double> Aeq;
    Eigen::VectorXd Beq;

    // Use active_set as a quadratic programming solver.
    igl::active_set(
    	A_T_A,
    	f,
	    bf,
	    fixed_verticies_values,
	    Aeq,
	    Beq,
	    C_I,
	    d,
	    Eigen::VectorXd(),
	    Eigen::VectorXd(),
	    igl::active_set_params(),
	    x
	);

	// Set the vertices of the bas-relief shape.
	for (int i = 0; i < num_vertices; i++) {
		// Scale back the vertices by 1 / mu_g.
		double mu_g_inv = (1.0 / x(num_vertices + mu[i] - 1));
		DV.row(i) = mu_g_inv * x(i) * V_tilde.row(i);
	}
   	DV.rowwise() += view.transpose();
}
