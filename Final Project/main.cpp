#include "appearance_mimicking_surfaces.h"
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>


int main(int argc, char *argv[])
{
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	// Keeps track of whether we are showing original vertices or deformed vertices.
	bool showing_deformed = false;

	// Load input meshes
	igl::read_triangle_mesh(
		(argc > 1 ? argv[1] : "../data/bunny.off"), V, F);

	int num_vertices = V.rows();
	// V.rowwise() += Eigen::Vector3d::Ones().transpose();
	// Deformed model vertices.
	Eigen::MatrixXd DV;

	// Set the lambda_min, lambda_max constant inequality constraints so
	// that each vertex is constrained between lambda_min and lambda_max.
	int num_constraints = num_vertices;
	Eigen::MatrixXd lambda_min(num_constraints, 2);
	Eigen::MatrixXd lambda_max(num_constraints, 2);
	for (int i = 0; i < num_vertices; i++) {
		lambda_min(i, 0) = i;
		lambda_min(i, 1) = .3;
		lambda_max(i, 0) = i;
		lambda_max(i, 1) = 1;
	}

	// Constrain the first vertex.
	int num_fixed_vertices = 1;
	Eigen::VectorXi bf(num_fixed_vertices);
	bf(0) = 0;

	Eigen::VectorXd weights(num_vertices);
	Eigen::VectorXi mu(num_vertices);
	for (int i = 0; i < num_vertices; i++) {
		weights(i) = 1;
		mu(i) = 1;
	}

	igl::opengl::glfw::Viewer viewer;

	// From https://github.com/libigl/libigl/issues/1201.
	auto &core = viewer.core();
	Eigen::Vector3f camera_pos;
	camera_pos << core.camera_eye;
	std::cout << "Camera position:  " << camera_pos << std::endl;

	// Set the viewpoint based on the camera's location.
	Eigen::Vector3d view;
	// view = (Eigen::Vector3d) camera_pos;
	view << (double) camera_pos[0], (double) camera_pos[1], (double) camera_pos[2];
	view << 1,0,1;
	appearance_mimicking_surfaces(V, F, view, lambda_min, lambda_max, bf, weights, mu, DV);

	std::cout << R"(
b,B 	Toggle bas-relief of the input mesh.
)";
	viewer.callback_key_pressed =
		[&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
	{
		switch(key)
		{
			case 'b':
			case 'B':
				if (showing_deformed) {
					viewer.data().set_vertices(V);
				} else {
					viewer.data().set_vertices(DV);
				}
				showing_deformed = !showing_deformed;
				break;
			default:
				return false;
		}
		return true;
	};

	viewer.data().set_mesh(V, F);

	// viewer.data().show_lines = false;
	// viewer.data().show_overlay = false;
	// viewer.data().face_based = false;
	viewer.launch();
	return EXIT_SUCCESS;
}
