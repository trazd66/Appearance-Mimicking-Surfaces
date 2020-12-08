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

	// Deformed model vertices.
	Eigen::MatrixXd DV;

	// Set the viewpoint to the origin.
	Eigen::Vector3d view;
	view << 0, 0, 0;

	Eigen::MatrixXd lambdaMin, lambdaMax;

	int num_fixed_vertices = 1;
	Eigen::VectorXi bf(num_fixed_vertices);
	Eigen::VectorXd weights(num_vertices);
	Eigen::VectorXd mu(num_vertices);
	for (int i = 0; i < num_vertices; i++) {
		weights(i) = 1;
		mu(i) = 1;
	}

	appearance_mimicking_surfaces(V, F, view, lambdaMin, lambdaMax, bf, weights, mu, DV);

	igl::opengl::glfw::Viewer viewer;

	std::cout << R"(
b,B 	Toggle bas-relief of the input mesh.
)";
	const auto update = [&]()
	{

	};
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
				break;
			default:
				return false;
		}
		update();
		return true;
	};

	viewer.data().set_mesh(V, F);

	update();
	// viewer.data().show_lines = false;
	// viewer.data().show_overlay = false;
	// viewer.data().face_based = false;
	viewer.launch();
	return EXIT_SUCCESS;
}
