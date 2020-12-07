#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>


int main(int argc, char *argv[])
{
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	// Load input meshes
	igl::read_triangle_mesh(
		(argc > 1 ? argv[1] : "../data/bunny.off"), V, F);

	// TODO:  Setup stuff here.

	igl::opengl::glfw::Viewer viewer;

	// TODO:  Replace these options with ones for us.
	std::cout<<R"(
<letter>	<action>
...
)";
	const auto update = [&]()
	{
		// TODO:  What should happen to the mesh each loop?
		// viewer.data().set_vertices(V);
	};
	viewer.callback_key_pressed =
		[&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
	{
		switch(key)
		{
			case 'l':
				break;
			default:
				return false;
		}
		update();
		return true;
	};

	viewer.data().set_mesh(V, F);

	// const Eigen::RowVector3d orange(1.0,0.7,0.2);
	// const Eigen::RowVector3d yellow(1.0,0.9,0.2);
	// const Eigen::RowVector3d blue(0.2,0.3,0.8);

	// viewer.data().set_edges(lP,lE,lC);

	update();
	// viewer.data().show_lines = false;
	// viewer.data().show_overlay = false;
	// viewer.data().face_based = false;
	viewer.launch();
	return EXIT_SUCCESS;
}
