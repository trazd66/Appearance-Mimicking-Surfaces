# Appearance Mimicking Surfaces
[Link to paper](https://cims.nyu.edu/gcl/papers/mimicking-2014.pdf)

## Summary:
A bas-relief is a model that depicts a recognizable shape, but with a limited volume. This work aims to tackle the problem of producing an "appearance-mimicking surface", a specific type of bas-relief which appears similar to the original shape when viewed from a specific angle.

In order to accomplish this, we would like to minimize the sum of squared distances between the normal of each point on the original surface and its corresponding point on the output bas-relief.
(See the paper for technical details.)

## Implementation:
Libraries used: Eigen, Libigl.
- We followed the matrix construction steps from the paper to construct the quadratic programming problem.
- For the quadratic programming solver, we used `igl::active_set`.
(The authors of the paper mentioned that it can also be turned into a conic programming problem and solved more efficiently with `Mosek`'s Conic solver.)

- We did not implement the Z-depth constraints due to limited time , hence this implementation is flawed and the results may contain unwanted artifacts

## Header
We implemented a single file implementation of the algorithm that can be included and used easily.

## Example usage

    // Load input meshes

    igl::read_triangle_mesh((argc > 1 ? argv[1] : "../data/bunny.off"), V, F);

    int  num_vertices = V.rows();

    // V.rowwise() += Eigen::Vector3d::Ones().transpose();

    // Deformed model vertices.

    Eigen::MatrixXd  DV;



    // Set the lambda_min, lambda_max constant inequality constraints so

    // that each vertex is constrained between lambda_min and lambda_max.

    int  num_constraints = num_vertices;

    Eigen::MatrixXd  lambda_min(num_constraints, 2);

    Eigen::MatrixXd  lambda_max(num_constraints, 2);

    for (int  i = 0; i < num_vertices; i++) {

    lambda_min(i, 0) = i;

    lambda_min(i, 1) = 1;

    lambda_max(i, 0) = i;

    lambda_max(i, 1) = 1;

    }



    // Constrain the first vertex.

    int  num_fixed_vertices = 1;

    Eigen::VectorXi  bf(num_fixed_vertices);

    bf(0) = 0;



    Eigen::VectorXd  weights(num_vertices);

    Eigen::VectorXi  mu(num_vertices);

    for (int  i = 0; i < num_vertices; i++) {

    weights(i) = 1;

    mu(i) = 1;

    }



    igl::opengl::glfw::Viewer  viewer;


    Eigen::Vector3d  view;

    view  <<  1,0,1;

    appearance_mimicking_surfaces(V, F, view, lambda_min, lambda_max, bf, weights, mu, DV);


## Compile and using the main file.
1. download the most recent `libigl` library into the directory
2. run `cmake .. -DCMAKE_BUILD_TYPE=Release` or `cmake .. -DCMAKE_BUILD_TYPE=Debug` for debug mode inside the `build` directory
3. run `make` inside the `build` directory
4. run `./appearance-mimicking-surfaces ` for the default bunny mesh, or add an argument to the your prefered mesh .e.g:  `./appearance-mimicking-surfaces ../data/cactus.obj`



## Sample results
1. with `lambdaMin` and `lambdaMax` set to the same, so it should squash the bunny into a thin sheat of paper while still preserving its shape.

<p align="center">
<img src="./sample pictures/bunny0.png"/>
<p align="center">
<img src="./sample pictures/bunny1.png"/>

2. with less extreme `lambdaMin` and `lambdaMax`

<p align="center">
<img src="./sample pictures/bunny2.png"/>
<p align="center">
<img src="./sample pictures/bunny3.png"/>
