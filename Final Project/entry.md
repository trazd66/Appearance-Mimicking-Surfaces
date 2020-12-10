# Appearance Mimicking Surfaces
[link to paper](https://cims.nyu.edu/gcl/papers/mimicking-2014.pdf)

## Summary:

A bas-relief is a model that depicts a recognizable shape, but with a limited volume. This work aims to tackle the problem of producing an "appearance-mimicking surface", a specific type of bas-relief which appears similar to the original shape when viewed from a specific angle.

  

In order to accomplish this, we would like to minimize the sum of squared distance between the normal of each point on the original surface and its corresponding point on the output bas-relief. 
(See details from the paper.)

## Implementation :
Libraries used : Eigen , libigl
We followed the matrix construction steps from the paper to construct the quadratic programming problem.
For the quadratic programming solver, we used `igl::active_set`
(The authors of the paper mentioned that it can also be turned into a conic programming problem and solved perfectly with `Mosek`'s Conic solver.)

## Header
We implemented a single file implementation of the algorithm that can be included and used easily.
The header looks like this. 

    // Given a mesh, a viewpoint, and linear constraints:
    
    // lambdaMin -- combined with lambdaMax this will "squeeze" the mesh into the defined thickness by projecting the vertices
    
    // lambdaMax
    
    // mu -- different parts of the shape can be given different thickness by segmenting the vertices into groups
    
    // fixed vertices
    
    // weights -- scales the difference of vertex normals so that more visible vertices are given preserved "more" than less visible
    
    // optional depth ordering
    
    // output a bas-relief deformed mesh whose appearance is preserved
    
    // from the viewpoint
    
    //
    
    //
    
    // Inputs:
    
    // V #V by 3 list of the vertex positions of the model
    
    // F #F by 3 list of triangle indices into V
    
    // view 3D vector of the coordinates of the viewpoint
    
    // lambdaMin #Constraints by 2 list of vertex indices and the lambdaMinValue: [vertexIdx,lambdaMinValue]
    
    // lambdaMax #Constraints by 2 list of [vertexIdx,lambdaMaxValue]
    
    // bf #fixedVertices vector of indices into V
    
    // weights #V length list of vertex weights
    
    // mu #V length list of mu indices
    
    // (separates mesh into independent regions each with their own thickness constraint)
    
    // (default should be all 1's)
    
    // Outputs:
    
    // DV #V by 3 list of the vertex positions of the deformed model
    



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

## Sample results
1. with `lambdaMin` and `lambdaMax` set to the same, so it should squash the bunny into a thin sheat of paper while still preserving its shape.



