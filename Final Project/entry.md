#Appearance Mimicking Surfaces

###Summary
A bas-relief is a model that depicts a recognizable shape, but with a limited volume.  This work aims to tackle the problem of producing an "appearance-mimicking surface", a specific type of bas-relief which appears similar to the original shape when viewed from a specific angle.

In order to accomplish this, we would like to minimize the sum of squared distance between the normals of each point on the original surface and its corresponding point on the output bas-relief.  Mathematically, we want to minimize:

<Include Equation (2) here>

where S^0 is the original shape, S is the bas-relief, and ...