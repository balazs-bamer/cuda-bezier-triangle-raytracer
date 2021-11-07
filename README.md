# Cuda-based Bezier triangle mesh raytracer

The aim of this project is to simulate image rendering of arbitrary-shaped "lenses" described by triangle meshes using Bezier triangle interpolation. It will be possible to simulate aspheric and anomorphic lenses, and of course spherical lenses as well. (Of course spherical lenses can be simulated much more efficiently.)

This work is divided into 3 parts:
1. C++ utilities and preprocessor. These parts check the input and produce the Bezier mesh appropriate to implement raytracing.
2. Reference single ray tracing implementation in C++. I use it actually to develop all geometry and raytracing algorithms, which will later be modified for GPU usage.
3. Initial implementation in NVIDIA Thrust, a quick proof-of-concept implementation based on the reference to produce images.
4. An iteratively developed CUDA version for best performance.

I use Eigen3 for vector calculations, because it is well supported under CUDA, and Google Test for unit testing. I also use Sebastian Reiter's STL reader utility to parse STL files.

[Here](https://db.bme.hu/~bamer/tdk/TDK.html) is my ancient naive implementation of spherical lens system simulator from 1994-1998 (in Hungarian). That work originates from the high school with most parts ready then. I've abandoned the work because a 80386 without an FPU was far too slow for any meaningful rendering. During the university I could complete it and use a Sun workstation for rendering the images.

I use At most C++17 features all thorough the code to let it interoperate with CUDA.

## C++ utilities and preprocessor

All classes and types and templated to be able to try them with float and double types. The aim is float, because it is much more efficient on GPU and more data fits in the cache. Double may be interesting for debugging if numerical problems would occur. For simplicity I omit the template parameters here.

For space and performance reasons, I use 3D Cartesian coordinates. Since most container sizes are known or calculable, I use `std::vector`s for simple storage. This will let me directly pass the data to GPU code.

### reference/util.h

This file contain common type declarations for all other sources, some geometry calculations like plane creatian and intersection, triangle altitudes. There is also a simple uniform triangle subdivision function template.

### reference/mesh.h

This is the home of the `Mesh` class template. I know there are several good STL and mesh processing libraries, but I've created this class for three reasons:
1. I wanted full control over mesh processing.
2. It contains important preprocessing steps for later stages, which belong to triangle mesh processing, not Bezier triangle meshes.
3. I enjoy working with geometry vectors.

The `Mesh` class exposes a set of STL vector-like interface functions for easier use. It also provides overloaded operators for displace and inflate operations, STL file I/O. Here I only write about the more interesting functions.
Note, this class can only handle a mesh consisting of **exactly one totally filled** shape.

#### Mesh::standardizeVertices

This function is used to make vertices of different triangles coincide when they are practically the same mesh vertex, but have sightly different values from different numeric errors. This will enable later functions to compare floating point values for identity check. Note, these comparisons only occur on output of this function without any subsequent calculations. Here first I calculate an epsilon value for proximity comparisons using the smallest triangle length. Since the naive implementation would be O(n^2), a more sophisticated algorithm is needed.

I solve the problem using 1D projections. For each dimension, I take the vertices in ascending order and group them based on if they are closer to each other than epsilon. Pairwise proximity check will then only be needed in each group. I take the dimension where the maximum groups size is minimal (less degenerate projection) and then perform the pairwise check. I adjust each vertex pair to the smaller one according to an artificial alphabetical ordering of the coordinates.

#### Mesh::standardizeNormals

This function makes all triangle normals point outwards if calculated as `(triangle[1] - triangle[0]).cross(triangle[2] - triangle[0])` (note: `using Triangle = std::array<Vertex, 3u>`). This funtion must be used after `Mesh::standardizeVertices`, and assumes that there is no triangle vertex on the side of an other triangle. The algorithm is roughly
1. Find a vertex with the smallest X coordinate.
2. Find a neighbouring triangle with the smallest x component (negative for all such triangles).
3. Make this triangle good by possibly swapping two vertices.
4. Create a set of triangles with the side neighbours of this initial triangle, marking their known good neighbour as well.
5. Take a triangle from the set, make it good using the info from its known neighbour.
6. Put its unprocessed neighbours to the set.
7. Loop to 5. if the sset is not empty.

This function also gathers auxiliary data for later stages, see below.

#### Mesh::makeEllipsoid(int32_t const aSectors, int32_t const aBelts, Vector const &aSize)

Deletes all current contents and creates an ellipsoid with the given parameters. The triangle vertices will be on the ellipsoid surface. The ellipsoid is approximated somewhat similarly as the longitude-latitude mesh on the Earth. There will be `aBelts` pieces of "belts" along each "parallel", each divided evenly by `aSectors` vertices. Between the belts there are 2 * `aSectors` pieces of triangles, while between the extreme belts and poles only `aSectors` pieces. The parameter aSize contains the principal semi-axes of the ellipsoid. There is also a special case called `makeUnitSphere` to create a unit sphere.

#### Mesh::getFace2neighbours

Returns an `std::vector<Neighbours>`, each item corresponding to a mesh triangle. Here the struct `Neighbours` holds information about the neighbouring triangles of the aactual one's each side (starting at vertex i within the triangle):
* mFellowTriangles for mesh index of the neighbour of the actual side (`i` -> `(i+1)%3`)
* mFellowCommonSideStarts for the vertex index `k` in each fellow triangle where the common side starts such that the side in the neighbouring triangle (`k` -> `(k+1)%3`) is the same side as (`i` -> `(i+1)%3`).
Must be run after `Mesh::standardizeNormals`.

#### Mesh::getVertex2averageNormals

Returns a map of vertices to the "vertice normals", which will be the Bezier triangle normals in the vertices. "Vertex normals" are the weighted average of triangle normals containing the vertex, weighted by the triangle angles at thaty vertex. Must be run after `Mesh::standardizeNormals`.

### reference/bezierTriangle.h

Here resides the `BezierTriangle` class which is capable of
* Creating a Bezier triangle over a regular triangle, considering its neighbours, in a way that the neighbouring Bezier surfaces will give a surface with C1 continuity. [[1]](#1) I chosed C1 because cubic Bezier triangles are sufficient for this, and the interpolation and thus the surface intersection calculations will be faster.
* Linear interpolation on the underlying planar triangle using barycentric coordinates whithin the underlying triangle.
* Bezier interpolation using barycentric coordinates whithin the underlying triangle.
* Ray and underlying planar triangle intersection with barycentric output. TODO implement.
* Ray and Bezier triangle intersection with barycentric output. TODO implement.
* Normal vector calculation on the barycentric triangle using barycentric coordinates. TODO implement.

The Bezier triangle will contain all the original triangle vertices. Control point calculations follow [[1]](#1) with the details figured out by myself where the paper was not specific. There are 3 `constexpr` parameters influencing control point placement, which were empirically estimated TODO see where
TODO perhaps write control point calculation in detail.

Apart of the constructor, important public interface is the `interpolateLinear` and `interpolate` (Bezier) function families, which are self-explanatory.

### reference/bezierMesh.h

The `BezierMesh` class is similar to and based on the `Mesh` class, and is responsible of
* Constructing the Bezier triangle mesh of `BezierTriangle` instances in Clough-Tocher subdivision [[1]](#1), see later. TODO docs.
* Obtaining a triangular mesh approximation (`interpolate`) with evenly splitting each subtriangle side into `aDivisor` parts.
* Split "thick" Bezier triangles into smaller, "thinner" ones.

#### BezierMesh::BezierMesh

TODO docs.

#### BezierMesh::splitThickBezierTriangles

TODO docs.

#### Ellipse approximation test

I've chosen an ellipsoid with principal semi-axes 1.0, 2.0 and 4.0 was used for these tests. I've chosen this shape because it somewhat resembles a lens, has curvatures with reasonably high variety if radii. I've measured the error of the interpolated mesh (`aDivisor` == 3) vertices relative to the ellipsoid surface point in the very same direction (where the line containing the center and the vertex intersects the ellipsoid). Average quadratic relative errors were between **1.9e-5** and **2.2e-3** for the same tuned parameters.
TODO add details.

#### Shortcomings

In theory, my algorithms can handle concave meshes, and the `Mesh` class even meshes with holes in it (with hole walls covered by triangles). In practice, I've tried so far only ellipsoids, and a fairly complicated machine part STL, which caused numeric errors, most probably because the extreme variation between sizes of neighbouring triangles and sharp edges. Such shapes are not quite practical for photographic image rendering.

## References

See [Wikipedia](https://wikipedia.org) for well-known or simpler stuff.

<a id="1">[1]</a> 
Shaoming Wang (2004). 
A smooth surface interpolation to 3D triangulations
Journal of Computational and Applied Mathematics, Volume 163, Issue 11, February 2004, pp 287–293

