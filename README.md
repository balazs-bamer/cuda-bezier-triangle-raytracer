# Cuda-based Bézier triangle mesh raytracer

The aim of this project is to simulate image rendering of arbitrary-shaped "lenses" described by triangle meshes using Bézier triangle interpolation. It will be possible to simulate aspheric and anomorphic lenses, and of course spherical lenses as well. (Of course spherical lenses can be simulated much more efficiently.) This project is just for personal learning, so I omit investigation of existing results.

This work is divided into 3 parts:
1. C++ utilities and preprocessor. These parts check the input and produce the Bézier mesh appropriate to implement raytracing.
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
2. I wanted to preserve original vertices. If the mesh creation algorithm is sophisticated enough, it knows better where and how to place vertices. Any postprocessing can only make it worse.
3. It contains important preprocessing steps for later stages, which belong to triangle mesh processing, not Bézier triangle meshes.
4. I enjoy working with geometry vectors.

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

Deletes all current contents and creates an ellipsoid with the given parameters. The triangle vertices will be on the ellipsoid surface. The ellipsoid is approximated somewhat similarly as the azimuth-inclination mesh on the Earth. There will be `aBelts` pieces of "belts" along each "parallel", each divided evenly by `aSectors` vertices. Between the belts there are 2 * `aSectors` pieces of triangles, while between the extreme belts and poles only `aSectors` pieces. The parameter aSize contains the principal semi-axes of the ellipsoid. There is also a special case called `makeUnitSphere` to create a unit sphere.

#### Mesh::getFace2neighbours

Returns an `std::vector<Neighbours>`, each item corresponding to a mesh triangle. Here the struct `Neighbours` holds information about the neighbouring triangles of the aactual one's each side (starting at vertex i within the triangle):
* mFellowTriangles for mesh index of the neighbour of the actual side (`i` -> `(i+1)%3`)
* mFellowCommonSideStarts for the vertex index `k` in each fellow triangle where the common side starts such that the side in the neighbouring triangle (`k` -> `(k+1)%3`) is the same side as (`i` -> `(i+1)%3`).
Must be run after `Mesh::standardizeNormals`.

#### Mesh::getVertex2averageNormals

Returns a map of vertices to the "vertice normals", which will be the Bézier triangle normals in the vertices. "Vertex normals" are the weighted average of triangle normals containing the vertex, weighted by the triangle angles at thaty vertex. Must be run after `Mesh::standardizeNormals`. Perhaps it would be important to implement a way to import these values along with the triangles when the original surface is known and dervatives are present. This would enable much more exact approximation.

### reference/bezierTriangle.h

Here resides the `BézierTriangle` class which is capable of
* Creating a Bézier triangle over a regular triangle, considering its neighbours, in a way that the neighbouring Bézier surfaces will give a surface with C1 continuity. [[1]](#1) I chosed C1 because cubic Bézier triangles are sufficient for this, and the interpolation and thus the surface intersection calculations will be faster.
* Linear interpolation on the underlying planar triangle using barycentric coordinates within the underlying triangle.
* Bézier interpolation using barycentric coordinates within the underlying triangle.
* Ray and underlying planar triangle intersection with barycentric output.
* Ray and Bézier triangle intersection with barycentric output.
* Normal vector calculation on the barycentric triangle using barycentric coordinates.

The Bézier triangle will contain all the original triangle vertices. Control point calculations follow [[1]](#1) with the details figured out by myself where the paper was not specific. There are 3 `constexpr` parameters influencing control point placement, which were empirically estimated TODO see where
TODO perhaps write control point calculation in detail.

The class contains some precomputed values useful for intersection calculations, like maximum distances from the underlying triangle inside and outside.

#### BezierTriangle::intersection

Theres is no closed formula for Bézier triangle and ray intersection, like there is for planes and spheres. There are existing ray - Bézier triangle intersection algorithms, such as [[2]](#2). However, this algorithm needs investigation of several cases, which is not well suitable for GPUs. Moreover, the article does not contain performance data, and it seemed to me a big effort to implement it. So I've found out a rather simple algorithm with one or two identical computation-intensive loops, which are easy to run parallel for many rays.

First I check the intersection with the underlying planar triangle, and if it intersects, calculate the lengths along the ray measured from its starting points
* to the planar intersection
* to the maximum possible surface point outside
* to the minimum possible surface point inside - these two come from the maximum distances inside and outside and the incidence to the underlying plane.

I determine whether the inside or the outside case is closer to the ray start point, and calculate first the intersection with the closer one. If there is none, I take the further one. For either case, I pass two of the above distance extremes which are parameters of a binary search.

For any case, the intersection calculation is performed in a fixed number of iterations. The number depends on the wanted accuracy. The interseciton calculation will only take place when the ray and Bézier surface distance from the underlying plane "change magnitude" along the ray in the given interval. If so, the following binary search yields the result:
1. Set bias = (0, 0, 0) 3D Cartesian vector.
2. Project the candidate (middle) point of the ray to the underlying plane, and translate it with the bias (see later why).
3. Calculate its barycentric coordinates.
4. Use these to calculate the corresponding surface point.
5. Project it back to the underlying plane. This is needed, because the projection is in general different from that we started from.
6. Set bias = last ray point projection - surface point projection.
7. Goto 2 when there is iteration left.

It is easy to see that the process converges. When ready, it is important to check if the intersection is within the domain of the current Bézier triangle. If not, the intersection would be more accurate when calculated on the appropriate neighbouring triangle. To enable this, in this case I return the side index to be considered for a similar calculation process for the neighbouring triangle. Since its planar intersection will definitely be outside the underlying triangle, I omit that check to let the function finish. TODO figures.

#### BezierTriangle::getNormal

This function is used in the previous one to calculate the surface normal in the intersection point. I used the principle described in [[3]](#3) to calculate closed-formula barycentric derivatives and use them to calculate two perpendicular surface tangent directions in that point. Their cross product gives the desired normal.

### reference/bezierMesh.h

The `BézierMesh` class is similar to and based on the `Mesh` class, and is responsible of
* Constructing the Bézier triangle mesh of `BézierTriangle` instances in Clough-Tocher subdivision [[1]](#1), see later. Input is plain triangle `Mesh`.
* Obtaining a triangular mesh approximation (`interpolate`) with evenly subdividing each subtriangle side into `aDivisor` parts.
* Split "thick" Bézier triangles into smaller, "thinner" ones.

#### BezierMesh::splitThickBezierTriangles

I've implemented this function because the raytracing will inspect mesh and ray intersections using the underlying triangle mesh (after Clough-Tocher subdivision). When we have the planar intersection, it will be used to calculate the ray intersection with the Bézier triangle above it. However, if the Bézier triangle forms a relatively too "tall" dome above the underlying triangle, it is more likely a ray can travel through it without intersecting any planar triangle. Of course it is still possible when the Bézier triangles are "close" to the underlying triangles, but much less likely. In such a rare case, the simulation will find the ray pass just right beside the object.

The idea here is to take all Bézier triangles, and subdivide each one proven to be too tall. Note, as we are in the Bézier triangle domain, each triangle is the reasult of the Clough-Tocher subdivision. I could have gone with taking that subdivision of the tall triangle, but it tends to create "thin" triangles, as the edges won't ever be split. I know there are well-known subdivision algorithm, but I wanted to make a very simple implementation and take advantage of the Bézier triangle control point information, which already takes into account the neighbouring triangles.

Curently I take the approximate Bézier triangle height (maximum at the original triangle centroid and quarter points of the original sides) and the triangle perimeter ratio, and if it is larger than a constant, I will split it. It might be interesting to incorporate an absolute value parameter too.

The algorithm is simple: I take the underlying side midpoints of a "tall" Bézier triangle, and divide it into 4 pieces. Now it creates vertices where the neighbouring triangle may not have on (no need to split it), so I propagate the subdivision info across each edge to the neighbouring triangle, and split them accordingly into 2 or 3 pieces, for 1 or 2 neighbours being subdivided, respectively. TODO figures.

New (dividing) vertices are calculated as a fixed linear combination of the side midpoint of the underlying triangle edge and the Bézier interpolation of that point (which is "above" the underlying triangle). The linear combination factor has been found out empirically, see below. The function gives back a new finer `Mesh` instance, which then needs to be preprocessed again and used in creation of a new, finer `BézierMesh`.

#### BezierMesh::intersect

This is now only a simple brute-force search for the intersection giving the shortest ray.

#### Ellipse approximation test

I've chosen an ellipsoid with principal semi-axes 1.0, 2.0 and 4.0 was used for these tests. I've chosen this shape because it somewhat resembles a lens, has curvatures with reasonably high variety in radii. I've measured the error of the interpolated mesh (`aDivisor` == 3) vertices relative to the ellipsoid surface point in the very same direction (where the line containing the center and the vertex intersects the ellipsoid). Average quadratic relative errors were between **1.9e-5** and **2.2e-3** for the same tuned parameters.
I find these results good enough to start with. Should the error be too big for some application, a more specialized parameter tuning is possible or I might use more sophiticated preprocessing algorithms.
There are general error limits for surface approximation using Bézier triangles like in [[4]](#4). Here for any triangle where the parametric function describing the surface to approximate is known and its second partial derivatives exist, an upper limit for the error can be expressed.

#### Shortcomings

In theory, my algorithms can handle concave meshes, and the `Mesh` class even meshes with holes in it (with hole walls covered by triangles, so holes in topological sense). In practice apart of ellipsoids, I have tried the algorithm only with a fairly complicated machine part STL. It caused numeric errors, most probably because the extreme variation between sizes of neighbouring triangles and sharp edges. Such shapes are not quite practical for photographic image rendering.

I have tested and adjusted algorithm parameters so far only with ellipsoids. TODO later construct some more complicated shapes, including concave ones.

The algorithm does not report mesh iontersection for large angles of incidence (above approximately 70 degrees).

## References

See [Wikipedia](https://wikipedia.org) for well-known or simpler stuff.

<a id="1">[1]</a> 
Shaoming Wang (2004). 
A smooth surface interpolation to 3D triangulations
Journal of Computational and Applied Mathematics, Volume 163, Issue 11 (February 2004), pp 287–293

<a id="2">[2]</a> 
S.H.M. Roth, P. Diezi, M.H. Gross (2000).
Triangular Bezier clipping
Proceedings the Eighth Pacific Conference on Computer Graphics and Applications (October 2000)

<a id="3">[3]</a>
Gerald Farin (1986).
Triangular Bernstein-Bézier patches
Computer Aided Geometric Design 3, pp 83-127

<a id="4">[4]</a>
Chang Geng-zhe, Feng Yu-yu (1983).
Error bound for Bernstein-Bézier triangulation approximation
Journal of Computational Mathematics Vol. 1, No. 4 (October 1983), pp. 335-340
