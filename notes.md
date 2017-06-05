The purpose of having lighting and shading is to allow us to generate things that look like something in the real world
There are several different approaches to lighting. Raytracing is one of them, but we won't be doing it
We'll be looking at a different model.
Color is determined by
* The reflective properties of objects. Clothing, for example, will have a different reflection from a metallic surface (matte vs glossy)
* The color of the light, generally we don't think about this because we usually use white light, but different color lights can have interesting effects too

2 Kinds of light sources
* Ambient light
	* No particular source, permeates image with equal intensity
* Point light
	* Lights that come from a particular pointNo


**Calculating Colors**
I: color value for a surface (illumination)

Lighting Equation
I = I ambient + I diffuse + I specular

Diffuse = dull, matte reflection

**Ambient Reflection**
* A: Ambient Light (0-255)
* kA: constant of ambient reflection (0-1)
* Iambient= A*kA done 3 times

**Diffuse Reflection (matte)**
* L: Point light source <x, y, z>, (0-255)
* kD: constant of diffuse reflection (0-1)
* Reflected light is evenly distributed in all directions, but the light is only coming from one source (in contrast to ambient reflection)
* When is the reflection the strongest? When the theta between the point light source and the normal of the surface is 0
* Idiffuse = L * kD * cos(theta)
* How do we find cos(theta)? cos(theta) = Normal Vector * Light Vector if vN and vL are unit vectors
* We will thus need to normalize the vectors, creating a unit vector
* NV = Vx/||V||, Vy/||V||, Vz/||V||
* ||V|| = (Vx^2+Vy^2+Vz^2)^(1/2)

**Specular Reflection**
* L: Point light source
* kS: constant of specular reflection
* In contrast to diffuse reflection specular reflection reflects the light strongly in one direction
* Since reflections are not scattered randomly for specular reflection, we must consider the viewing angle
* The angle between the viewing vector and the reflection vector is what determines the reflection intensity this time around
* Ispecular = L * kS * cos(phi)
* cos(phi) = vR * vV
* R = 2P-L
* R = 2N(N*L) - L
* To find R we have to first find P, which is the projection of L onto N (the direction of N but the y endpoint of L)
* We also need vector S, which starts from the endpoint of vector P to the endpoint of Vector R
* We thus need to find vector S
* R = P+S
* S = P-L
* P = N(N*L)
* R = P+P-L
* R=2P-L
* R=2N(N*L)-L
* cos(phi)=(2N(N*L)-L)*V
* Ispecular = L * kS * [(N * 2(N*L) - L) * V] ^ p
* p determines how quickly the reflection fades
* Make sure that the dot products should be >=0, if it is negative, change it to 0


Add Ispecular, Idiffuse, and Iambient to calculate our final color

**Shading Model**
How often do we do these calculations?

* Flat shading
	* Calculate I once per polygon
	* Will need high polygon count to look good
	* Certain objects just won't work well, like boxes
* Goroud
	* Generate a new color for every pixel, but we do not calculate an I for every pixel
 	* Create a list of *vertex normal* values
 	* Instead of relying on surface normals alone, we will use vertex normals
 	* Vertex normals are the sum of the surface normals of the polygons that share the same vertex
 	* Vertex Normal: Normalized sum of all surface normals that share a vertex
 	* Calculate I for each vertex normal of a polygon
 	* Enables smoother colors in shading without increasing polygon count
 	* Generate new I values in scanline and draw_line (similar to finding x and z)

* Phong
	* Calculating I for each pixel
	* High computation intensity
	* Lookup vertex normals for each polygon vertex
	* Generate new normal values in scanline and drawline
	* Generate new normal values in scannline and draw_line

**draw.c**
    * flat_shade
    * phong_shade
    
** write mdl commands into my_main.c**
    * LIGHT
    * CONSTANTS
    * AMBIENT
    * //SHADING 