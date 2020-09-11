#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <exception>

#include "CommandLine.h"
#include "Image.h"
#include "Ray.h"
#include "Camera.h"
#include "LightSource.h"
#include "Mesh.h"
#include "Scene.h"
#include "RayTracer.h"
#include "phm.h"

using namespace std;



int main (int argc, char ** argv) {
	CommandLine args;
	if (argc > 1) {
		try {
			args.parse(argc, argv);
		} catch (const std::exception & e) {
			std::cerr << e.what() << std::endl;
			args.printUsage (argv[0]);
			exit(1);
		}
	}
	// Initialization

	Image image (args.width (), args.height ());
	Scene scene;
	

	Camera camera(Vec3f(0.2f, 0.5f, 10.f),
	//Camera camera(Vec3f(0.3f, 4.f, 0.f),
				  Vec3f(),
				  Vec3f(0.f, 1.f, 0.f),
				  60.f,
				  float(args.width()) / args.height());

	scene.camera() = camera;
	
	// Loading a mesh
	Mesh mesh;
	try {
		//mesh.loadOFF("../2_cubs.off");
		mesh.loadOFF("../bigcube_check.off");
	}
	catch (const std::exception & e) {
		std::cerr << e.what() << std::endl;
		exit(1);
	}
	cout << "mesh loaded " << endl;


	 int numVerts = mesh.vertexPositions().size();
    mesh.vertexPositions().insert( mesh.vertexPositions().end(),
                               {Vec3f(10.f,-1.5f,10.f),
                                Vec3f(10.f,-1.5f,-10.f),
                                Vec3f(-10.f,-1.5f,10.f),
                                Vec3f(-10.f,-1.5f,-10.f)});
    mesh.vertexNormals().insert( mesh.vertexNormals().end(),
                                {Vec3f(0.f,1.f,0.f),
                                 Vec3f(0.f,1.f,0.f),
                                 Vec3f(0.f,1.f,0.f),
                                 Vec3f(0.f,1.f,0.f)});
    mesh.indexedTriangles().push_back(Vec3i(numVerts,numVerts+1, numVerts+3));
 mesh.indexedTriangles().push_back(Vec3i(numVerts,numVerts+2, numVerts+3));



 	//adding walls
 	numVerts = mesh.vertexPositions().size();
 	mesh.vertexPositions().insert( mesh.vertexPositions().end(),
                               {Vec3f(6.f,-10.f,10.f),
                                Vec3f(6.f,-10.f,-10.f),
                                Vec3f(6.f,10.f,10.f),
                                Vec3f(6.f,10.f,-10.f)});
    mesh.vertexNormals().insert( mesh.vertexNormals().end(),
                                {Vec3f(1.f,0.f,0.f),
                                 Vec3f(1.f,0.f,0.f),
                                 Vec3f(1.f,0.f,0.f),
                                 Vec3f(1.f,0.f,0.f)});
    mesh.indexedTriangles().push_back(Vec3i(numVerts,numVerts+1, numVerts+3));
 	mesh.indexedTriangles().push_back(Vec3i(numVerts,numVerts+2, numVerts+3));

 	numVerts = mesh.vertexPositions().size();
 	mesh.vertexPositions().insert( mesh.vertexPositions().end(),
                               {Vec3f(-6.f,-10.f,10.f),
                                Vec3f(-6.f,-10.f,-10.f),
                                Vec3f(-6.f,10.f,10.f),
                                Vec3f(-6.f,10.f,-10.f)});
    mesh.vertexNormals().insert( mesh.vertexNormals().end(),
                                {Vec3f(-1.f,0.f,0.f),
                                 Vec3f(-1.f,0.f,0.f),
                                 Vec3f(-1.f,0.f,0.f),
                                 Vec3f(-1.f,0.f,0.f)});
    mesh.indexedTriangles().push_back(Vec3i(numVerts,numVerts+1, numVerts+3));
 	mesh.indexedTriangles().push_back(Vec3i(numVerts,numVerts+2, numVerts+3));
 
    scene.meshes ().push_back (mesh);

    //creating area light source. The code supports multiple light sorces passed as a vector
	Vec3f light_pos = Vec3f(-1.3, 5.0, 0.f);
 	LightSource lgtsr(light_pos, Vec3f(1.f,1.f,1.f), 0.f);
 	lgtsr.init_area();
 	lgtsr.initAreaSource();
	scene.lightsrc().push_back(lgtsr); 

	Vec3f albedo = Vec3f(167,40,40); 
	Material matr(albedo, 0.8);
	scene.material() = matr;



	RayTracer rayTracer;

	// Rendering
	image.fillBackground ();
	std::cout << "Ray tracing: starts";
	//rayTracer.render (scene, image);
	rayTracer.render_MC(scene, image);
	std::cout << "ends." << std::endl;
	image.savePPM (args.outputFilename ());

	

	return 0;
}