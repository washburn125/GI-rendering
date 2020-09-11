#pragma once

#include <random>
#include <cmath>
#include <algorithm>
#include <limits>

#include "Vec3.h"
#include "Image.h"
#include "Camera.h"
#include "Scene.h"
#include "AABB.h"
#include "phm.h"

std::default_random_engine ran;
/*

*/

using namespace std;

class RayTracer {
public:
	RayTracer () {}
	virtual ~RayTracer() {}
	   
	inline bool rayTrace (const Ray & ray, 
						  const Scene & scene, 
						  size_t & meshIndex, 
						  size_t & triangleIndex, 
						  float & u, 
						  float & v, 
						  float & d) {
		const auto& meshes = scene.meshes();
		float closest = std::numeric_limits<float>::max();
		bool intersectionFound = false;
		for (size_t mIndex = 0; mIndex < meshes.size(); mIndex++) {
			const auto& P = meshes[mIndex].vertexPositions();
			const auto& T = meshes[mIndex].indexedTriangles();
			for (size_t tIndex = 0; tIndex < T.size(); tIndex++) {
				const Vec3i& triangle = T[tIndex];
				float ut, vt, dt;
				if (ray.triangleIntersect(P[triangle[0]], P[triangle[1]], P[triangle[2]], ut, vt, dt) == true) {
					if (dt > 0.f && dt < closest) {
						intersectionFound = true;
						closest = dt;
						meshIndex = mIndex;
						triangleIndex = tIndex;
						u = ut;
						v = vt;
						d = dt;
					}
				}
			}
		}
		return intersectionFound;
	}


	inline Vec3f shade(const Scene & scene, size_t meshIndex, size_t triangleIndex, float u, float v) {
		uniform_real_distribution<float> unif(0.f,1.f);
		//cout << meshIndex << endl;
		const auto& mesh = scene.meshes()[meshIndex];
		const auto& P = mesh.vertexPositions();
		const auto& N = mesh.vertexNormals();
		const Vec3i& triangle = mesh.indexedTriangles()[triangleIndex];
		Vec3f hitNormal = normalize((1.f - u - v) * N[triangle[0]] + u * N[triangle[1]] + v * N[triangle[2]]);
		
		Vec3f center_point = (1.f - u - v) * P[triangle[0]] + u * P[triangle[1]] + v * P[triangle[2]];
		int num_ls = scene.lightsrc().size();
		Vec3f cum_color = Vec3f();
		for (int k = 0 ; k < num_ls; k++){
			
			Vec3f light_pos = scene.lightsrc()[k].randSourcePoint(unif(ran),unif(ran));
			Vec3f color =  scene.material().evaluateColorResponse(hitNormal, normalize(light_pos - center_point), scene.camera().getPosition() - center_point, 0.9);
			Vec3f light_color = scene.lightsrc()[k].getColor();
			for (int i = 0; i < 3; i++ ){
				color[i] =5* color[i]*light_color[i]/0.1/length(light_pos - center_point)/length(light_pos - center_point);  //constant are found experimentally
				if(color[i] > 1)
					color[i] = 1.f;
				if (color[i] < 0)
					color[i] = 0.f;
			
			}

			Ray ray(center_point,  light_pos);
			float u2, v2, d2;
			bool inter = rayTrace(ray, scene,  meshIndex, triangleIndex, u2, v2, d2) ;
			if (inter) {
				color = Vec3f(0.f, 0.f, 0.f);
			}

			
			cum_color += color;

		}


		return  cum_color;
	}


		//shading routine using photon maps

		inline Vec3f shade_map(const Scene & scene, size_t meshIndex, size_t triangleIndex, float u, float v, 
		KDTree& tree, vector<Photon>& photons) {
			
			uniform_real_distribution<float> unif(0.f,1.f);
			//cout << meshIndex << endl;
			const auto& mesh = scene.meshes()[meshIndex];
			const auto& P = mesh.vertexPositions();
			const auto& N = mesh.vertexNormals();
			const Vec3i& triangle = mesh.indexedTriangles()[triangleIndex];
			Vec3f hitNormal = normalize((1.f - u - v) * N[triangle[0]] + u * N[triangle[1]] + v * N[triangle[2]]);
			Vec3f center_point = (1.f - u - v) * P[triangle[0]] + u * P[triangle[1]] + v * P[triangle[2]];
			
			Vec3f cum_color = Vec3f();
			double rad = 0.5;

			vector<size_t> neighbours = tree.neighborhood_indices(center_point, rad);
			int num_neigh = neighbours.size();
			float cum_intense = 0;
			Vec3f av_dir = Vec3f();
			for (int i = 0 ; i < num_neigh; i++){
				//cum_intense += max(0.0, 1- dist(center_point, photons[neighbours[i]].getPos())/rad);
				cum_intense += photons[neighbours[i]].getPower();
				av_dir += photons[neighbours[i]].getDir();
			}

			av_dir = av_dir/num_neigh;

			cum_intense /= M_PI* rad*rad;
			cum_intense/= photons.size();
			cum_intense *= 5000;
			

			Vec3f radiance = Vec3f(1.f, 1.f, 1.f) * cum_intense; 
			Vec3f res  = scene.material().evaluateColorResponse(hitNormal, normalize(av_dir), 
				scene.camera().getPosition() - center_point, 0.9); 
			res *= radiance;

			for (int i = 0; i < 3; i++ ){
				if(res[i] > 1)
					res[i] = 1.f;
				if (res[i] < 0)
					res[i] = 0.f;
			
			}
			return  res;
	}

	inline Vec3f rand_ray_hemisphere(Vec3f direction, int ray_ind, int ray_num){
		 uniform_real_distribution<float> unif(0.f,1.f);
		 Vec3f b3 = normalize(direction);
		 Vec3f different = (std::abs(b3[0]) < 0.5f) ? Vec3f(1.0f, 0.0f, 0.0f) : Vec3f(0.0f, 1.0f, 0.0f);
		 Vec3f b1 = normalize(cross(b3, different));
		 Vec3f b2 = cross(b1, b3);

		 
		 float z = unif(ran);
		 float r = std::sqrt(1.0f - z * z);
		 //stratified sampling
		 float theta = (- 0.5f + (ray_ind + unif(ran))/ray_num)* 2* 3.14;
		 //cout << "ray ind is " << ray_ind << " and ray_num is " << ray_num << " so theta is " << theta /2/3.14 << endl;
		 float x = r * std::cos(theta);
		 float y = r * std::sin(theta);
		 
		 // Construct the vector that has coordinates (x,y,z) in the basis formed by b1, b2, b3
		 return x * b1 + y * b2 + z * b3;

	}

 
 	inline Vec3f path_trace(int level, Ray ray, const Scene& scene, bool& inter,int ray_ind, 
 		int ray_num, KDTree& tree, vector<Photon>& photons){
 		
		size_t meshIndex, triangleIndex;
		float u, v, d;
		bool intersectionFound = rayTrace (ray, scene, meshIndex, triangleIndex, u, v, d);

		if (intersectionFound && d > 0.f){
			inter = true;
			Vec3f res = shade_map(scene, meshIndex, triangleIndex, u, v, tree, photons); //shading with photon maps
			//Vec3f res = shade(scene, meshIndex, triangleIndex, u, v); //uncomment  to use a shading without photon maps
			if (level == 3)
				return res;
 
		const auto& P = scene.meshes()[meshIndex].vertexPositions();
		const auto& N = scene.meshes()[meshIndex].vertexNormals();
		const Vec3i& triangle = scene.meshes()[meshIndex].indexedTriangles()[triangleIndex];
		Vec3f hitNormal = normalize((1.f - u - v) * N[triangle[0]] + u * N[triangle[1]] + v * N[triangle[2]]);
		Vec3f center_point = normalize((1.f - u - v) * P[triangle[0]] + u * P[triangle[1]] + v * P[triangle[2]]);
		int m_level = level;
		
		return res + path_trace(level+1, Ray(center_point, rand_ray_hemisphere(hitNormal, ray_ind, ray_num)), 
		scene, inter, ray_ind, ray_num,tree, photons)/(m_level +1)/2;
		}
		else 
			return Vec3f();

 	}



 	inline void render_MC (const Scene& scene, Image& image) {

		uniform_real_distribution<float> unif(0.f,1.f);
		int num_rays =10;
		size_t w = image.width();
		size_t h = image.height();
		const Camera& camera = scene.camera();

		vector<Photon> photons;
		createPhotonMap(scene, photons);
		KDTree tree(photons);
		cout << "   photon map built" << endl;
		for (int y = 0; y < h; y++) {
			static int progress = 0;
			progress++;
			if (progress % 10 == 0)
				std::cout << ".";
#pragma omp parallel for
			for (int x = 0; x < w; x++) {
				Vec3f colorResponse;
				int num_inter = 0;
				Vec3f cum_color(0, 0, 0);

					for (int i = 0; i < num_rays; i++){
						Ray ray = camera.rayAt ((float (x) + unif(ran) ) / w, 1.f - (float (y) +unif(ran) )/ h);
						bool inter = false;
						Vec3f res = path_trace(0, ray, scene, inter, i, num_rays,tree, photons);
						if (inter){
							num_inter++;
							cum_color += res;
						}
					}
					if (num_inter > 0){
						Vec3f av_color = cum_color/num_rays;
						image (x,y)   = image(x,y)*(num_rays -num_inter )/num_rays + av_color;
					}	
			}
		}
	}


	


	inline void render (const Scene& scene, Image& image) {

		uniform_real_distribution<float> unif(0.f,1.f);

		int num_rays =5;
		size_t w = image.width();
		size_t h = image.height();
		const Camera& camera = scene.camera();
		for (int y = 0; y < h; y++) {
			static int progress = 0;
			progress++;
			if (progress % 10 == 0)
				std::cout << ".";
#pragma omp parallel for
			for (int x = 0; x < w; x++) {
				Vec3f colorResponse;
				int num_inter = 0;
				Vec3f cum_color(0, 0, 0);

					for (int i = 0; i < num_rays; i++){
						Ray ray = camera.rayAt ((float (x) + unif(ran) ) / w, 1.f - (float (y) +unif(ran) )/ h);
						size_t meshIndex, triangleIndex;
						float u, v, d;
						bool intersectionFound = rayTrace (ray, scene, meshIndex, triangleIndex, u, v, d);
						if (intersectionFound && d > 0.f){
							num_inter++;
							Vec3f res = shade(scene, meshIndex, triangleIndex, u, v);
							cum_color += res;
							
						}
					}
					if (num_inter > 0){
						Vec3f av_color = cum_color/num_rays;
						image (x,y)   = image(x,y)*(num_rays -num_inter )/num_rays + av_color;
						//if (image (x,y)[0] >=1 || image (x,y)[1]>= 1 || image (x,y)[2] >= 1 )
						//	cout << "attention, color value out of boders!  " << image(x,y) << endl;
					}	
			}
		}
	}

	void tracePhoton(Photon photon, const Scene& scene, vector<Photon>& photons){
		uniform_real_distribution<float> unif(0.f,1.f);
		size_t meshIndex, triangleIndex;
		float u, v, d;
		Ray ray(photon.getPos(), photon.getDir());
		bool intersectionFound = rayTrace (ray, scene, meshIndex, triangleIndex, u, v, d);
		if (intersectionFound){
			const auto& mesh = scene.meshes()[meshIndex];
			const auto& P = mesh.vertexPositions();
			const auto& N = mesh.vertexNormals();
			const Vec3i& triangle = mesh.indexedTriangles()[triangleIndex];
			Vec3f hitNormal = normalize((1.f - u - v) * N[triangle[0]] + u * N[triangle[1]] + v * N[triangle[2]]);
			Vec3f center_point = (1.f - u - v) * P[triangle[0]] + u * P[triangle[1]] + v * P[triangle[2]];

			photon.setPos(center_point);
			photons.push_back(photon);

			if (unif(ran) > 0.3 && photon.getPower() > 0.01){
					Vec3f ref_ray = ray.direction() - 2* dot(ray.direction(), hitNormal)*hitNormal;
					tracePhoton(Photon(center_point, photon.getPower()/2, ref_ray), scene, photons);
			}
		}
	}



	void createPhotonMap(const Scene& scene, vector<Photon>& photons){
		uniform_real_distribution<float> unif(0.f,1.f);
		int num_ls = scene.lightsrc().size();
		for (int i = 0; i < num_ls; i++){
			Vec3f light_pos = scene.lightsrc()[i].randSourcePoint(unif(ran),unif(ran));
			Vec3f norm = normalize(light_pos)*(-1);

			for (int j = 0; j < 500000 ; j++){
				Vec3f ray = rand_ray_hemisphere(norm, 0, 1);
				Photon photon(light_pos, 1.0, ray);
				tracePhoton(photon, scene, photons);
			}
		}
		//writing in file
		ofstream myfile;
 		myfile.open ("../map2.txt");
		myfile << photons.size() << endl;
		for (int i = 0; i < photons.size(); i++){
			Vec3f coord= photons[i].getPos();
			myfile << coord[0] << " "<< coord[1]<<" " << coord[2] << endl;
		}
		myfile.close();
	}

};

