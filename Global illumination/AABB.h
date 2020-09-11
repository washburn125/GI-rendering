#pragma once

#include "Vec3.h"
#include "Ray.h"


using namespace std;

class AABB {
public:

	AABB() {}
	AABB(const Vec3f down, const Vec3f up ){
		m_downleft =  down;
		m_upright = up;
		
	}
/*
	vector<Vec3f> rayInterTest(const Ray &ray){
		std::vector<Vec3f> res;
		Vec3f in = Vec3f();
		Vec3f out = Vec3f();
		float tmin = - 10000000;
		float tmax = 10000000;


		for (int i = 0; i < 3; i++) {
			//float invD = 1.0f / ray.direction[a];
			float t0 = (m_downleft[i] - ray.origin[i]) / ray.direction[i];
			float t1 = (m_upright[i] - ray.origin[i]) / ray.direction[i];
			if (ray.direction[i] < 0.0f) {
				float temp = t1;
				t1 = t0;
				t0 = temp;
			}

			tmin = t0 > tmin ? t0 : tmin;
			tmax = t1 < tmax ? t1 : tmax;

			if (tmax <= tmin)
				return NULL;
		}


		in = ray.origin + ray.direction*tmin
		out = ray.origin + ray.direction*tmax


		res.push_back(in);
		res.push_back(out);

		return res;


	}
*/
	

	/// Generate a ray for a given (u,v) coordinate on the image plane.
	

private:
	Vec3f m_downleft;
	Vec3f m_upright;
};