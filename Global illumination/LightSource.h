#pragma once

#include "Vec3.h"
#include "Ray.h"

class LightSource {
public:
	LightSource() {} 
	LightSource (const Vec3f pos, const Vec3f color,
			float intensity ) {
		m_intensity =  intensity;
		m_position = pos;
		m_color = color;
	}

	Vec3f getPosition() const {
		return m_position;
	}

	Vec3f get_hor() const{

		return m_horizontal;
	}

	Vec3f get_ver()const{
		return m_vertical;
	}

	Vec3f getColor() const {
		return m_color;
	}

	Vec3f init_area(){
		coord.push_back(Vec3f(-0.1, 2.0, 0.2));
		coord.insert(coord.end(), {Vec3f(0.1, 2.0, 0.2),Vec3f(-0.1, 2.0, -0.1),Vec3f(0.1, 2.0, -0.1)});
		
	}

	void initAreaSource ( const Vec3f& up = Vec3f(0.f, 1.f, 0.f),  float size = 0.2){
		 	//size = 0.2

		Vec3f lookAt = Vec3f();
		Vec3f point = Vec3f(m_position[0], m_position[1], 0);

		float halfHeight = size /2;
		float halfWidth =size/2;

		Vec3f atFrom = normalize(m_position);

		Vec3f u = normalize(cross(up, atFrom));

		Vec3f v = normalize( cross(atFrom, u));

		m_lowerLeftCorner = m_position - halfWidth * u - halfHeight * v ; //- atFrom; - wtf, we dont<t need it 
		m_horizontal = 2.f * halfWidth * u;
		m_vertical = 2.f * halfHeight * v;
	}

	Vec3f randSourcePoint(float u, float v) const {
		return Vec3f(m_lowerLeftCorner + u * m_horizontal + v * m_vertical ); //\ position?
	}
	
	/// Generate a ray for a given (u,v) coordinate on the image plane.

private:
	float m_intensity;
	Vec3f m_position;
	Vec3f m_color;

	std::vector<Vec3f> coord;
	Vec3f m_lowerLeftCorner;
	Vec3f m_horizontal;
	Vec3f m_vertical;



};