#pragma once

#include <vector>

#include "Camera.h"
#include "Mesh.h"
#include "LightSource.h"
#include "Material.h"

#include "BVH.h"

class Scene {
public:
	inline Scene () {}
	virtual ~Scene() {}

	inline const Camera& camera() const { return m_camera; }

	inline Camera& camera() { return m_camera; }
	
	inline const std::vector<Mesh> & meshes () const { return m_meshes;  }
	
	inline std::vector<Mesh> & meshes () { return m_meshes; }

	inline const std::vector<LightSource>& lightsrc() const{return m_lightsrc;}
	inline std::vector<LightSource>& lightsrc() {return m_lightsrc;}

	inline const  Material& material() const{return m_material;}

	inline  Material& material() {return m_material;}

	inline  BVH& bvh() {return m_bvh;}

private:
	Camera m_camera;
	std::vector<Mesh> m_meshes;
	std::vector<LightSource> m_lightsrc;
	
	Material m_material;
	BVH m_bvh;
};