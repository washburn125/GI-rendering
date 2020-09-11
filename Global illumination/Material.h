#include "Vec3.h"
#include "Ray.h"

class Material {
public:
	LightSource (const Vec3f pos, const Vec3f color,
			float intensity ) : 
		m_intensity =  intensity;
		m_position = pos;
		m_color = color;
	}

	/// Generate a ray for a given (u,v) coordinate on the image plane.
	Ray rayAt (float u, float v) const {
		return Ray (m_position, normalize (m_lowerLeftCorner + u * m_horizontal + v * m_vertical - m_position));
	}

private:
	float m_intensity;
	
};