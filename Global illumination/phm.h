#pragma once


#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

#include "Vec3.h"
#include "Ray.h"


/// The majotrity of the code for implemetation of kd-tree was taken from https://github.com/crvs/KDTree
/// and then adapted to the project

class Photon {
public:
	Photon(Vec3f pos, float power, Vec3f dir, short flag){
		m_pos = pos;
		m_p = power;  //correspond to intensity
		m_dir = dir;
		m_flag = flag;
	}

	Photon(Vec3f pos, float power, Vec3f dir){
		m_pos = pos;
		m_p = power; 
		m_dir = dir;
		
	}

	Vec3f getPos() { return m_pos; }
	Vec3f getDir() {return m_dir;}
	float getPower(){return m_p;}
	void setPos(Vec3f pos){m_pos = pos;}
	void setDir(Vec3f dir){m_dir = dir;}
	void setPower(float power){m_p = power;}

	

private:
    Vec3f m_pos;       
 	float m_p;           
    Vec3f m_dir;
    short m_flag;         
};




/*
 * file: KDTree.hpp
 * author: J. Frederico Carvalho
 *
 * This is an adaptation of the KD-tree implementation in rosetta code
 *  https://rosettacode.org/wiki/K-d_tree
 * It is a reimplementation of the C code using C++.
 * It also includes a few more queries than the original
 *
 */


//using point_t = std::vector< double >;
using indexArr = std::vector< size_t >;
using pointIndex = typename std::pair< Vec3f, size_t >;

class KDNode {
   public:
    using KDNodePtr = std::shared_ptr< KDNode >;
    size_t index;
    Vec3f x;
    KDNodePtr left;
    KDNodePtr right;

    // initializer
    KDNode(){}

    KDNode(const Vec3f &pt, const size_t &idx_, const KDNodePtr &left_,
               const KDNodePtr &right_) {
    x = pt;
    index = idx_;
    left = left_;
    right = right_;
	}

    KDNode(const pointIndex &pi, const KDNodePtr &left_,
               const KDNodePtr &right_) {
    x = pi.first;
    index = pi.second;
    left = left_;
    right = right_;
	}



    // getter
    float coord(const size_t &idx){
		return x[idx];

    }

    // conversions
   //***************** explicit operator bool() { return !(x.size()==0); } //explicit??? do we need it& can VEC3F BE EMPTY
	//operator Vec() { return x; }
	operator size_t() { return index; }
	operator pointIndex() { return pointIndex(x, index); }
};

using KDNodePtr = std::shared_ptr< KDNode >; //do we need it one more time


KDNodePtr NewKDNodePtr(){
	KDNodePtr mynode = std::make_shared< KDNode >();
    return mynode;
}

//inline float dist(const Vec3f &, const Vec3f &); //not sur we need it
inline float dist(const KDNodePtr &nd1, const KDNodePtr &nd2){
	return dist(nd1->x, nd2->x); //this pointers?

}

// Need for sorting
class comparer {
   public:
    size_t idx;
    explicit comparer(size_t idx_) : idx{idx_} {};
    inline bool compare_idx(const pointIndex &a,  //
                                  const pointIndex &b   //
	) {
    return (a.first[idx] < b.first[idx]);  //
}
};

using pointIndexArr = typename std::vector< pointIndex >;
//using pointIndexArrP = typename std::vector<std::pair< Vec3f, size_t> >;


inline void sort_on_idx(const pointIndexArr::iterator &begin,  //
                        const pointIndexArr::iterator &end,    //
                        size_t idx) {
    comparer comp(idx);
    comp.idx = idx;

    using std::placeholders::_1;
    using std::placeholders::_2;

    std::sort(begin, end, std::bind(&comparer::compare_idx, comp, _1, _2));
}

using pointVec = std::vector<Vec3f>;

class KDTree {
	public:
    KDTree(){}


    KDNodePtr root;
    KDNodePtr leaf;

    KDNodePtr make_tree(const pointIndexArr::iterator &begin,  //
                        const pointIndexArr::iterator &end,    //
                        const size_t &length,                  //
                        const size_t &level                    //
    ){
    	 if (begin == end) {
        	return NewKDNodePtr();  // empty tree
    	}

    	size_t dim = 3;//begin->first.size(); // *** just 3?

    	if (length > 1) {
        	sort_on_idx(begin, end, level);
    	}

    	auto middle = begin + (length / 2);

    	auto l_begin = begin;
		auto l_end = middle;
		auto r_begin = middle + 1;
		auto r_end = end;
		KDNodePtr left;

	    size_t l_len = length / 2;
	    size_t r_len = length - l_len - 1;

	    if (l_len > 0 && dim > 0) {
        	left = make_tree(l_begin, l_end, l_len, (level + 1) % dim); //whyyy
	    } else {
	        left = leaf;
	    }

	    KDNodePtr right;
	 
	    if (r_len > 0 && dim > 0) {
	        right = make_tree(r_begin, r_end, r_len, (level + 1) % dim);
	    } else {
	        right = leaf;
	    }


    	return std::make_shared< KDNode >(*middle, left, right);	
    }


        explicit KDTree(vector<Photon> point_array) {
		    leaf = std::make_shared< KDNode >();
		    // iterators
		    pointIndexArr arr;

		    for (size_t i = 0; i < point_array.size(); i++) {
		        arr.push_back(pointIndex(point_array[i].getPos(), i));  //point and its unique code? in unordered structure
		    }

		    auto begin = arr.begin();
		    auto end = arr.end();

		    size_t length = arr.size();
		    size_t level = 0;  // starting

		    root = KDTree::make_tree(begin, end, length, level);
	
	}

    explicit KDTree(pointVec point_array) {
    leaf = std::make_shared< KDNode >();
    // iterators
    pointIndexArr arr;

    for (size_t i = 0; i < point_array.size(); i++) {
        arr.push_back(pointIndex(point_array[i], i));  //point and its unique code? in unordered structure
    }

    auto begin = arr.begin();
    auto end = arr.end();

    size_t length = arr.size();
    size_t level = 0;  // starting

    root = KDTree::make_tree(begin, end, length, level);
	
	}


   private:
    KDNodePtr nearest_(           //
        const KDNodePtr &branch,  //
        const Vec3f &pt,        //
        const size_t &level,      //
        const KDNodePtr &best,    //
        const float &best_dist   //
    ){
    	float d, dx, dx2;

    	if (branch == NULL) { //if (!bool(*branch)) {
        	return NewKDNodePtr();  // basically, null
    	}

    	Vec3f branch_pt = branch->x;// branch_pt(*branch); //how it works&&&
    	size_t dim =3 ;// branch_pt.size();

    	d = dist(branch_pt, pt);
    	dx = branch_pt[level] - pt[level];
    	dx2 = dx * dx;

    	KDNodePtr best_l = best;
    	float best_dist_l = best_dist;

    	 if (d < best_dist) {
	        best_dist_l = d;
	        best_l = branch;
    	}


	    	size_t next_lv = (level + 1) % dim;
		    KDNodePtr section;
		    KDNodePtr other;


		 if (dx > 0) {
       		 section = branch->left;
	        other = branch->right;
	    } else {
	        section = branch->right;
	        other = branch->left;
	    }

	    KDNodePtr further = nearest_(section, pt, next_lv, best_l, best_dist_l);
	    if (further!= NULL) {
	        float dl = dist(further->x, pt);
	        if (dl < best_dist_l) {
	            best_dist_l = dl;
	            best_l = further;
	        }
	    }

		     if (dx2 < best_dist_l) {
	        further = nearest_(other, pt, next_lv, best_l, best_dist_l);
	        if (further!= NULL) {
	            float dl = dist(further->x, pt);
	            if (dl < best_dist_l) {
	                best_dist_l = dl;
	                best_l = further;
	            }
	        }
	    }

    	return best_l;


    }

    // default caller
    KDNodePtr nearest_(const Vec3f &pt){
    	size_t level  = 0;

 		float branch_dist = dist(root->x, pt);
    	return nearest_(root,          // beginning of tree
                    pt,            // point we are querying
                    level,         // start from level 0
                    root,          // best is the root
                    branch_dist);  // best_dist = branch_dist
    }


    Vec3f nearest_point(const Vec3f &pt){
    	KDNodePtr ptn =  nearest_(pt);
    	return ptn->x;
	}
    size_t nearest_index(const Vec3f &pt) {
    return size_t(*nearest_(pt));
	}


    pointIndex nearest_pointIndex(const Vec3f &pt){
    	   KDNodePtr Nearest = nearest_(pt);
    return pointIndex(Nearest->x, size_t(*Nearest));
    }


    pointIndexArr neighborhood_(  //
        const KDNodePtr &branch,  //
        const Vec3f &pt,        //
        const double &rad,        //
        const size_t &level       //
    ){
    	
    	    double d, dx, dx2;

    if (branch->x[0] == 0.0 && branch->x[1] == 0.0 && branch->x[2] == 0.0) {
        // branch has no point, means it is a leaf,
        // no points to add
       // cout << "leaf dound " << endl;
        return pointIndexArr();
    }
    


    size_t dim = 3;//pt.size();

    double r2 = rad * rad;

    d = dist(branch->x, pt);
    dx = branch->x[level] - pt[level];
    dx2 = dx * dx;

    pointIndexArr nbh, nbh_s, nbh_o;
    if (d <= r2) {
        nbh.push_back(pointIndex(*branch));
    }

    
    KDNodePtr section;
    KDNodePtr other;
    if (dx > 0) {
        section = branch->left;
        other = branch->right;
    } else {
        section = branch->right;
        other = branch->left;
    }

    nbh_s = neighborhood_(section, pt, rad, (level + 1) % dim);

    
    nbh.insert(nbh.end(), nbh_s.begin(), nbh_s.end());
    if (dx2 < r2) {
        nbh_o = neighborhood_(other, pt, rad, (level + 1) % dim);
        nbh.insert(nbh.end(), nbh_o.begin(), nbh_o.end());
    }

    return nbh;


    }

   public:
    pointIndexArr neighborhood(  //
        const Vec3f &pt,       //
        const double &rad){

    	  size_t level = 0;
    return neighborhood_(root, pt, rad, level);
    }

    pointVec neighborhood_points(  //
        const Vec3f &pt,         //
        const double &rad){

    	size_t level = 0;
    pointIndexArr nbh = neighborhood_(root, pt, rad, level);
    pointVec nbhp;
    nbhp.resize(nbh.size());
    std::transform(nbh.begin(), nbh.end(), nbhp.begin(),
                   [](pointIndex x) { return x.first; });
    return nbhp;
    }

    indexArr neighborhood_indices(  //
        const Vec3f &pt,          //
        const double &rad){
    	//cout << "in neighborhood_indices  " << endl;
    	 size_t level = 0;
    pointIndexArr nbh = neighborhood_(root, pt, rad, level);
   // cout << "neigh found " << endl;
    indexArr nbhi;
    nbhi.resize(nbh.size());
    std::transform(nbh.begin(), nbh.end(), nbhi.begin(),
                   [](pointIndex x) { return x.second; });
    return nbhi;

    }
};