/*
LI Yunsheng, 2015.04.20
RRT algorithm implementation
*/


#ifndef RRT_H
#define RRT_H

#include <CCStateSpace.h>

class Node {
public:
	CCStateSpace::State state;
	int spiltdim;
	double spiltval;
	Node* left;
	Node* right;
	// Node* parent; // parent in kd-tree
	Node* parent_rrt; // parent in RRTree

	// constructor
	Node(const CCStateSpace::State& s) :state(s), spiltdim(0), spiltval(s[0]), left(nullptr), right(nullptr), parent_rrt(nullptr){}
	//
	~Node();
};

// Rapidly Random Exploring Tree. Its internal data structure is like kd-tree.
class RRTree {
public:
	Node* root;

	// constructor
	RRTree() :root(nullptr) {}
	RRTree(CCStateSpace::State * s): root(new Node(*s)) {}
	// destructor
	~RRTree();

	// member functions
    // insert a Node
	Node* insert(const CCStateSpace::State& s);
	//Node* insert(Node* n);
	// the nearst node to s in the tree
	Node* nearest(const CCStateSpace::State& s);
	Node* nearest(Node* n);
	
};

// build rrt from start to goal configuration
RRTree* buildRRTree(CCStateSpace::State* start, CCStateSpace::State* goal, double** costmap, double critical_val);

#endif
