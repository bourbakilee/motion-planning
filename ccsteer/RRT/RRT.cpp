#include "RRT.h"
#include "evaluation.h"
#include <cstdlib>

Node::~Node()
{
	this->left = nullptr;
	this->right = nullptr;
	//this->parent = nullptr;
	this->parent_rrt = nullptr;
}

void delete_node(Node* n)
{
	if (n == nullptr)
	{
		return;
	}
	else if (n->left == nullptr && n->right == nullptr)
	{
		delete n;
		return;
	}
	else
	{
		if (n->left != nullptr)
			delete_node(n->left);
		if (n->right != nullptr)
			delete_node(n->right);
	}
}

RRTree::~RRTree()
{
	delete_node(this->root);
}

Node* recursive_insert(const CCStateSpace::State & s, Node* node, int spiltdim)
{
	if (node == nullptr)
	{
		node = new Node(s);
	}
	else if(s[spiltdim]<node->spiltval)
	{
		node->left = recursive_insert(s, node->left, (spiltdim + 1) % 3);
	}
	else
	{
		node->right = recursive_insert(s, node->right, (spiltdim + 1) % 3);
	}
	return node;
}

Node* RRTree::insert(const CCStateSpace::State & s)
{
	return recursive_insert(s, this->root, 0);
}


void nearest_neighbor(const CCStateSpace::State& s, Node* node, int spiltdim, double min_dist, Node* nn)
{
	if (node == nullptr)
		return;
	double distance = CCStateSpace::distance(s, node->state, CCStateSpace::CCStateSpace());
	if (distance < min_dist)
	{
		nn = node;
		min_dist = distance;
	}
	if (s[spiltdim] < node->spiltval)
	{
		nearest_neighbor(s, node->left, (spiltdim + 1) % 3, min_dist, nn);
		nearest_neighbor(s, node->right, (spiltdim + 1) % 3, min_dist, nn);
	}
	else
	{
		nearest_neighbor(s, node->right, (spiltdim + 1) % 3, min_dist, nn);
		nearest_neighbor(s, node->left, (spiltdim + 1) % 3, min_dist, nn);
	}
}

Node * RRTree::nearest(const CCStateSpace::State & s)
{
	Node* nn = nullptr;
	double min_dist = CCStateSpace::distance(s, this->root->state, CCStateSpace::CCStateSpace());
	nearest_neighbor(s, this->root, 0, min_dist, nn);
	return nn;
}

Node * RRTree::nearest(Node * n)
{
	return this->nearest(n->state);
}

// generate random number
inline double random(double start, double end)
{
	return start + (end - start)*rand() / (RAND_MAX + 1.0);
}

// build a rrt from start state to goal state
RRTree* buildRRTree(CCStateSpace::State * start, CCStateSpace::State * goal, double** costmap, double critical_val)
{
	RRTree* rrt = new RRTree(start);
	Node* goal_node = new Node(*goal);
	const CCStateSpace space = CCStateSpace();
	Node *rand, *near;
	int i = 0;
	for (;;)
	{
		// sampling region:[0,100]X[0,100]X[0,2pi]
		if (i % 10 == 0)
		{
			// goal-biasing, 10% probility
			rand = goal_node;
		}
		else
		{
			rand = new Node(CCStateSpace::State(space, random(0., 100.), random(0., 100.), random(0., two_pi)));
		}
		i++;
		// the nearest node in the tree to the sampling node
		near = rrt->nearest(rand);
		CCStateSpace::CCPath path(near->state, rand->state, space);
		// check if the computed path has too high cost
		bool res = evaluation(path, costmap, critical_val);
		if (res)
		{
			Node* tmp = rrt->insert(rand->state);
			tmp->parent_rrt = near;
			if (rand == goal_node)
			{
				break;
			}
		}
	}
	return rrt;
}
