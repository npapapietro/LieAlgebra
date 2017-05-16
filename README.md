# Lie_Algebra
This project allows for calculation of irreducible representations and tensor product decomposition of Lie Algebras.  So far the classical groups: A, B, C, and D are included. The formalism of the bases follow those found in [1]. A lot of the algorithms follow the usage in [2]. A much more thorough walkthrough of the topics can be found in [3].

## Getting Started

The use of this header is rather simple. Only two additional libraries are needed.

### Prerequisite Libraries
Eigen Library found at http://eigen.tuxfamily.org  This handles most of the linear algebra. Several helper functions are used to handle rational numbers (finding the inverse of a matrix for example).

Boost Library http://www.boost.org/ Specifically the Rational Number Library.

## Running Lie_Algebra
To choose a group's algebra simply initialize the class with the desired dimension.

Example:
```
#include "RootStructure.h"


using namespace Eigen;
using namespace Representation;
using std::cout;
using std::endl;
int main() {
	
	LieBase<A> G(2);
	
	Eigen::Matrix<boost::rational<int>,2, 1> x,y;

	x << boost::rational<int>(2), boost::rational<int>(0);
	y << boost::rational<int>(1), boost::rational<int>(0);

	auto sum = G.tensorProductDecomp(weight(x), weight(y));

	//Tensor product between (2,0) and (1,0)
	cout << " (2,0)x(1,0) = ";
	for (size_t i = 0; i<sum.size(); i++)
	{
		cout << int_cast(sum[i].omega).transpose();
		if (i != sum.size() - 1)
			cout << " + ";
	}cout << endl;
	
	return 0
  }
```

## Sources
[1] https://www.math.stonybrook.edu/~kirillov/mat552/liegroups.pdf

[2] https://arxiv.org/abs/1206.6379

[3] Georgi, Howard (1999) Lie Algebras in Particle Physics. Reading, Massachusetts: Perseus Books. ISBN 0-7382-0233-9.

[4] R. Slansky, Group theory for unified model building

##Authors:

* **Nathan Papapietro** 

## Current Update

V2 

## Next Update
Redo the Qt application with the new backend. Possibly add exceptional groups
