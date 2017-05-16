# Lie_Algebra
This project allows for calculation of irreducible representations and tensor product decomposition of Lie Algebras.  So far the classical groups: A, B, C, and D are included. The formalism of the bases follow those found in [1]. A lot of the algorithms follow the usage in [2]. A much more thorough walkthrough of the topics can be found in [3].

## Getting Started

The use of this header is rather simple. Only two additional libraries are needed.

### Prerequisite Libraries
Eigen Library found at http://eigen.tuxfamily.org  This handles most of the linear algebra. The pseudoinverse of a matrix needed to be added.

Boost Library http://www.boost.org/ Specifically the math/special_functions/round.hpp because certain bases were not integer valued.

## Running Lie_Algebra
To choose a group's algebra simply initialize the class with the desired dimension.

Example:
```
typedef typename Eigen::VectorXd vec;
A_n A(2); // Initialize A_2

A.print(A.funamental_weights_omega); // Prints the fundamental weights in the dynkin (omega) basis 

vec v(2,0);

A.print(A.weight_tower(v)); //Prints all the states of the weight (2,0) in a linear fashion

A.print(A.tower_layout(A.weight_tower(v))); //Prints the states in a 2d format that has all states that are the same distance from highest weight on the same line

vec w(2,0);

A.print(A.TensorProductDecomp_weights(v,w));//Prints out the irreducible tensor sum of the tensor product in weight form

std::cout<<A.dim(v)<<std::endl; //Prints out the dimenions of (2,0) = 6
```

## Sources
[1] https://www.math.stonybrook.edu/~kirillov/mat552/liegroups.pdf

[2] https://arxiv.org/abs/1206.6379

[3] Georgi, Howard (1999) Lie Algebras in Particle Physics. Reading, Massachusetts: Perseus Books. ISBN 0-7382-0233-9.

[4] R. Slansky, Group theory for unified model building

##Authors:

* **Nathan Papapietro** 

## Current Update

Qt Released

## Next Update
The program is at a point where I feel nothing more to add at the moment. I will begin to review it and rewrite algorithms for greater effecicency
