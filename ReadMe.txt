Learnin Koopman Operator

Ref: M. Khosravi, "Representer Theorem for Learning Koopman Operators," in IEEE Transactions on Automatic Control, doi: 10.1109/TAC.2023.3242325.
arXiv link: https://arxiv.org/abs/2208.01681

The repository contains example of learning Koopman operator for the EDMD method and 
the regularized learning approaches discussed in the reference above, i.e., the following 
learning problem 
   min_K   E(K) + lambda * R(K),
   s.t.    K in C,
where E(.) is the sum squared error loss function, C is a given constraint, 
R(.) is the regularization term. 

Mohammad Khosravi
Email: mohammad.khosravi@tudelft.nl
Delft Center for Systems and Control (DCSC)
Delft University of Technology (TU Delft) 
August 2022
