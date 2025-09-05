# System reliability analysis using linear programming
System reliability analysis of general sets can be done reformulating the calculation of failure probability in terms of linear programming. The bounds of the system failure probability is calculated by solving an optimization problem. These Matlab codes developed here uses the cvx (http://cvxr.com/cvx/) toolbox to solve the convex optimization problem. The bounds are computed by solving a linear programming (LP) problem. 

The description of the examples are provided in the paper: 
- Song, J., & Kiureghian, A. der. (2003). Bounds on System Reliability by Linear Programming. Journal of Engineering Mechanics, 129(6), 627–636. https://doi.org/10.1061/~ASCE!0733-9399~2003!129:6~627!
