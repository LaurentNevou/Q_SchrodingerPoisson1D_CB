# Q_SchrodingerPoisson1D_CB
Full Schrodinger-Poisson solver in 1D in the conduction band

![Results_InGaAs-AlInAs_300K](https://user-images.githubusercontent.com/35040499/147890985-b855b9c5-e2d3-4f99-a411-0dfa9527cf16.PNG)

This program solves the Schrodinger-Poisson equations in the conduction band for any heterostructures.

2 versions are available:
1) one is using the Kane model to take into account the non-parabolicity. Two algorithm are available, the shooting method and the diagonalisation of the Hamiltonian (FEM). Both algorithm lead to the same results but the FEM is faster. 
2) The second is only using the shooting method and the non-parabolicity is implemented via the alpha parameter, alpha=1/Egap; and meff(E)=meff(0)*(1+alpha*E)

Both models follow the book of Paul Harrison and actually makes the 2d density of states not constant
https://onlinelibrary.wiley.com/doi/book/10.1002/9781118923337

DOI:10.1002/9781118923337

A strain model is included. It basically shifts the conduction band edge
The strain is mainly interesting for InGaAs/GaAs heterostructures
The non-parabolicity is also included into the density of states for the Poisson solver.

-> Additionnal material can be added in the "materialDB_ZB.csv" file

-> II-VI and cubic nitride material parameters are available but should be grabt in the "Library.m" file

Enjoy! If you like it, don t forget the star!
