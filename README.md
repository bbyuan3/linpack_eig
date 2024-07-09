# linpack_eig
linpack eig matlab version

when we need a eig calculation, linpack eig is good reference

i transfer it into the matlab (maybe C another day...) version

dgeev : compute all eigenvalues and, optionally, 
the right and/or left eigenvectors of a real, nonsymmetric matrix.

DEGBAL: used to balance a general real matrix A

DGEHD2: reduces a general real matrix A to upper Hessenberg form H using 
an orthogonal similarity transformation: A = Q * H * Q^T. This function is
typically used as a preprocessing step for computing eigenvalues and eigenvectors.

DLAHQR:
 used to compute the eigenvalues and, optionally, the Schur factorization of 
 an upper Hessenberg matrix. This function is typically used after a matrix 
 has been reduced to upper Hessenberg form by DGEHRD.