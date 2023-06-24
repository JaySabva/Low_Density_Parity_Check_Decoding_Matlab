# Low_Density_Parity_Check Decoding using Tanner Graph in MATLAB

LDPC (Low-Density Parity-Check) decoding is an error-correction technique widely used in communication systems. This implementation focuses on LDPC decoding using a Tanner graph representation and is implemented in MATLAB.

The MATLAB code utilizes a given H-matrix to construct the Tanner graph. The degree of check nodes and variable nodes is calculated based on the H-matrix. The code performs Monte Carlo simulations to evaluate the success probability of decoding LDPC codes over a Binary Erasure Channel (BEC). It also provides an analysis of algorithmic convergence.

The LDPC decoding process involves passing messages between variable nodes and check nodes. At each iteration, variable nodes update their beliefs based on received erasure or codeword values. Check nodes then update their estimates based on the information received from variable nodes. The process continues iteratively until a stopping criterion, such as the maximum number of iterations, is reached.

The code plots the probability of successful decoding as a function of the error probability of the BEC. It also generates a graph illustrating the algorithmic convergence, showing the error probability over iterations for different error probabilities. Theoretical convergence curves are also provided for comparison.

This MATLAB implementation provides a useful tool for studying LDPC decoding performance and understanding the impact of different parameters on the decoding process. It can serve as a valuable resource for researchers and practitioners working with LDPC codes and communication systems.

# GitHub Repository:
