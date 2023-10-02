# Solving-Linear-Systems-Study-NA
A two-part study experiment on solving linear pentadiagonal systems. We make a comparison between Cramer's algorithm, gaussian elimination, PTRANSII, and Matlab's algorithm for solving linear systems. 
In the first study experiment, we compare the execution time of these 4 algorithms for solving pentadiagonal systems of dimensions 4x4 up to 35x35. 
For each dimension, we take as a sample the average resolution time of 10 random systems and then we plot the results of every algorithm in a Time=F(dimension) figure.
In the second study experiment, we only compare the "champions" of the first experiment whose difference in performance is not clear enough. These are Matlab's algorithm and PTRANSII. So to acknowledge a clear difference in their performances we take samples for dimensions 4x4 up to 2000x2000 with a 50-step increment (e.g. 4x4 - 54x54, etc.). The samples for each dimension are taken in the same way as before.
