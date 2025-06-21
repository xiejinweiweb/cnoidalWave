# cnoidalWave
This is a wave generation procedure for wavemake, see more details in the ref. (2003, A note on estimation of the Jacobian elliptic parameter in cnoidal wave theory)
By the way, in ref., there is a type error found on page 1918 of the document, after formula (8) "in which 𝑚′ is the complementary elliptic parameter, i.e. 𝑚′=m−1”. It should be 𝑚′=1−m.

Main functions:
1. The Jacobi parameter m is solved by the Newton iterative method
2. Obtain the free water surface using the Runge-Kutta 4 Integration, and solve the trajectory movement of the wavemaker for the cnoidal wave
3. visualization
