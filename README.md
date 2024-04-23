# CE-42
 Implementation of CE-43 paper (in Matlab) - THE CODES ARE NOT MINE


THIS IS THE ORIGINAL README FROM THE AUTHORS:  

This simulation code package is mainly used to reproduce the results of the following paper [1]:

[1] M. Cui and L. Dai, “Channel Estimation for Extremely Large-Scale MIMO: Far-Field or Near-Field?,” IEEE Transactions on Communications, vol. 70, no. 4, pp. 2663-2677, April 2022.

*********************************************************************************************************************************
If you use this simulation code package in any way, please cite the original paper [1] above. 
 
The author in charge of this simulation code pacakge is: Mingyao Cui (email: cmy20@mails.tsinghua.edu.cn).

Reference: We highly respect reproducible research, so we try to provide the simulation codes for our published papers 
( more information can be found at: 
http://oa.ee.tsinghua.edu.cn/dailinglong/publications/publications.html )

Please note that the MATLAB R2020b is used for this simulation code package,  
and there may be some imcompatibility problems among different MATLAB versions. 

Copyright reserved by the Broadband Communications and Signal Processing Laboratory (led by Dr. Linglong Dai), 
Beijing National Research Center for Information Science and Technology (BNRist), 
Department of Electronic Engineering, Tsinghua University, Beijing 100084, China. 

*********************************************************************************************************************************
Abstract of the paper: 

Extremely large-scale multiple-input-multiple-output (XL-MIMO) is promising to meet the high rate requirements for future 6G. 
To realize efficient precoding, accurate channel state information is essential. 
Existing channel estimation algorithms with low pilot overhead heavily rely on the channel sparsity in the angular domain, 
which is achieved by the classical far-field planar-wavefront assumption. 
However, due to the non-negligible near-field spherical-wavefront property in XL-MIMO, 
this channel sparsity in the angular domain is not achievable. 
Therefore, existing far-field channel estimation schemes will suffer from severe performance loss. 
To address this problem, in this paper, we study the near-field channel estimation by exploiting the polar-domain sparsity. 
Specifically, unlike the classical angular-domain representation that only considers the angular information, 
we propose a polar-domain representation, which simultaneously accounts for both the angular and distance information. 
In this way, the near-field channel also exhibits sparsity in the polar domain, based on which, 
we propose on-grid and off-grid near-field XL-MIMO channel estimation schemes. 
Firstly, an on-grid polar-domain simultaneous orthogonal matching pursuit (P-SOMP) algorithm is proposed to efficiently estimate the near-field channel. 
Furthermore, an off-grid polar-domain simultaneous iterative gridless weighted (P-SIGW) algorithm is proposed to improve the estimation accuracy. 
Finally, simulations are provided to verify the effectiveness of our schemes.
*********************************************************************************************************************************
How to use this simulation code package?

The simulation results an be obtained by the running the file 'main_NMSE_vs_distance.m'.

*********************************************************************************************************************************
Enjoy the reproducible research!