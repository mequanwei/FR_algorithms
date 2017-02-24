This code is for our CVPR2011 paper: 

Meng Yang, Lei Zhang, Jian Yang, and David Zhang. Robust Sparse Coding for Face Recognition. CVPR 2011.

Copyright: Meng Yang, Lei Zhang. BRC, PolyU, Hong Kong.
Contact:   csmyang@comp.polyu.edu.hk;cslzhang@comp.polyu.edu.hk

For simplity, in the demo, we fix iter=1 in step 4 of IRLS. We use the maximal iterative number
to stop iteration.

For the sparse coding, we used l1_ls toolbox to solve it. l1_ls toox is not the fastest tool to solve l1-norm
minimization, but it often has stable and good performance.

For the database:
In the paper, we used Extended YaleB, AR, MPIE databases.

Extended Yale B could be download in the Homepage of Extended Yale B. You should also divide the database 
into subset1,subset2 and subset 3 when you do the experiment with random pixel corruption and random block occlusion.
(How to divide the database into the three subsets can also be found in the homepage of Extended Yale B.) 

For the AR dataset used in the paper, you could download it at 'http://www4.comp.polyu.edu.hk/~csmyang/Publication.html' 
through the link to our ECCV paper. Once you download it, you could directly run Demo_RSC_AR_disguise.m and Demo_RSC_AR_disguise2.m.

For MPIE, you need to get it from CMU.


If you have any question about my paper and my code, welcome to contact us.


Thanks.


YANG Meng
Ph.D candidate
epartment of computing
HK PolyU.