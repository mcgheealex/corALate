# Augmented Lagrangian Digital Image Correlation (2D_ALDIC)
AL-DIC(Augmented Lagrangian DIC) is a fast, parallel-computing hybrid DIC algorithm, which combines advantages of local subset DIC method (fast computation speed, and parallel computing) and finite-element-based global DIC method (guarantee global kinematic compatibility and decrease noise).  

Welcome to give the ALDIC code ratings and comments in the MATLAB File Exchange community: [![View Augmented Lagrangian Digital Image Correlation and Tracking on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/70499-augmented-lagrangian-digital-image-correlation-and-tracking)
[![DOI](https://data.caltech.edu/badge/174234539.svg)](https://data.caltech.edu/badge/latestdoi/174234539)

## Advantages of AL-DIC algorithm
* [1] It’s a fast algorithm using distributed parallel computing.  
* [2]	Global kinematic compatibility is added as a global constraint in the form of augmented Lagrangian, and solved using Alternating Direction Method of Multipliers scheme.
* [3]	Both displacement fields and affine deformation gradients are correlated at the same time.
* [4]	No need of much manual experience about choosing displacement smoothing filters.
* [5]	It works well with compressed DIC images and adaptive mesh. See our paper: Yang, J. & Bhattacharya, K. Exp Mech (2019). https://doi.org/10.1007/s11340-018-00459-y;
* [6]	Both accumulative and incremental DIC modes are implemented to deal with image sequences, which is especially quite useful for very large deformations. 
* [7]	ALDIC application example -- uniaxial compression experiment:
https://github.com/jyang526843/2D_ALDIC_v3/blob/master/Example_aldic_foam_compression_strain_eyy.gif
* [8]	ALDIC is extended with adaptive quadtree mesh to solve complex geometry. Some examples: https://uwmadison.box.com/s/4n5hmf04rzp4la96bt2rcjk4f6o5d5nf

## Prerequisites & Installation
AL-DIC MATLAB code was tested on MATLAB versions later than R2018a. Both single thread and parallel computing features are included in AL-DIC code. Please download and unzip the code to the MATLAB working path. Then, execute the mail file: main_ALDIC.m.

## Code manual 
Full size code manual is available at:
https://www.researchgate.net/publication/344796296_Augmented_Lagrangian_Digital_Image_Correlation_AL-DIC_Code_Manual


## Code demo videos
ALDIC Matlab code demo:
(Youtube) https://www.youtube.com/watch?v=JctudMfO-7w
(Bilibili) https://www.bilibili.com/video/BV1hf4y1i7bK/


I also attach my EASF webinar to introduce AL-DIC/DVC algorithm and review other DIC/DVC methods:
(Youtube) https://www.youtube.com/watch?v=-t61WrVagZ4
(Bilibili) https://www.bilibili.com/video/BV1ff4y1B71L/



## Citation
* [1] For full details, and to use this code, please cite our paper:
Yang, J. and Bhattacharya, K. Augmented Lagrangian Digital Image Correlation. Exp.Mech. 59: 187, 2018. https://doi.org/10.1007/s11340-018-00457-0. Full text can be requested at: www.researchgate.net/publication/329456141_Augmented_Lagrangian_Digital_Image_Correlation  
* [2] Yang, J. (2019, March 6). 2D_ALDIC (Version 3.3). CaltechDATA. https://data.caltech.edu/records/1443  
% =========================================
* [3] Yang, J. and Bhattacharya, K. Combining Image Compression with Digital Image Correlation. Exp.Mech. 59: 629-642, 2019. https://doi.org/10.1007/s11340-018-00459-y. Full text can be requested at: https://www.researchgate.net/publication/330489954_Combining_Image_Compression_with_Digital_Image_Correlation
* [4] Finite-element-based Global DIC code is also available at:
https://www.mathworks.com/matlabcentral/fileexchange/82873-2d-finite-element-global-digital-image-correlation-fe-dic
* [5] Besides 2D-DIC, our new code "ALDVC" (augmented Lagrangian Digital Volume Correlation) to track deformations in volumetric images is also available:
https://www.mathworks.com/matlabcentral/fileexchange/77019-augmented-lagrangian-digital-volume-correlation-aldvc

## Contact and support
Jin Yang (Caltech solid mechanics, PhD '19): jyang526@wisc.edu  -or-  aldicdvc@gmail.com
I appreciate your comments and ratings to help me further improve this code. If you have other questions, feel free to email me.


##
 
<p align="center">
  <img width="538" height="301" src="https://github.com/jyang526843/2D_ALDIC_v3/blob/master/logo_aldic.png"></p>
  <p align="center">
  <img width="245" height="176" src="https://github.com/jyang526843/2D_ALDIC_v3/blob/master/Example_aldic_foam_compression_strain_eyy.gif"></p>
  <p align="center"><img width="600" height="338" src="https://github.com/jyang526843/2D_ALDIC/blob/master/results_ALDIC_Quadtree_demo/Demo1.gif"></p>
  <p align="center">
 <img width="600" height="338" src="https://github.com/jyang526843/2D_ALDIC/blob/master/results_ALDIC_Quadtree_demo/Demo2.gif"></p>
  <p align="center">
 <img width="600" height="338" src="https://github.com/jyang526843/2D_ALDIC/blob/master/results_ALDIC_Quadtree_demo/Demo3.gif">
</p>


 

