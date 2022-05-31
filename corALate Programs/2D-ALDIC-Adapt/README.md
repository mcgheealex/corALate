# 2D_ALDIC 
AL-DIC(Augmented Lagrangian DIC) is a fast, parallel-computing hybrid DIC algorithm, which combines advantages of local subset DIC method (fast computation speed, and parallel computing) and finite-element-based global DIC method (guarantee global kinematic compatibility and decrease noise).  

## Advantages of AL-DIC algorithm
* [1] It’s a fast algorithm using distributed parallel computing.  
* [2] Global kinematic compatibility is added as a global constraint in the form of augmented Lagrangian, and solved using Alternating Direction Method of Multipliers scheme.
* [3] Both displacement fields and affine deformation gradients are correlated at the same time.
* [4] No need of much manual experience about choosing displacement smoothing filters.
* [5] It works well with compressed DIC images and adaptive mesh. See our paper: Yang, J. & Bhattacharya, K. Exp Mech (2019). https://doi.org/10.1007/s11340-018-00459-y;
* [6] Both accumulative and incremental DIC modes are implemented to deal with image sequences, which is especially quite useful for very large deformations.

## Prerequisites & Installation
AL-DIC MATLAB code was tested on MATLAB versions later than R2018a. Both single thread and parallel computing features are included in AL-DIC code. Please download and unzip the code to the MATLAB working path. Then, execute the mail file main_ALDIC.m.

## Citation
* [1] Yang, J. (2019, March 6). 2D_ALDIC (Version 3.3). CaltechDATA. https://data.caltech.edu/records/1443
* [2] Yang, J. and Bhattacharya, K. Augmented Lagrangian Digital Image Correlation. Exp.Mech. 59: 187, 2018. https://doi.org/10.1007/s11340-018-00457-0.   or 
* Full text can be requested at: www.researchgate.net/publication/329456141_Augmented_Lagrangian_Digital_Image_Correlation  


## Contact and support
Jin Yang (Caltech solid mechanics, PhD '19): jyang526@wisc.edu  -or-  aldicdvc@gmail.com

Welcome to give the ALDIC code ratings and comments in the MATLAB File Exchange community: [![View Augmented Lagrangian Digital Image Correlation and Tracking on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/70499-augmented-lagrangian-digital-image-correlation-and-tracking)

##
 
<p align="center">
  <img width="538" height="301" src="https://github.com/jyang526843/2D_ALDIC_v3/blob/master/logo_aldic.png">
  <img width="245" height="176" src="https://github.com/jyang526843/2D_ALDIC_v3/blob/master/Example_aldic_foam_compression_strain_eyy.gif">
</p>


 

