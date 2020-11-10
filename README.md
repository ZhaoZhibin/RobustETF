# RobustETF

# Robust Enhanced Trend Filtering with Unknown Noise


This repository contains the implementation details of our paper: [Signal Processing]
"[**Robust Enhanced Trend Filtering with Unknown Noise**](https://www.sciencedirect.com/science/article/pii/S0165168420304333)" 
by [Zhibin Zhao](https://zhaozhibin.github.io/). 


## About
One important step in time series analysis is the extraction of an underlying trend. However, the true trend is often submerged by complex background noise, espcially non-Gaussian noise or outliers. Accurate trend extraction against outliers from a raw signal is a challenging task. To address this challenge, this paper extends $l_1$ trend filtering to a robust enhanced trend filtering called RobustETF by combining mix of Gaussian (MoG) and non-convex sparsity-inducing functions. We first model the noise as a MoG distribution to allow RobustETF to be robust in the presence of any type of non-Gaussian noise or outliers. After that, to handle the biased estimation of the $l_1$ norm, we use the Gibbs distribution embedding smoothed and non-convex sparsity-inducing functions to faithfully preserve the amplitude of the trend. Furthermore, we design an extended EM algorithm to solve the resulting non-convex optimization problem. Finally, we show the results of experiments on both real-world and synthetic data to compare the performance of the proposed algorithm against other state-of-the-art methods.

## Dependencies
- Matlab R2016b

## Pakages

This repository is organized as:
- [funs](https://github.com/ZhaoZhibin/RobustETF/tree/master/funs) contains the main functions of the algorithm.
- [util](https://github.com/ZhaoZhibin/RobustETF/tree/master/util) contains the extra functions of the test.
- [Results](https://github.com/ZhaoZhibin/RobustETF/tree/master/Results) contains the results of the algorithm.
- [Data](https://github.com/ZhaoZhibin/RobustETF/tree/master/Data) contains the data generated by NSF I/UCR Center on Intelligent Maintenance Systems (IMS) with support from Rexnord Corp. in Milwaukee, WI.
- [tqwt_matlab_toolbox](https://github.com/ZhaoZhibin/RobustETF/tree/master/tqwt_matlab_toolbox) contains the TQWT toolbox copied from I. W. Selesnick.
In our implementation, **Matlab R2016b** is used to perform all the experiments.

## Implementation:
Flow the steps presented below:
-  Clone this repository.
```
git clone https://github.com/ZhaoZhibin/RobustETF.git
open it with matlab
```
-  Test Simulation: Check the parameters setting of simulation in `Config.m` and run `Test_simulaton.m`. 
-  Test Your Own Data: Run `Test_Your_Own_Data.m`. 


## Citation
If you feel our RobustETF is useful for your research, please consider citing our paper: 

```
@article{zhao2020Robust,
  title={Robust Enhanced Trend Filtering with Unknown Noise},
  author={Zhao, Zhibin and Wang, Shibin and Wong, David and Sun, Chuang and Ruqiang Yan and Chen, Xuefeng},
  journal={Signal Processing},
  year={2020},
  publisher={Elsevier}
}
```
## Contact
- zhibinzhao1993@gmail.com


