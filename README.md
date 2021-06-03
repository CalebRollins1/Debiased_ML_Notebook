# Debiased_ML_Notebook

This guide will walk through R code for the debiased machine learning algorithm. The example used looks at the impact of access to a 401k on lifetime savings. This guide has two goals: 

1) Help the reader understand the technical details of the debiased machine learning algorithm
2) Provide a framework the reader can use to implement debiased machine learning with their own data and specifications

In the PCR folder, there are the following csv files with results:

coverage_exp100: Simulation with 100 iterations comparing biased vs. unbiased using measurement error

dropout100: Simulation with 500 iterations where X is 100 x 100 using missingness

dropout500: Simulation with 500 iterations where X is 500 x 500 using missingness

final_exp100: Simulation with 500 iterations where X is 100 x 100 using measurement error

final_exp500: Simulation with 500 iterations where X is 500 x 500 using measurement error

PCR_results: Results running PCR on simulated data with different number of principal components kept



