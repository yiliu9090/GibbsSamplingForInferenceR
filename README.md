# Gibbs Sampling For Inference R
 
## Summary 

This is a software for inferences using Bayesian non-parametric on the problem on Transcriptional Bursting. 

This software is exists as a docker container. 

## Installation 

1. Install docker using the following instructions https://docs.docker.com/get-docker/
2. Run `git clone https://github.com/yiliu9090/GibbsSamplingForInferenceR.git`
3. Run `cd GibbsSamplingForInferenceR`
4. Run `docker build . -t [tags]`
5. After it is built, congrats! 

## Running the algorithm 

To run this algorithm, one needs to have a `.json` file and `.txt` file. All the 

    { 
    "NAME":["Example0_1AND1"],
    "DATA_LOCATION":[
        "Data/ExampleData/Example_0_1AND1.txt"
    ],
    "MAXN":50,
    "GAMMAPRIORA": 0.1,
    "GAMMAPRIORB": 1, 
    "ALPHA": [0.0001,0.0003,0.0005,0.0007,
            0.001,0.005,0.007,
            0.01,0.05],
    "MCMCITER": 1000,
    "BURNIN":200,
    "DUMP_LOCATION":[
        "Data/Output/"],
    "NSD":1.00,
    "SEPARATION_FACTOR":[1.00],
    "SEED":2021
    }



This is the data file. 
    1.648249
    11.86002
    0.278794
    2.419705
    7.066532
    8.294041
    ...

