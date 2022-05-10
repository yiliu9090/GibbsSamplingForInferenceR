# Gibbs Sampling For Inference R
 
## Summary 

This is a software for inferences using Bayesian non-parametric on the problem on Transcriptional Bursting. 

This software is exists as a docker container. 

## Installation 

1. Install docker using the following instructions https://docs.docker.com/get-docker/
2. Run `git clone https://github.com/yiliu9090/GibbsSamplingForInferenceR.git`
3. Run `cd GibbsSamplingForInferenceR`
4. Run `docker build . -t [name of container]`
5. After it is built, congrats! 

## Running the algorithm 

### Inputs 

To run this algorithm, one needs to have a `.json` file and `.txt` file. All the configurations are in `.json` file. 
In this system, we have an example `DataExample.json` file. 
We will illustrate how to configure this using the following example. 

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

We will define each term here. 
- `"NAME"`: list of Names, multiple names 
- `"DATA_LOCATION"`: list of data location of data with respect to the container 
- `"MAXN"`: Maximumn N_0
- `"GAMMAPRIORA"`: a from prior Gamma(a,b)
- `"GAMMAPRIORB"`: b from prior Gamma(a,b)
- `"ALPHA"`: list of alpha to scan through 
- `"MCMCITER"`: Number of MCMC iteration to run through
- `"BURNIN"`: Number of samples thrown away as burn in 
- `"DUMP_LOCATION"`: list of folders to put all the output data 
- `"NSD"`: Number of standard deviation before overfitting
- `"SEPARATION_FACTOR"`: the minimum value for $\frac{\lambda_i}{\lambda_j}$ if $\lambda_i > \lambda_j$
- `"SEED"`: Random Seed



This is the data file with the timings. `Data/ExampleData/Example_0_1AND1.txt` is an example

    1.648249
    11.86002
    0.278794
    2.419705
    7.066532
    8.294041
    ...

### Running the container 
With this two files ready, one can then run the code using the following docker command 

    docker container run -d --rm -v [data location in your pc]:[data location in the container] [name of the container] [json file location in the container]

For example, lets call our container `gitgibbstrial` (in step 4 of Installation, run `docker build . -t gitgibbstrial`)

Run `docker container run -d --rm -v /[your file location]/GibbsSamplingForInferenceR/Data:/Workspace/Data gitgibbstrial Data/JSON/DataExample.json`

If one is not familar with docker commands, what one can do is to just use this line to run and change the file content accordingly.

### Outputs 

- Plot of likelihood against iteration number for each $\alpha$
- Plot of the histogram of the posterior $\lambda$ for each $\alpha$
- Plot of the $\lambda$ estimates against $N$
- Plot of $\lambda$ estimates against $\alpha$ and $\log(\alpha)$

One can find an example outputs in Data/Output