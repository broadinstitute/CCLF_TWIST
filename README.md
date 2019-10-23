# CCLF_TWIST
A pipeline to compute the CN, SNV, ... on TWIST sequences from the CCLF program 

## Usage

0. `git clone https://github.com/broadinstitute/CCLF_TWIST`

1. (Recommended in a [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) or virtual environment) `conda create -n cclf python=3.6`

2. Run `cd CCLF_TWIST && pip install -r requirements.txt` to load all the dependencies

3. information on the pipeline running on Terra is given [here](https://cclf.gitbook.io/tsca/)

4. you also have to get the JKBio package by running (from within the same parent folder as where CCLF_TWIST is): `cd .. && git clone https://github.com/jkobject/JKBIO`

5. you might also have to load the dependencies here: `cd JKBio && pip install -r requirements.txt`


@BroadInstitute

Jérémie Kalfon @jkobject
Gwen Miller
Aniket Shetty