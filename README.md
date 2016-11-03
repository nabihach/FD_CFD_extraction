# FD_CFD_extraction

This repository contains the implementation of two algorithms, TANE and CTANE, corresponding to the following publications:

1. "TANE: An Efficient Algorithm for Discovering Functional and Approximate Dependencies" (link: https://www.lri.fr/~pierres/donn%E9es/save/these/articles/lpr-queue/huhtala99tane.pdf)

2. "Discovering Conditional Functional Dependencies" (link: http://homepages.inf.ed.ac.uk/fgeerts/pdf/CFDdiscovery.pdf)

We have also provided several CSV files as test data.

This code was used in the following work:
"Automatic Discovery of Functional Dependencies and Conditional Functional
Dependencies: A Comparative Study" (link: https://cs.uwaterloo.ca/~nasghar/848.pdf)

##Running the code

To run tane.py on a particular csv file (e.g. adult.csv), execute the following command in your terminal:
```
python tane.py adult.csv
```
To run ctane.py on the same data, execute:
```
python ctane.py adult.csv
```
To run ctane.py and obtain k-frequent CFDs, execute:
```
python ctane.py adult.csv k
```
where k is your integer of choice.

