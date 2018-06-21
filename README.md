This repo holds all code that is used to simulate and replicate *all* results
reported in our arXiv paper: https://arxiv.org/abs/1711.05833

The code is individual Julia scripts with self-evident names and a description
comment at the head of each file.

The code was written in `julia 0.6`, `OrdinaryDiffEq 3.15.1`
and `DynamicalBilliards 1.6.3`. However, I am very confident that almost everything
runs on `julia 0.7` and `DynamicalBilliards 2.0.0`. The latest versions
give massive performance benefits, especially for the billiards, so it is advisable
to use those instead.

For any questions, contact me at my institute mail address (see arXiv paper).
