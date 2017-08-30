# rpws

Gridded evaluation of
samples from the ensemble of Gaussian random plane waves in
two and three dimensions. MATLAB code driving a NUFFT library.

### Alex H. Barnett, (c) 2006-2017

## Dependencies

- (MATLAB)[http://mathworks.com] or octave
- (FINUFFT)[http://github.com/ahbarnett/finufft]

## Description

This code evaluates on a regular square or cubical grid of points
a sample drawn at random from the set of eigenfunctions of the
Laplacian with eigenvalue 1.

the function
```math
u(x) = \mbox{Re } \lim_{n\to\infty} \sum_{j=1}^n a_n e^{in_j \cdot x}
```

The mean of $u$ is 0 and the variance is 1.

See https://math.dartmouth.edu/~ahb/rpws/

