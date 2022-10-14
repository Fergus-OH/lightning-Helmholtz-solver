# Lightning Helmholtz Solver

<img src="https://user-images.githubusercontent.com/77754538/194128251-c13a9ae0-7c2e-4f30-80af-4c061708221e.gif" width="100%" height="100%">


[https://user-images.githubusercontent.com/77754538/194108726-112bab3b-30c4-443d-9f5c-154a378ba39d.mp4]::


## Usage

### Required
- `wavenum`

    (integer) The attributes that every function declared with this
    keyword should have (in the form of source code, with a leading `:`).
- `P`
    (integer, vector, or string) vector of corners as complex numbers z = x+iy.
     $$\left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right)$$

### Optional

If you don't specify a particular option, its default value is used. The
available configuration options are:

- `attributes`

    (string) The attributes that every function declared with this
    keyword should have (in the form of source code, with a leading `:`).

    Default: nothing

- `attributes`

    (string) The attributes that every function declared with this
    keyword should have (in the form of source code, with a leading `:`).

    Default: nothing
    
### Examples
