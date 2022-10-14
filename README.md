# Lightning Helmholtz Solver

<img src="https://user-images.githubusercontent.com/77754538/194128251-c13a9ae0-7c2e-4f30-80af-4c061708221e.gif" width="100%" height="100%">


[https://user-images.githubusercontent.com/77754538/194108726-112bab3b-30c4-443d-9f5c-154a378ba39d.mp4]::


## Usage

U = HELMHOLTZ(wavenum,P,G) solves the Helmholtz equation with
Dirichlet boundary data on the simply-connected region Omega bounded by P, which may be a polygon or circular polygon.
      

### Required
- `wavenum`

    (integer) The attributes that every function declared with this
    keyword should have (in the form of source code, with a leading `:`).
- `P`
    (integer, vector, or string) vector of corners as complex numbers $z = x+iy$ in counterclockwise order to specify a polygon or cell array of corner $v$ and pairs $[v r]$ to specify a circular polgon: $r = $ radius of curvature of arc from this $v$ to the next or one of the following specifird strings 'sqr'[square], 'rec'[tangle], 'snow'[flake], ...
pent[agaon],
            'hex'[agon], 'L', 'circleL', or 'C'
or integer ≥ 3, the number of corners of a random polygon
or integer ≤ -3, -1 x no. of corners of a random ... circular polygon].

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
