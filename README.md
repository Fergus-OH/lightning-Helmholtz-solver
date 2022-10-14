# Lightning Helmholtz Solver

<img src="https://user-images.githubusercontent.com/77754538/194128251-c13a9ae0-7c2e-4f30-80af-4c061708221e.gif" width="100%" height="100%">


[https://user-images.githubusercontent.com/77754538/194108726-112bab3b-30c4-443d-9f5c-154a378ba39d.mp4]::


## Usage

`U = helmholtz(wavenum, P, g)` solves the Helmholtz equation with
Dirichlet boundary data on the simply-connected region $\Omega$ bounded by `P`, which may be a polygon or circular polygon.




Repeating positional arguments


      
### Input arguments

#### Required positional arguments


- `wavenum`

    (integer) The attributes that every function declared with this
    keyword should have (in the form of source code, with a leading `:`).
- `P`
    (integer, vector, or string) 
    vector of corners as complex numbers $z = x+iy$ in counterclockwise order to specify a polygon or cell array of corners `v` and pairs `[v r]` to specify a circular polygon: $r =$ radius of curvature of arc from this v to the next or one of the following specified strings `'sqr'`[square], `'rec'`[tangle], `'snow'`[flake], `'pent'`[agaon], `'hex'`[agon], `'L'`, `'circleL'`, or `'C'` or integer $\ge 3$, the number of corners of a random polygon or integer $\le -3, -1$ $\times$ no. of corners of a random circular polygon.



#### Optional positional arguments
- `g` function handle for Dirichlet boundary data that satisfies helm(g) or cell array of function handles for sides P1-P2, P2-P3, (default `@(z) exp(-1i*real(wavenum*exp(-1i*z0ang)*z)))` for `wavenumber`$>0$ `@(z) @(z) besselh(0,-wavenum*abs(z-(z0_pt)))` for `wavenumber`$<0$)

#### Optional name-value arguments

Name-value arguments
Optional pairs of arguments as Name1=Value1,...,NameN=ValueN, where Name is the argument name and Value is the corresponding value. Name-value arguments must appear after other arguments, but the order of the pairs does not matter.


- `tol`       float - 1e-6 - tolerance etc
- `z0`        complex number
- `fs`        float+ -  12? - set font size for plots

<!-- | Parameter   | Type | Default | Description |
| :---------- | :--: | :------:| :-----------|
| `tol`       | float | 1e-6 | tolerance |
| `z0`        | complex number |
| `noplots`   | flag |
| `noplots3d` | flag |
| `steps`     | flag |
| `scat`      | flag |
| `slow`      | flag |
| `fs`        | float+ | 12? | set font size for plots | -->



#### flags

the following flag parameters can be specified



| Flag        | Description |
| :---------- | :-----------|
| `noplots`   | surpresses plotting |
| `noplots3d` | surpresses 3D surface plotting |
| `steps`     | for step-by-step plots of errors on boundary and poles |
| `scat`      | to plot only the scattered field |
| `slow`      | to turn off adaptive mode for cleaner root-exponential convergence curves |








If you don't specify a particular option, its default value is used. The
available configuration options are:

- `g`

    (string) The attributes that every function declared with this
    keyword should have (in the form of source code, with a leading `:`).

    Default: nothing

- `tol`

    (string) The attributes that every function declared with this
    keyword should have (in the form of source code, with a leading `:`).

    Default: nothing
- `tol`
- `z0`
- `noplots`
- `noplots3d`
- `steps`
- `scat`
- `slow`
- `fs`





    
### Examples

Examples:

`helmholtz(50,'sqr');`                              % plane with square

`helmholtz(-50,'sqr');`                             % point with square

`helmholtz(-20,'pent','tol',1e-10);`                % pentagon

`helmholtz(-20,'circleL','z0',2+3i);`               % circular L-shape

`helmholtz(20,'bullet','z0',1);`                    % bullet

`helmholtz(-30,'snow','steps')`                     % snowflake

`helmholtz(50,[1/2*exp(2i*pi*([1:3])/3)],'z0',1i)`  % triangle

% two point sources:
```matlab
wavenum = -30; z0_pt = .5+1i;
g = @(z) besselh(0,-wavenum*abs(z-(z0_pt))) + besselh(0,-wavenum*abs(z-(-z0_pt')));
helmholtz(wavenum,'sqr',g,'noplot3d');
```
