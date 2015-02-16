## Parallel Computing Toolbox (PCT)

When working with numerical problems long enough, you tend to run into at least one of two ubiquitous problems:  something takes too long to run, or it is too large to process easily.  In either case, _parallelism_ can be used to split the program into more manageable pieces which run side-by-side to achieve a solution.

We will discuss the [Parallel Computing Toolbox](http://www.mathworks.com/products/parallel-computing) today.  Primarily, this is the way that you can exploit multicore and GPU parallelism on a single workstation.

### Parallel Model

MATLAB parallelizes tasks using constructs like `parfor`, a parallel `for` loop.  These constructs variously split control or data between several processes (workers) and then reunify the data in a usable format after the calculation has been carried out.

You can also break up tasks by data, using the idea of SPMD—"Single Program, Multiple Data".  In this case, data across a domain are distributed between processes which then execute the same algorithm on each piece of the array.

To start a group of workers for use with MATLAB, type
    parpool

#### Parallel loops • `parfor`, `parpool`
`parfor` is the workhorse of PCT.  Naïvely, you can expect to just replace your `for`s with `parfor`s and gain the benefit of multiple cores.  In practice, you often have to tweak your algorithms a bit to remove loop dependencies.

As an example, we may want to calculate the result of some complicated nonvector function on a vector.  We'll use `sin` as a stand-in:
    for i = 1:1024
      A(i) = sin(i*2*pi/1024);
    end
    plot(A)

This could be solved in parallel because every loop cycle is independent of all others:
    parfor i = 1:1024
      A(i) = sin(i*2*pi/1024);
    end
    plot(A)

(The first time we run this, it may be rather slow, since the "parallel pool" of worker threads has to be initialized (with `parpool`).  Subsequent runs will be lightning fast, but `parfor` will only be faster than `for` in this case if the function being used is computationally intensive.)

Requirements for `parfor` loops:
-   Task independent
-   Order independent
-   Cannot introduce new variables (_e.g._ `eval`, `load`, `global`, etc.)
-   Cannot change the control flow:  no `break` or `return`
-   Cannot nest `parfor` loops

It's worth noting at this point that you can't exploit parallelism with the built-in ODE solvers like `ode45` to speed up a single DE solution.  These iteratively solve forward from an initial condition and thus each loop is not independent of the others.  However, we could (for instance) perform a parameter sweep search by solving with `ode45` in parallel for a set of parameters:

<blockquote>
A _parameter sweep_ consists of an equation (in this case, a differential equation) and a set of parameters which will be varied in order to search for some salient features of the system.
</blockquote>

**Problem**.  We wish to search through the parameter space of the equation
$$\ddot{x} + a \dot{x} + b x = 0$$
in $a$ and $b$ in order to find the maximum resulting values of the differential equation's solution for each parameter pair.

**Solution**.  The second-order linear ODE
$$\ddot{x} + a \dot{x} + b x = 0$$
converts to the system of coupled first-order ODEs
$$x' = y$$
$$y' = -a y - b x$$

The function required converts this problem into a system of equations:
    %%  Parameter Sweep ODE
    function [ X ] = dx( t, X0, a, b )
        T = [-a -b ;
              1  0];
        X = T * X0;
    end

A representative solution can be found as follows:
    % Initial condition
    X0 = [1 1]';
    [t x] = ode45(@(t, x) dx(t, x, 0.5, 0.25), [0 10], X0);
    plot(t,x)

**Serial**.  In this case, we will search for the maximum value of the output and store that as the key variable in a solution array.  A serial code is included as `sweep.m`:
    %%  Parameter Sweep Search Script
    % Parameters to test
    A = linspace(-1,1,150);
    B = linspace(0.5,2,100);

    % Solution array
    X_max = zeros(length(A), length(B));

    % Initial condition
    X0 = [1 1]';

    for i = 1:length(A)
        for j = 1:length(B)
            [t x] = ode45(@(t, x) dx(t, x, A(i), B(j)), [0 10], X0);
            X_max(i,j) = max(x(:,1));
        end
    end

To examine the output, take a look at the surface plot:
    surf(B,A,X_max)

**Parallel**.  Since the loops of this code are independent, it can be trivially parallelized simply by dropping a `parfor` in for the inner `for`.

Now compare the relative timings of the two codes for a large parameter space:
    A = linspace(-1,1,150);
    B = linspace(0.5,2,100);
    
    tic; sweep; toc
    tic; sweep_parallel; toc
    surf(B,A,X_max)

Conveniently, if there is no parallel pool available, `parfor` simply reduces to `for`-like serial behavior.

If `parfor` isn't adequate to your problem, there are several other options, like batch scripting and distributed array computing.


#### Batch scripts • `batch`


#### Distributed data & control • `distributed`, `spmd`

(also some good stuff here:  http://www.bu.edu/tech/support/research/training-consulting/online-tutorials/matlab-pct/distributed-array-examples/)

**Data**.  What we have seen thus far is best described by the notion of "control-level" parallelism (although a simple subset of this idea).  But we can also split a program over data sets.  The MATLAB construct for this is the `distributed` matrix.

How does a `distributed` array differ from a regular one?  Consider this:
    spmd
        % Each worker creates its own 80x100 array.
        A = zeros(80, 1000);
        % One 80x100 array is created, and it is partitioned between workers.
        D = codistributed(A)
    end
    whos

(Note that you can _see_ the data which were created, but you can't access them outside of the `spmd` loop.  Also, `codistributed` arrays inside of `spmd` statements are `distributed` arrays outside of them.)
    D = gather(D);
    
    n = 10;
    F = distributed(magic(n)); % Distribute array to workers
    spmd
        Q = F .* 2;
    end
    M = gather(F)              % Return array to client

**Control**.  In an `spmd` block, each worker has access to the unique index `labindex` and the total number of workers 'numlabs'.  We can use this to decompose the problem.

For instance, let's integrate
$$ \frac{x + \sqrt{1+x}}{\sqrt[3]{1+x}} $$
over the range [-1,1] across the number of workers we have on our system.
    % Integrand:
    f = @(x) (x+sqrt(1+x))/((1+x).^(1/3));
    
    % Enter an SPMD block to run the enclosed code in parallel on a number
    % of MATLAB workers:
    spmd
        lhs = (labindex - 1) / numlabs * 2 - 1;
        rhs = labindex / numlabs * 2 - 1;
        
        this_integral = integral(f, lhs, rhs, 'ArrayValued', true);
        total_integral = gplus(this_integral, 1);
    end

(Again, the examples we are forced to use in a workshop are not computationally intensive enough to really show the power of parallel computing.)

Another example ([source](https://stackoverflow.com/questions/13071485/matlab-slow-parallel-processing-with-distributed-arrays)):
    clear all;
    
    D = rand(512, 512, 3);
    S = size(D);
    [fx, fy, fz] = gradient(D);
    
    % this part could also be parallelized - at least a bit.
    tic;
    DHess = zeros([3 3 S(1) S(2) S(3)]);
    [DHess(1,1,:,:,:), DHess(1,2,:,:,:), DHess(1,3,:,:,:)] = gradient(fx);
    [DHess(2,1,:,:,:), DHess(2,2,:,:,:), DHess(2,3,:,:,:)] = gradient(fy);
    [DHess(3,1,:,:,:), DHess(3,2,:,:,:), DHess(3,3,:,:,:)] = gradient(fz);
    toc
    
    % a sequential implementation
    d = zeros([3, S(1) S(2) S(3)]);
    disp('sequential')
    tic
    for i = 1 : S(1)
        for ii = 1 : S(2)
            for iii = 1 : S(3)
                d(:,i,ii,iii) = eig(squeeze(DHess(:,:,i,ii,iii)));
            end
        end
    end
    toc
    
    % a parallel implementation
    disp('parallel')
    tic
    spmd
        % just for information
        disp(['lab ' num2str(labindex)]);
        
        % distribute the input data along the third dimension
        % This is the dimension of the outer-most loop, hence this is where we
        % want to parallelize!
        DHess_dist  = codistributed(DHess, codistributor1d(3));
        DHess_local = getLocalPart(DHess_dist);
        
        % create an output data distribution - 
        % note that this time we split along the second dimension
        codist = codistributor1d(2, codistributor1d.unsetPartition, [3, S(1) S(2) S(3)]);
        localSize = [3 codist.Partition(labindex) S(2) S(3)];
        
        % allocate local part of the output array d
        d_local = zeros(localSize);
        
        % your ordinary loop, BUT! the outermost loop is split amongst the
        % threads explicitly, using local indexing. In the loop only local parts
        % of matrix d and DHess are accessed
        for i = 1:size(d_local,2)
            for ii = 1 : S(2)
                for iii = 1 : S(3)
                    d_local(:,i,ii,iii) = eig(squeeze(DHess_local(:,:,i,ii,iii)));
                end
            end
        end
        
        % assemble local results to a codistributed matrix
        d_dist = codistributed.build(d_local, codist);
    end
    toc
    
    isequal(d, d_dist)


#### GPU programming • `gpuArray`, `arrayfun`, `feval`


#### MapReduce programming • `mapreduce`


#### Profiling • `mpiprofile`



## Distributed Compute Server (DCS)

DCS is used on clusters to allow large-scale solution of MATLAB problems.

## Tips

    pools = matlabpool('size');
    cpus = feature('numCores'); %undocumented feature
    if pools ~= (cpus - 1)
        if pools > 0
            matlabpool('close');
        end
        matlabpool('open', cpus - 1);
    end
