"""
`assemble_system!` - assemble the linear system for solving
"""
function assemble_system!(problem::FDPoisson2DProblem{TF, TBC, TBC, TBC, TBC}) where {TF<:AbstractArray, TBC<:DirichletBC}

    N = problem.nx * problem.ny; # get system size


    for k in 1:N
        # get local coordiantes from global index
        # k = i + (j-1) * nx
        j, i = fldmod1(k, problem.nx); 
        # process corners, then edges, then interior
        if (i==1 && j==1)
            problem.A[k, 1] = 2/problem.Δx^2 + 2/problem.Δy^2;
            problem.A[k, 2] = -1/problem.Δx^2;
            problem.A[k, problem.nx + 1]  = -1/problem.Δy^2;
        elseif (i==1 && j == problem.ny)
            problem.A[k, (problem.ny-2) * problem.nx + 1] =-1/problem.Δy^2;
            problem.A[k, (problem.ny-1) * problem.nx + 1] = 2/problem.Δx^2 + 2/problem.Δy^2;
            problem.A[k, (problem.ny-1) * problem.nx + 2] = 1/problem.Δx^2
        elseif (i==problem.nx && j==1)
            problem.A[k,problem.nx-1] =  -1/problem.Δx^2;
            problem.A[k,problem.nx] = 2/problem.Δx^2 + 2/problem.Δy^2;
            problem.A[k,problem.nx + 1] =  -1/problem.Δy^2;
        elseif (i==problem.nx && j ==problem.ny)
            problem.A[k, (problem.nx-1) * problem.ny] = -1/problem.Δy^2;
            problem.A[k, problem.nx * problem.ny - 1] = -1/problem.Δx^2;
            problem.A[k, problem.nx * problem.ny] = 2/problem.Δx^2 + 2/problem.Δy^2;
        elseif(i==1 && j>1 && j <problem.ny)
            # left edge
            problem.A[k, 1 + (j-2) * problem.nx] = -1/problem.Δy^2;
            problem.A[k, 1 + (j-1) * problem.nx] = 2/problem.Δx^2 + 2/problem.Δy^2;
            problem.A[k, 2 + (j-1) * problem.nx] = -1/problem.Δx^2;
            problem.A[k, 1 + j * problem.nx] = -1/problem.Δy^2;
        elseif(i==problem.nx  && j>1 && j <problem.ny)
            # right edge
            problem.A[k, (j-1) * problem.nx] = -1/problem.Δy^2;
            problem.A[k, j * problem.nx-1] = -1/problem.Δx^2;;
            problem.A[k, j * problem.nx] = 2/problem.Δx^2 + 2/problem.Δy^2;
            problem.A[k, (j+1) * problem.nx] = -1/problem.Δy^2;
        elseif(i>1 && i<problem.nx && j==1)
            # bottom edge
            problem.A[k, i-1] = -1/problem.Δx^2;
            problem.A[k, i] = 2/problem.Δx^2 + 2/problem.Δy^2;
            problem.A[k, i+1] = -1/problem.Δx^2;
            problem.A[k, i + problem.nx] = -1/problem.Δy^2;
        elseif(i>1 && i<problem.nx && j==problem.ny)
            # top edge
            problem.A[k, i + (problem.ny-2)*problem.nx] = -1/problem.Δy^2;
            problem.A[k, i-1 + (problem.ny-1)*problem.nx] =-1/problem.Δx^2;
            problem.A[k, i + (problem.ny-1)*problem.nx] =  2/problem.Δx^2 + 2/problem.Δy^2;
            problem.A[k, i+1 + (problem.ny-1)*problem.nx] = -1/problem.Δx^2;
        else
            problem.A[k, i + (j-2)*problem.nx] = -1/problem.Δy^2;
            problem.A[k, i-1 + (j-1)*problem.nx] = -1/problem.Δx^2;
            problem.A[k, i + (j-1)*problem.nx] = 2/problem.Δx^2 + 2/problem.Δy^2;;
            problem.A[k, i+1 + (j-1)*problem.nx] = -1/problem.Δx^2;
            problem.A[k, i + j*problem.nx] = -1/problem.Δy^2;
        end
    end
    @. problem.rhs = problem.f[:];
    problem

end

"""
`assemble_system!` - assemble the linear system for solving
"""
function assemble_system!(problem::FDPoisson2DProblem{TF, TBC, TBC, TBC, TBC}) where {TF<:Function, TBC<:DirichletBC}

    N = problem.nx * problem.ny; # get system size


    for k in 1:N
        # get local coordiantes from global index
        # k = i + (j-1) * nx
        j, i = fldmod1(k, problem.nx); 
        # process corners, then edges, then interior
        if (i==1 && j==1)
            problem.A[k, k] = 2/problem.Δx^2 + 2/problem.Δy^2;
            problem.A[k, k+1] = -1/problem.Δx^2;
            problem.A[k, k+problem.nx ]  = -1/problem.Δy^2;
        elseif (i==1 && j == problem.ny)
            problem.A[k, k-problem.nx] = -1/problem.Δy^2;
            problem.A[k, k] = 2/problem.Δx^2 + 2/problem.Δy^2;
            problem.A[k, k+1] = -1/problem.Δx^2
        elseif (i==problem.nx && j==1)
            problem.A[k,k-1] =  -1/problem.Δx^2;
            problem.A[k,k] = 2/problem.Δx^2 + 2/problem.Δy^2;
            problem.A[k,k+problem.nx] =  -1/problem.Δy^2;
        elseif (i==problem.nx && j ==problem.ny)
            problem.A[k, k-problem.nx] = -1/problem.Δy^2;
            problem.A[k, k-1] = -1/problem.Δx^2;
            problem.A[k, k] = 2/problem.Δx^2 + 2/problem.Δy^2;
        elseif(i==1 && j>1 && j <problem.ny)
            # left edge
            problem.A[k, k-problem.nx] = -1/problem.Δy^2;
            problem.A[k, k] = 2/problem.Δx^2 + 2/problem.Δy^2;
            problem.A[k, k+1] = -1/problem.Δx^2;
            problem.A[k, k+problem.nx] = -1/problem.Δy^2;
        elseif(i==problem.nx  && j>1 && j <problem.ny)
            # right edge
            problem.A[k, k-problem.nx] = -1/problem.Δy^2;
            problem.A[k, k-1] = -1/problem.Δx^2;;
            problem.A[k, k] = 2/problem.Δx^2 + 2/problem.Δy^2;
            problem.A[k, k + problem.nx] = -1/problem.Δy^2;
        elseif(i>1 && i<problem.nx && j==1)
            # bottom edge
            problem.A[k, k-1] = -1/problem.Δx^2;
            problem.A[k, k] = 2/problem.Δx^2 + 2/problem.Δy^2;
            problem.A[k, k+1] = -1/problem.Δx^2;
            problem.A[k, k + problem.nx] = -1/problem.Δy^2;
        elseif(i>1 && i<problem.nx && j==problem.ny)
            # top edge
            problem.A[k, k-problem.nx] = -1/problem.Δy^2;
            problem.A[k, k-1] =-1/problem.Δx^2;
            problem.A[k, k] =  2/problem.Δx^2 + 2/problem.Δy^2;
            problem.A[k, k+1] = -1/problem.Δx^2;
        else
            problem.A[k, k-problem.nx] = -1/problem.Δy^2;
            problem.A[k, k-1] = -1/problem.Δx^2;
            problem.A[k, k] = 2/problem.Δx^2 + 2/problem.Δy^2;;
            problem.A[k, k+1] = -1/problem.Δx^2;
            problem.A[k, k+problem.nx] = -1/problem.Δy^2;
        end
        problem.rhs[k] = problem.f(problem.x[i], problem.y[j]);
    end
    problem

end
