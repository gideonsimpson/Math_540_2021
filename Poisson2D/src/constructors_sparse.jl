"""
`FDPoisson2DProblem` - Initalize the data structure
"""

function FDPoisson2DProblem(a, b, c, d, n, f)

    zero_function = (x,y)-> 0;
    zero_bc = DirichletBC(zero_function);
    return FDPoisson2DProblem(a, b, c, d, n, n, f, zero_bc, zero_bc, zero_bc, zero_bc);
end

function FDPoisson2DProblem(a, b, c, d, nx, ny, f, bct::TBC, bcb::TBC, bcl::TBC, bcr::TBC) where {TBC<:DirichletBC} 
    Δx = (b-a)/(nx+1)
    Δy = (d-c)/(ny+1)
    x = LinRange(a, b, nx+2)[2:end-1];
    y = LinRange(c, d, ny+2)[2:end-1];

    # allocate the matrix and the rhs
    A = spzeros(nx * ny, nx * ny);
    rhs = zeros(nx*ny);


    return FDPoisson2DProblem(a, b, c, d, Δx, nx, Δy, ny, x, y, f, A, rhs, bct, bcb, bcl, bcr)
end
