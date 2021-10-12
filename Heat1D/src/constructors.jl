"""
`FDHeatProblem` - Initalize the data structure
"""
function FDHeatProblem(a, b, nx, Δt, tmax, u0, bcl::PeriodicBC, bcr::PeriodicBC, ts::BEuler)
    Δx = (b-a)/(nx+1);
    x = LinRange(a, b, nx+2)[1:end-1];
    nΔt = Int(tmax/Δt)
    A = spzeros(nx+1, nx+1);
    rhs = zeros(nx+1);

    return FDHeatProblem(a, b, Δx, nx, x, Δt, nΔt, tmax, u0, A, rhs, bcl, bcr, ts)
end
