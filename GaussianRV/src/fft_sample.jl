function sample_field(N)
    uhat = zeros(ComplexF64,2*N); # preallocate space

    ξ = randn(N); # generate the random variables

    # construct the eigenvalues
    k = 1:N;
    λ = @. 1/(π*k)^2;

    # fill in the nonzero entries
    # NOTE we need to multiply by 2 *N for FFT scaling
    @. uhat[2:N+1] = 2 * N * sqrt(λ) * sqrt(2) * ξ;

    # invert and get the relevant imaginary part
    u = imag.(ifft(uhat))[N+2:end];
    return u
end