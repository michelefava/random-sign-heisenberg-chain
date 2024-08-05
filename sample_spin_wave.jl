using LinearAlgebra
using Random
using DelimitedFiles
using Distributions

# given a NxN matrix A, returns a N-element vector.
# Its k-th component is the average of the entries 
# of A on the two diagonals shifted by k and -k
# w.r.t. the main diagonal
function diagonal_means(A)
    N = size(A,1)
    sums = zeros(N)
    for offset in 0:N-1
        v1 = [A[j, j+offset] for j in 1:N-offset]
        v2 = [A[j+offset, j] for j in 1:N-offset]
        sums[offset+1] = (sum(v1) + sum(v2)) / (2 * (N - offset))
    end
    return sums
end

# return a random vector of signs cs
# s.t. sum_j c_j ≠ 0 and not all c_j are equal
function random_cj_vector(L, random_number_generator)
    cs = zeros(L) 
    while true
        cs = rand(random_number_generator, Bernoulli(0.5), L)
        cs = 2cs .- 1
        #exlude special cases where the regularization might be ill-defined
        if sum(cs)!=0 && abs(sum(cs))<L
            break
        end
    end
    return cs
end

# return sites j where c_j=+1
function A_sites(cs)
    return findall(x->x==1, cs)
end

# return sites j where c_j=-1
function B_sites(cs)
    return findall(x->x==-1, cs)
end

# given the disorder realization cs,
# return the square matrix D defined in Eq. (S60)
function quadratic_bosonic_form(cs, uniform_field=0, 
            boundary_field=0, periodic_boundary_conditions=False)

    PBC = periodic_boundary_conditions #shorthand

    boundary_field = abs(boundary_field)
    uniform_field = abs(uniform_field)
    if PBC
        @assert boundary_field==0
    end

    Asites = A_sites(cs)
    Bsites = B_sites(cs)

    NA = length(Asites)
    NB = length(Bsites)

    L = length(cs)

    #construct the matrix D defined in Eq. (D3) of the supplementary material
    D = zeros(L, L)

    shuffled_sites = [Asites; Bsites]

    for j_ind in eachindex(shuffled_sites)

        j = shuffled_sites[j_ind]

        D[j_ind, j_ind] = 2

        if !PBC && j==1
            D[j_ind, j_ind] = 1 + boundary_field
        end

        if !PBC && j==L
            D[j_ind, j_ind] = 1
            continue
        end

        next_j = j % L + 1

        i_ind = findfirst(i->i==next_j, shuffled_sites)

        if cs[j] == cs[next_j]
            D[j_ind, i_ind] = -1
            D[i_ind, j_ind] = -1
        else
            D[j_ind, i_ind] = 1
            D[i_ind, j_ind] = 1 
        end

    end

    # add a uniform field
    sigma = (NA>NB ? +1 : -1) #set the direction so that it stabilizes the spin-glass order
    D +=  uniform_field * sigma * diagm([ones(NA); zeros(NB)])
    D += -uniform_field * sigma * diagm([zeros(NA); ones(NB)])
    
    return D
end

# given a quadratic form D, and a matrix J0 as in Eq. (S64),
# return the vector Delta and the matrix T^(-1) in Eq. (S66)
# if J0_hat is specified, then J0 is the sign matrix
# expressed in the alpha basis, whereas J0_hat is
# the sign matrix in the beta basis.
# Eq. (S66) of the paper then becomes T J0 T^dag = J0_hat
function diagonalize_quadratic_form(D, J0, J0_hat=J0)

    C = cholesky(D)
    #Note: C.L * C.U ≈ D

    # Note: in terms of the notation introduced above Eq. (S65),
    # K=C.U

    M = C.U * J0 * C.L
    @assert norm(M-M')<1e-10
        
    F = eigen(Symmetric(M)) #see Eq. (S66)

    # reorder so that the eigenvalues are in decreasing order
    Lambda = reverse(F.values)
    O = reverse(F.vectors, dims=2)

    Delta = J0_hat * Lambda

    # ensure that negative elements of Delta are compatible with numerical errors
    @assert minimum(Delta)>-1e-12
    # and set them to zero
    Delta[Delta.<0] .= 0

    T = diagm(1 ./ sqrt.(Delta)) * O' * C.U # Eq. (S66)
    #@show norm(T' * J0 * T - J0) #should be very small

    Tinv = inv(T)
    #@show norm(Tinv' * D * Tinv - diagm(E)) #should be very small

    return Delta, Tinv
end

function entropy_from_n(n_mode, n_mode_cutoff)
    # use n_mode to compute the entropy
    n_mode[ n_mode .< n_mode_cutoff ] .= n_mode_cutoff
    exp_minus_e_mode = n_mode ./ (n_mode .+ 1)
    exp_minus_e_mode[ exp_minus_e_mode .< n_mode_cutoff ] .= 0
    e_mode = - log.(exp_minus_e_mode)
    a = log.(1 ./ (n_mode .+ 1))
    b = e_mode .* n_mode
    b[ n_mode .<= n_mode_cutoff ] .= 0
    return sum(b-a)
end

let

    L = 1000
    uniform_field = 1e-7
    boundary_field = 0
    periodic_boundary_conditions = true

    N_samples = 100 #number of distinct disorder realizations

    ### parameters for entanglement calculation
    cut_positions = Vector(10:10:L-1)
    n_mode_cutoff = 1e-7 #ignore modes with occupation smaller than cutoff
    additive_regulator = 1e-8
    # added to the diagonal of the reduced correlator matrix to ensure its positivity
    # (in spite of numerical errors)

    # suffix of the output files
    file_suffix = "_L_$(L)_uni_$(uniform_field)_bound_$(boundary_field).txt"

    # if using open boundaries, remove edge sites when averaging correlators
    if !periodic_boundary_conditions
        edge = L÷5
    end


    rng = MersenneTwister()

    for sample_ind in 1:N_samples

        @show sample_ind

        cs = random_cj_vector(L, rng)

        D = quadratic_bosonic_form(cs, uniform_field, boundary_field,
                periodic_boundary_conditions)

        NA = length(A_sites(cs))
        NB = length(B_sites(cs))
        J0 = diagm([ones(NA); -ones(NB)]) #see Eq. (S64)

        Delta, Tinv = diagonalize_quadratic_form(D, J0)

        #compute the average of Gs Eq. (S69)
        B0 = diagm([zeros(NA); ones(NB)])
        corr_matrix =  Tinv * B0 * Tinv' - B0
        permutation = sortperm([A_sites(cs); B_sites(cs)])
        corr_matrix = corr_matrix[permutation,:]
        corr_matrix = corr_matrix[:,permutation]
        if periodic_boundary_conditions
            Gs_av = diagonal_means(corr_matrix)
        else
            Gs_av = diagonal_means(corr_matrix[edge+1:L-edge,edge+1:L-edge])
        end


        #compute the average of GN Eq. (S70)
        corr_matrix .*= (cs * cs')
        if periodic_boundary_conditions
            GN_av = diagonal_means(corr_matrix)
        else
            GN_av = diagonal_means(corr_matrix[edge+1:L-edge,edge+1:L-edge])
        end
        corr_matrix = log.(corr_matrix)
        if periodic_boundary_conditions
            GN_typ = exp.(diagonal_means(corr_matrix))
        else
            GN_typ = exp.(diagonal_means(corr_matrix[edge+1:L-edge,edge+1:L-edge]))
        end
        # Note: n_j in Eq. (S67) is stored along the diagonal of the corr_matrix
        # and can thus be read from GN or Gs

        # save correlators to file
        open("cs" * file_suffix,"a") do f
            writedlm(f, cs')
        end
        open("GN_typ" * file_suffix,"a") do f
            writedlm(f, GN_typ')
        end
        open("GN_av" * file_suffix,"a") do f
            writedlm(f, GN_av')
        end
        open("Gs_av" * file_suffix,"a") do f
            writedlm(f, GN_av')
        end

        lowest_modes = sort(Delta)[1:5]

        open("modes_energy" * file_suffix,"a") do f
            writedlm(f, lowest_modes')
        end


        #  compute entanglement

        # compute <alpha_i alpha_j^\dag> matrix
        G0 = diagm([ones(NA); zeros(NB)])
        alpha_corr_matrix =  Symmetric(Tinv * G0 * Tinv')
        permutation = sortperm([A_sites(cs); B_sites(cs)])
        alpha_corr_matrix = alpha_corr_matrix[permutation,:]
        alpha_corr_matrix = alpha_corr_matrix[:,permutation]

        J0 = diagm(cs)

        vN_entropy = zeros(length(cut_positions))

        for cut_ind in eachindex(cut_positions)
            cut_position = cut_positions[cut_ind]

            J0_tilde = J0[1:cut_position, 1:cut_position]
            corr_tilde = alpha_corr_matrix[1:cut_position, 1:cut_position]
            #add a small component proportional to the identity to ensure positivity
            corr_tilde += additive_regulator * diagm(ones(cut_position))


            NA_tilde = Int64(sum((cs[1:cut_position] .+ 1) / 2))
            NB_tilde = Int64(cut_position - NA_tilde)
            J0_hat = diagm([ones(NA_tilde); -ones(NB_tilde)])

            n_mode, _ = diagonalize_quadratic_form(corr_tilde, J0_tilde, J0_hat)
            n_mode .-= [ones(NA_tilde); zeros(NB_tilde)]
            vN_entropy[cut_ind] = entropy_from_n(n_mode, n_mode_cutoff)
        end

        open("entropy" * file_suffix,"a") do f
            writedlm(f, vN_entropy')
        end

    end

end