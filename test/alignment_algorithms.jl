function check_doubly_stochastic(D::AbstractMatrix{U}) where {U<:Real}
    n = size(D)[1]
    for i in 1:n
        @test sum(D[i, :]) ≈ 1  # check sum of elements in a row == 1
        @test sum(D[:, i]) ≈ 1 # check sum of elements in a column == 1
    end
end

function check_permutation_matrix(A::AbstractMatrix{U}) where {U<:Real}
    n = size(A)[1]
    for j in 1:n
        for i in 1:n
            @test (A[i, j] == one(U) || A[i, j] == zero(U)) # check that entries are 1 or 0
        end
    end
    return check_doubly_stochastic(A) # check that rows and columns sum to 1
end

@testset "alignment algorithms" begin
    I3 = [1 0 0; 0 1 0; 0 0 1]
    @testset "flat doubly stochastic" begin
        n = 5
        DS = GraphsOptim._flat_doubly_stochastic(n)
        @test size(DS) == (n, n) # check dimensions
        check_doubly_stochastic(DS)
    end

    @testset "gradient" begin
        A = rand(5, 5)
        A = A + A'
        B = rand(5, 5)
        B = B + B'
        P = float.(Matrix(I(5))) # float because matrix of bool
        @test GraphsOptim._gradient(A, B, P) == -2 * A * P * B
        @test GraphsOptim._gradient(zeros(5, 5), zeros(5, 5), zeros(5, 5)) == zeros(5, 5)
    end

    @testset "distance" begin
        A = [0 1 1; 1 0 1; 1 1 0]
        B = [0 1 0; 1 0 1; 0 1 0]
        @test GraphsOptim._distance(A, B, I3) == sqrt(2)
        @test GraphsOptim._distance(A, A, I3) == 0
    end

    @testset "assignment problem" begin
        A = [7 4 5; 6 9 8; 9 9 11]
        @test GraphsOptim._solve_assignment_problem(A; optimizer=HiGHS.Optimizer) == I3
        check_permutation_matrix(
            GraphsOptim._solve_assignment_problem(rand(5, 5); optimizer=HiGHS.Optimizer)
        )
    end

    @testset "update P" begin
        P = [0 1 1; 1 0 1; 1 1 0]
        Q = [0 1 0; 1 0 1; 0 1 0]
        @test GraphsOptim._update_P(P, Q, 1.0) == P
        @test GraphsOptim._update_P(P, Q, 0.0) == Q
    end

    @testset "FAQ" begin
        P, _, _ = faq(rand(5, 5), rand(5, 5); optimizer=HiGHS.Optimizer)
        check_doubly_stochastic(P)
        @testset "empty graphs" begin
            A = adjacency_matrix(SimpleGraph(5))
            B = adjacency_matrix(SimpleGraph(5))
            _, norm, converged = faq(A, B; optimizer=HiGHS.Optimizer)
            @test norm == 0.0
            @test converged
        end

        @testset "two path and complete graph" begin
            A = adjacency_matrix(path_graph(3))
            B = adjacency_matrix(complete_graph(3))
            _, norm, converged = faq(A, B; optimizer=HiGHS.Optimizer)
            @test norm == sqrt(2)
            @test converged
        end
    end
end;
