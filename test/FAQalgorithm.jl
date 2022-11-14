function check_doubly_stochastic(D::AbstractMatrix{U}) where {U<:Real}
    n=size(D)[1]
    for i in 1:n
        @test 0.99<sum(D[i,:])<1.01 #check sum of elements in a row == 1
        @test 0.99<sum(D[:,i])<1.01 #check sum of elements in a column == 1
    end
end

function check_permutation_matrix(A::AbstractMatrix{U}) where {U<:Real}
    n=size(A)[1]
    for j in 1:n
        for i in 1:n
            @test (A[i,j]==one(U) || A[i,j]==zero(U)) #check that entries are 1 or 0
        end
    end
    check_doubly_stochastic(A) #check that rows and colums sum to 1
end

@testset "FAQ algorithm" begin
    I3=[1 0 0; 0 1 0; 0 0 1]
    @testset "flat doubly stochastic"   begin
        n=5
        DS=_flat_doubly_stochastic(n)
        @test size(DS)==(n,n) #check dimensions
        check_doubly_stochastic(DS)
    end;

    @testset "gradient" begin
        A=rand(5,5)
        A=A+A'
        B=rand(5,5)
        B=B+B'
        P=float.(Matrix(I(5))) #float because matrix of bool
        @test _gradient(A,B,P)==-2*A*P*B
        @test _gradient(zeros(5,5),zeros(5,5),zeros(5,5))==zeros(5,5)
    end;

    @testset "distance" begin
        A=[0 1 1; 1 0 1; 1 1 0]
        B=[0 1 0; 1 0 1; 0 1 0]
        @test _distance(A,B,I3)==sqrt(2)
        @test _distance(A,A,I3)==0
    end;

    @testset "transportation problem" begin
        A=[1 4 5; 6 2 8; 9 9 3]
        @test _solve_transportation_problem(A,optimizer=GLPK.Optimizer)==I3
        check_doubly_stochastic(_solve_transportation_problem(rand(5,5);optimizer=GLPK.Optimizer))
    end;

    @testset "assignment problem" begin
        A=[7 4 5; 6 9 8; 9 9 11]
        @test _solve_assignment_problem(A,optimizer=GLPK.Optimizer)==I3
        check_permutation_matrix(_solve_assignment_problem(rand(5,5),optimizer=GLPK.Optimizer))
    end;

    @testset "update P" begin
        P=[0 1 1; 1 0 1; 1 1 0]
        Q=[0 1 0; 1 0 1; 0 1 0]
        @test _update_P(P,Q,1.)==P
        @test _update_P(P,Q,0.)==Q
    end;

    @testset "FAQ" begin
        P,_,_=faq(rand(5,5),rand(5,5),optimizer=GLPK.Optimizer)
        check_doubly_stochastic(P)
        @testset "empty graphs" begin
            A=adjacency_matrix(SimpleGraph(5))
            B=adjacency_matrix(SimpleGraph(5))
            _,norm,converged=faq(A,B,optimizer=GLPK.Optimizer)
            @test norm==0.
            @test converged
        end
        @testset "path and its permuted" begin
            A=adjacency_matrix(path_graph(5))
            P=zeros(5,5)
            P[1,1]=1
            P[2,3]=1
            P[3,2]=1
            P[4,5]=1
            P[5,4]=1
            B=P*A*P'
            P_res,norm,converged=faq(A,B,optimizer=HiGHS.Optimizer)
            @test P_res==P
            @test norm==0.
            @test converged
        end
        @testset "two path and complete graph" begin
            A=adjacency_matrix(path_graph(3))
            B=adjacency_matrix(complete_graph(3))
            _,norm,converged=faq(A,B,optimizer=HiGHS.Optimizer)
            @test norm==sqrt(2)
            @test converged
        end
    end

end;
