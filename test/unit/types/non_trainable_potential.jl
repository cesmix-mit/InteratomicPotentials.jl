@testset "Non-Trainable Potential Default Implmentation Unit Tests" begin
    p = MockNonTrainablePotential()

    @test get_parameters(p) == (;)
    @test set_parameters(p, (;)) === p
    @test serialize_parameters(p) == []
    @test deserialize_parameters(p, []) === p

    @test get_hyperparameters(p) == (;)
    @test set_hyperparameters(p, (;)) === p
    @test serialize_hyperparameters(p) == []
    @test deserialize_hyperparameters(p, []) === p
end
