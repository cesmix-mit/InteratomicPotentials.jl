@testset "Trainable Potential Default Implmentation Unit Tests" begin
    p = MockTrainablePotential(1.0, 2.0, 3.0)

    @test get_parameters(p) == (; A=1.0, B=2.0)
    @test set_parameters(p, (; A=2.0, B=4.0)) == MockTrainablePotential(2.0, 4.0, 3.0)
    @test serialize_parameters(p) == [1.0, 2.0]
    @test deserialize_parameters(p, [2.0, 4.0]) == MockTrainablePotential(2.0, 4.0, 3.0)

    @test get_hyperparameters(p) == (; C=3.0)
    @test set_hyperparameters(p, (; C=6.0)) == MockTrainablePotential(1.0, 2.0, 6.0)
    @test serialize_hyperparameters(p) == [3.0]
    @test deserialize_hyperparameters(p, [6.0]) == MockTrainablePotential(1.0, 2.0, 6.0)
end
