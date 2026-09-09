using Test
using NeutralAtoms

function motion_error_options(;
    atom_motion=true,
    free_motion=false,
    xy_motion=true,
    z_motion=true,
    Doppler=true,
)
    return Dict(
        "laser_noise" => false,
        "spontaneous_decay_intermediate" => false,
        "spontaneous_decay_rydberg" => false,
        "blockade"=>true,
        "atom_motion" => atom_motion,
        "free_motion" => free_motion,
        "xy_motion" => xy_motion,
        "z_motion" => z_motion,
        "Doppler" => Doppler,
    )
end

@testset "Motion regressions" begin
    trap_params = [500.0, 1.2, w0_to_z0(1.2, 0.813)]
    atom_params = [86.9, 10.0]
    ωr, ωz = trap_frequencies(atom_params, trap_params)

    @testset "Zero-temperature sampling" begin
        zero_temp_atom_params = [86.9, 0.0]
        samples, acc_rate = samples_generate(trap_params, zero_temp_atom_params, 3; harmonic=true)

        @test acc_rate == 1
        @test length(samples) == 3
        @test all(sample -> all(iszero, sample), samples)
    end

    @testset "Single-atom static shifts survive disabled axes" begin
        sample = [0.01, -0.02, 0.05, 0.1, -0.08, 0.03]
        center = [1.5, -0.25, 0.4]

        xy_off_opts = motion_error_options(; xy_motion=false, z_motion=true)
        X, Y, Z, Vx, Vy, Vz = NeutralAtoms.get_atom_trajectories(sample, center, ωr, ωz, xy_off_opts)

        @test X(0.0) == center[1]
        @test Y(0.0) == center[2]
        @test Z(0.0) ≈ center[3] + sample[3]

        z_off_opts = motion_error_options(; xy_motion=true, z_motion=false)
        X, Y, Z, Vx, Vy, Vz = NeutralAtoms.get_atom_trajectories(sample, center, ωr, ωz, z_off_opts)

        @test X(0.0) ≈ center[1] + sample[1]
        @test Y(0.0) ≈ center[2] + sample[2]
        @test Z(0.0) == center[3]

        static_opts = motion_error_options(; atom_motion=false, xy_motion=true, z_motion=true)
        X, Y, Z, Vx, Vy, Vz = NeutralAtoms.get_atom_trajectories(sample, center, ωr, ωz, static_opts)

        @test X(0.0) == center[1]
        @test Y(0.0) == center[2]
        @test Z(0.0) == center[3]
        @test Vx(0.0) == 0.0
        @test Vy(0.0) == 0.0
        @test Vz(0.0) == 0.0
    end

    @testset "CZ blockade keeps trap separation when xy motion is disabled" begin
        sample1 = [0.01, -0.02, 0.05, 0.1, -0.08, 0.03]
        sample2 = [-0.015, 0.01, -0.04, -0.09, 0.07, -0.02]
        center1 = [-1.35, 0.0, 0.0]
        center2 = [1.35, 0.0, 0.0]
        c6 = 1.0

        xy_on_opts = motion_error_options(; xy_motion=true, z_motion=true)
        xy_off_opts = motion_error_options(; xy_motion=false, z_motion=true)

        V_true = NeutralAtoms.get_V(sample1, sample2, center1, center2, ωr, ωz, xy_on_opts, c6)
        V_false = NeutralAtoms.get_V(sample1, sample2, center1, center2, ωr, ωz, xy_off_opts, c6)

        X1, Y1, Z1, Vx1, Vy1, Vz1 = NeutralAtoms.get_atom_trajectories(sample1, center1, ωr, ωz, xy_off_opts)
        X2, Y2, Z2, Vx2, Vy2, Vz2 = NeutralAtoms.get_atom_trajectories(sample2, center2, ωr, ωz, xy_off_opts)

        @test X1(0.0) == center1[1]
        @test X2(0.0) == center2[1]
        @test X1(0.0) != X2(0.0)
        @test Y1(0.0) == center1[2]
        @test Y2(0.0) == center2[2]
        @test Vx1(0.0) == 0.0
        @test Vy1(0.0) == 0.0
        @test Vx2(0.0) == 0.0
        @test Vy2(0.0) == 0.0
        @test Z1(0.0) ≈ center1[3] + sample1[3]
        @test Z2(0.0) ≈ center2[3] + sample2[3]
        @test Vz1(0.0) ≈ sample1[6]
        @test Vz2(0.0) ≈ sample2[6]

        expected_r2 = (center1[1] - center2[1])^2 +
                      (center1[2] - center2[2])^2 +
                      (center1[3] + sample1[3] - center2[3] - sample2[3])^2

        @test isfinite(V_false(0.0))
        @test V_false(0.0) ≈ c6 / (1e-18 + expected_r2^3)
        @test 0.1 <= V_false(0.0) / V_true(0.0) <= 10.0
    end
end
