@info "Preparing environment"
import Pkg
Pkg.activate(".")
Pkg.instantiate()

import CSV
using DataPipes
using FlexiMaps
using StaticArrays
using StatsBase: mean, geomean
using DirectionalStatistics
using Accessors
using ThreadsX


function main()
    @info "Reading data"
    events = @p	CSV.File("icecube_he_events.csv") |>
        map((_.datetime, coords=CoordsEllipse(
            (_.ra, _.dec),
            1.3 .* (_.ra_cderr⁻, _.ra_cderr⁺),
            1.3 .* (_.dec_err⁻, _.dec_err⁺),
        ),))
    rfc = @p CSV.File("rfc.csv") |> map((_.name, _.flux, coords=(_.ra, _.dec)))

    @info "Computing simulations"
    computed = map([0, 0.45]) do add°
        real_avg = average_flux_nearby(events, rfc, add°)
        rand_avgs = ThreadsX.map(1:5*10^4) do _
            average_flux_nearby(randomize(events), rfc, add°)
        end
        (; add°, real_avg, rand_avgs)
    end

    @info "P-values:"
    map(computed) do c
        @info "" c.add° p=mean([c.real_avg; c.rand_avgs] .>= c.real_avg)
    end
end


struct CoordsEllipse
    radec::NTuple{2, Float64}
    racos⁻⁺::NTuple{2, Float64}
    dec⁻⁺::NTuple{2, Float64}
end

average_flux_nearby(evts, rfc, add°) = @p let
    evts
    flatmap() do evt
        @p rfc |>
            filter(separation_minus_error(evt.coords, _.coords) <= deg2rad(add°))
    end
    map(_.flux)
    geomean
end

randomize(evts) = @p evts |> map(@set _.coords.radec[1] = 2π * rand())

function separation_minus_error(x::CoordsEllipse, y::NTuple{2, Float64})
    Δ = separation_relradec(x.radec,  y)
    a = Δ.Δracosδ > 0 ? x.racos⁻⁺[2] : x.racos⁻⁺[1]
    b = Δ.Δdec > 0 ? x.dec⁻⁺[2] : x.dec⁻⁺[1]
    return distance_to_ellipse(SVector(a, b), SVector(Δ...))
end

function separation_relradec(x, y)
    cosδ = cos((y[2] + x[2]) / 2)
    Δracosδ = Circular.center_angle(y[1] - x[1]) * cosδ
    Δdec = y[2] - x[2]
    (;Δracosδ, Δdec)
end

using LinearAlgebra: norm_sqr, normalize, norm

# taken/adapted from https://github.com/0xfaded/ellipse_demo/issues/1
function distance_to_ellipse(ab, xy, maxiter=2)
    xyabs = abs.(xy)
    ab² = ab.^2
    ab⁻¹ = 1 ./ ab
    abp = SVector(ab².x - ab².y, ab².y - ab².x) .* ab⁻¹
    T = SVector(1 / √2, 1 / √2)
    for _ in 1:maxiter
        exy = abp .* T.^3
        qxy = xyabs - exy
        rqxy = qxy * sqrt(norm_sqr(ab .* T - exy) / norm_sqr(qxy))
        T = max.((rqxy + exy) .* ab⁻¹, 0) |> normalize
    end
    nearest_abs = ab .* T
    dist = norm(xyabs - nearest_abs)
    return norm_sqr(xyabs) < norm_sqr(nearest_abs) ? -dist : dist
end

main()
