using Optim

function rosenbrock(x::Vector)

    return (20 - x[1])^2 + (15 - x[2])^2 + (1e-4 - x[3])^2
end

res = optimize(rosenbrock,[0.0, 0.0, 0.0])

show(res)