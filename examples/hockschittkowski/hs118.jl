using JuMP
using Base.Test

m = Model()

L = zeros(15)
L[1] =  8.0
L[2] = 43.0
L[3] =  3.0

U = zeros(15)
U[1] = 21.0
U[2] = 57.0
U[3] = 16.0
for k in 1:4
    U[3*k+1] =  90.0
    U[3*k+2] = 120.0
    U[3*k+3] =  60.0
end

@defVar(m, L[i] <= x[i=1:15] <= U[i])

@setNLObjective(m, Min, 
    sum{2.3     * x[3*k+1]   +
        0.0001  * x[3*k+1]^2 + 
        1.7     * x[3*k+2]   +
        0.0001  * x[3*k+2]^2 +
        2.2     * x[3*k+3] +
        0.00015 * x[3*k+3]^2, k=0:4})

# setObjective(m, :Min, 
#   sum([(2.3     * x[3*k+1]   +
#         0.0001  * x[3*k+1]^2 + 
#         1.7     * x[3*k+2]   +
#         0.0001  * x[3*k+2]^2 +
#         2.2     * x[3*k+3]^2 +
#         0.00015 * x[3*k+3]^2) for k in 0:4 ]))


# constr1
for j in 1:4
    @addConstraint(m, x[3*j+1] - x[3*j-2] + 7 <= 13)
    @addConstraint(m, x[3*j+1] - x[3*j-2] + 7 >=  0)
end

# constr2
for j in 1:4
    @addConstraint(m, x[3*j+2] - x[3*j-1] + 7 <= 14)
    @addConstraint(m, x[3*j+2] - x[3*j-1] + 7 >=  0)
end

# constr3
for j in 1:4
    @addConstraint(m, x[3*j+3] - x[3*j  ] + 7 <= 13)
    @addConstraint(m, x[3*j+3] - x[3*j  ] + 7 >=  0)
end

@addConstraint(m, x[1] + x[2] + x[3]    >= 60)
@addConstraint(m, x[4] + x[5] + x[6]    >= 50)
@addConstraint(m, x[7] + x[8] + x[9]    >= 70)
@addConstraint(m, x[10] + x[11] + x[12] >= 85)
@addConstraint(m, x[13] + x[14] + x[15] >= 100)

# Initial solution
setValue(x[1], 20.0)
setValue(x[2], 55.0)
setValue(x[3], 15.0)
setValue(x[4], 20.0)
setValue(x[5], 60.0)
setValue(x[6], 20.0)
setValue(x[7], 20.0)
setValue(x[8], 60.0)
setValue(x[9], 20.0)
setValue(x[10], 20.0)
setValue(x[11], 60.0)
setValue(x[12], 20.0)
setValue(x[13], 20.0)
setValue(x[14], 60.0)
setValue(x[15], 20.0)

JuMP.solveIpopt(m)
# solve(m)

println(getValue(x))

@test_approx_eq_eps getValue(x[1]) 8.0  1e-5
@test_approx_eq_eps getValue(x[2]) 49.0 1e-5
@test_approx_eq_eps getValue(x[3]) 3.0  1e-5
@test_approx_eq_eps getValue(x[4]) 1.0  1e-5
