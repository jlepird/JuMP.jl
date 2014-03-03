#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# Hock-Schittkowski Nonlinear Test Suite  -  HS110
# This file is JuMP implementation of the model described in 
#  W. Hock, K. Schittkowski, Test Examples for Nonlinear Programming
#  Codes, Lecture Notes in Economics and Mathematical Systems, 
#  Springer, No, 187, 1981 
# More information, including original model description, at
# http://www.ai7.uni-bayreuth.de/downloads.htm
#
# This problem has an objective with squared logarithms and a product
# of variables squared.
#############################################################################

using JuMP
using Base.Test

m = Model()
@defVar(m, -2.001 <= x[1:10] <= 9.999)
for i = 1:10
    setValue(x[i], 9)
end

@setNLObjective(m, Min, 
    sum{ log(x[j] - 2)^2 + log(10 - x[j])^2, j=1:10} -
    (x[1]*x[2]*x[3]*x[4]*x[5]*x[6]*x[7]*x[8]*x[9]*x[10]) ^ 0.2
)

JuMP.solveIpopt(m)

@test_approx_eq_eps getObjectiveValue(m) -45.77846971 1e-5