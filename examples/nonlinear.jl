using JuMP

# rosenbrock

let

    m = Model()

    @defVar(m, x)
    @defVar(m, y)

    @setNLObjective(m, Min, (1-x)^2 + 100(y-x^2)^2)

    JuMP.solveIpopt(m)

    println("x = ", getValue(x), " y = ", getValue(y))

end

# clnlbeam
let
    N     = 1000
    ni    = N
    h     = 1/ni
    alpha = 350

    m = Model()

    @defVar(m, -1 <= t[1:(ni+1)] <= 1)
    @defVar(m, -0.05 <= x[1:(ni+1)] <= 0.05)
    @defVar(m, u[1:(ni+1)])

    @setNLObjective(m, Min, sum{ 0.5*h*(u[i+1]^2 + u[i]^2) + 0.5*alpha*h*(cos(t[i+1]) + cos(t[i])), i = 1:ni})
    
    #cons = Array(NLExpr,0)
    # cons1
    #for i in 1:ni
    #    @addNLConstr(cons, x[i+1] - x[i] - (0.5h)*(sin(t[i+1])+sin(t[i])) == 0)
    #end
    # cons2
    #for i in 1:ni
    #    @addNLConstr(cons, t[i+1] - t[i] - (0.5h)*u[i+1] - (0.5h)*u[i] == 0)
    #end

    JuMP.solveIpopt(m)

end


# hs071
let
    # min x1 * x4 * (x1 + x2 + x3) + x3
    # st  x1 * x2 * x3 * x4 >= 25
    #     x1^2 + x2^2 + x3^2 + x4^2 = 40
    #     1 <= x1, x2, x3, x4 <= 5
    # Start at (1,5,5,1)
    # End at (1.000..., 4.743..., 3.821..., 1.379...)

    m = Model()

    @defVar(m, 1 <= x[1:4] <= 5)

    @setNLObjective(m, Min, x[1]*x[4]*(x[1]+x[2]+x[3]) + x[3])

    #@addNLConstr(m??, x[1]*x[2]*x[3]*x[4] >= 25)
    #@addNLConstr(m??, sum{x[i]^2,i=1:4} == 40)

    JuMP.solveIpopt(m)

    #getValue(x[1]) == 1.00000000
    #getValue(x[2]) == 4.74299963
    #getValue(x[3]) == 3.82114998
    #getValue(x[4]) == 1.37940829
end