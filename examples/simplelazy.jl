#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# simplelazy.jl
#
# Solve a simple problem using lazy constraints. See documentation for more
# details about callbacks.
#############################################################################

using JuMP
using Gurobi

# We will use Gurobi
m = Model(solver=GurobiSolver())

# Define our variables to be inside a box, and integer
@defVar(m, 0 <= x <= 2, Int)
@defVar(m, 0 <= y <= 2, Int)

@setObjective(m, Max, y)

# We now define our callback function that takes one argument,
# the callback handle. Note that we can access m, x, and y because
# this function is defined inside the same scope
function corners(cb)
    x_val = getValue(x)
    y_val = getValue(y)
    println("In callback function, x=$x_val, y=$y_val")

    # We have two constraints, one cutting off the top
    # left corner and one cutting off the top right corner, e.g.
    # (0,2) +---+---+ (2,2)
    #       |xx/ \xx|
    #       |x/   \x|
    #       |/     \|
    #       +       +
    #       |       |
    #       |       |
    #       |       |
    # (0,0) +---+---+ (2,0)

    # Allow for some impreciseness in the solution
    TOL = 1e-6

    # Check top left, allowing some tolerance
    if y_val - x_val > 1 + TOL
        # Cut off this solution
        println("Solution was in top left, cut it off")
        # Use the original variables
        @addLazyConstraint(cb, y - x <= 1)
    # Check top right
    elseif y_val + x_val > 3 + TOL
        # Cut off this solution
        println("Solution was in top right, cut it off")
        # Use the original variables
        @addLazyConstraint(cb, y + x <= 3)
    end
end  # End of callback function

# Tell JuMP/Gurobi to use our callback function
setLazyCallback(m, corners)

# Solve the problem
solve(m)

# Print our final solution
println("Final solution: [ $(getValue(x)), $(getValue(y)) ]")
