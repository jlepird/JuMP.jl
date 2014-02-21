using Ipopt

function solveIpopt(m::Model)
    # check that there are no integer variables
    for j = 1:m.numCols
        if m.colCat[j] == INTEGER
            error("Integer variables present in nonlinear problem")
        end
    end

    @assert isa(m.nlobj, ReverseDiffSparse.SymbolicOutput)
    @assert length(m.obj.qvars1) == 0 && length(m.obj.aff.vars) == 0
    @assert length(m.quadconstr) == 0
    @assert m.objSense == :Min
    fg = genfgrad_simple(m.nlobj)

    _, rowlb, rowub = prepProblemBounds(m)
    A = prepConstrMatrix(m)

    tmpvec = Array(Float64, m.numCols)
    function eval_f(x)
        return fg(x,tmpvec)
    end

    function eval_grad_f(x, g)
        fg(x,g)
    end

    function eval_g(x, g)
        A_mul_B!(g,A,x)
    end

    # Copied from IpoptSolverInterface.jl
    function eval_jac_g(x, mode, rows, cols, values)
        if mode == :Structure
            # Convert column wise sparse to triple format
            idx = 1
            for col = 1:size(A,2)
                for pos = A.colptr[col]:(A.colptr[col+1]-1)
                    rows[idx] = A.rowval[pos]
                    cols[idx] = col
                    idx += 1
                end
            end
        else
            # Values
            idx = 1
            for col = 1:n
                for pos = A.colptr[col]:(A.colptr[col+1]-1)
                    values[idx] = A.nzval[pos]
                    idx += 1
                end
            end
        end
    end

    hmat, hfunc = gen_hessian_sparse_mat(m.nlobj)

    function eval_h(
        x::Vector{Float64},         # Current solution
        mode,                       # Either :Structure or :Values
        rows::Vector{Int32},        # Sparsity structure - row indices
        cols::Vector{Int32},        # Sparsity structure - column indices
        obj_factor::Float64,        # Lagrangian multiplier for objective
        lambda::Vector{Float64},    # Multipliers for each constraint
        values::Vector{Float64})    # The values of the Hessian

        if mode == :Structure
            # Convert column wise sparse to triple format
            idx = 1
            for col = 1:size(hmat,2)
                for pos = hmat.colptr[col]:(hmat.colptr[col+1]-1)
                    rows[idx] = hmat.rowval[pos]
                    cols[idx] = col
                    idx += 1
                end
            end
        else
            hmat.nzval = values
            hfunc(x, hmat)
            scale!(hmat.nzval, obj_factor)
        end
    end

    prob = createProblem(m.numCols, m.colLower, m.colUpper, length(m.linconstr),
        rowlb, rowub, length(A.nzval), length(hmat.nzval),
        eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)

    if !any(isnan(m.colVal))
        prob.x = m.colVal
    else
        # solve LP to find feasible point
        # do we need an iterior point?
        lpsol = linprog(zeros(m.numCols), A, rowlb, rowub, m.colLower, m.colUpper)
        @assert lpsol.status == :Optimal
        prob.x = lpsol.sol
    end

    status = solveProblem(prob)
    println("STATUS: ", Ipopt.ApplicationReturnStatus[status])
    m.colVal = prob.x
    m.objVal = prob.obj_val

end




