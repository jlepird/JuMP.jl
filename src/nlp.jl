type NLPData
    nlobj
    nlconstr
end

NLPData() = NLPData(nothing, Dict())

function initNLP(m::Model)
    if !haskey(m.ext, :NLP)
        m.ext[:NLP] = NLPData()
    end
end


typealias NonlinearConstraint GenericRangeConstraint{ReverseDiffSparse.SymbolicOutput}

using Ipopt

function solveIpopt(m::Model)
    # check that there are no integer variables
    for j = 1:m.numCols
        if m.colCat[j] == INTEGER
            error("Integer variables present in nonlinear problem")
        end
    end
    tic()
    @assert haskey(m.ext, :NLP)
    nldata = m.ext[:NLP]
    @assert isa(nldata.nlobj, ReverseDiffSparse.SymbolicOutput)
    @assert length(m.obj.qvars1) == 0 && length(m.obj.aff.vars) == 0
    @assert length(m.quadconstr) == 0
    @assert m.objSense == :Min
    fg = genfgrad_simple(nldata.nlobj)

    _, linrowlb, linrowub = prepProblemBounds(m)
    A = prepConstrMatrix(m)

    nlrowlb = Float64[]
    nlrowub = Float64[]
    n_nlconstr = 0
    constr_hessian = Dict()
    constr_grad = Dict()
    nnz_jac = length(A.nzval)
    nnz_hess = 0
    for k in keys(nldata.nlconstr)
        for c in nldata.nlconstr[k]
            n_nlconstr += 1
            push!(nlrowlb, c.lb)
            push!(nlrowub, c.ub)
            nnz_jac += length(c.terms.mapfromcanonical)
        end
        constr_hessian[k] = gen_hessian_sparse_color_parametric(nldata.nlconstr[k][1].terms)
        nnz_hess += length(nldata.nlconstr[k])*length(constr_hessian[k][1])
        constr_grad[k] = genfgrad_parametric(nldata.nlconstr[k][1].terms)
    end
    tprep = toq()
    #println("Prep time: $tprep")
    tmpvec = Array(Float64, m.numCols)
    tf, tgf, tg, tjg, th = zeros(5)
    function eval_f(x)
        #print("x = ");show(x);println()
        #println("f(x) = ", fg(x,tmpvec))
        tic()
        v = fg(x,tmpvec)
        tf += toq()
        return v
    end

    function eval_grad_f(x, g)
        tic()
        fg(x,g)
        tgf += toq()
        #print("x = ");show(x);println()
        #println("gradf(x) = ");show(g);println()
    end

    function eval_g(x, g)
        tic()
        A_mul_B!(sub(g,1:size(A,1)),A,x)
        #g[1:size(A,1)] = A*x
        pos = size(A,1)
        for k in keys(nldata.nlconstr)
            for c in nldata.nlconstr[k]
                pos += 1
                g[pos] = constr_grad[k](x, ReverseDiffSparse.IdentityArray(), tmpvec, ReverseDiffSparse.IdentityArray(), c.terms.inputvals)
            end
        end
        tg += toq()
        #print("x = ");show(x);println()
        #println(size(A,1), " g(x) = ");show(g);println()
    end

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
            row = size(A,1)
            for k in keys(nldata.nlconstr)
                for c in nldata.nlconstr[k]
                    row += 1
                    for v in c.terms.mapfromcanonical
                        rows[idx] = row
                        cols[idx] = v
                        idx += 1
                    end
                end
            end
            @assert idx-1 == nnz_jac
            #print("I ");show(rows);println()
            #print("J ");show(cols);println()


        else
            # Values
            tic()
            fill!(values,0.0)
            idx = 1
            for col = 1:size(A,2)
                for pos = A.colptr[col]:(A.colptr[col+1]-1)
                    values[idx] = A.nzval[pos]
                    idx += 1
                end
            end
            row = size(A,1)
            for k in keys(nldata.nlconstr)
                fgrad = constr_grad[k]
                for c in nldata.nlconstr[k]
                    l = length(c.terms.mapfromcanonical)
                    fgrad(x, ReverseDiffSparse.IdentityArray(), sub(values, idx:(idx+l-1)), c.terms.maptocanonical, c.terms.inputvals)
                    idx += l
                end
            end
            @assert idx-1 == nnz_jac
            tjg += toq()
            #print("x = ");show(x);println()
            #print("V ");show(values);println()
        end
    end

    hI, hJ, hfunc = gen_hessian_sparse_color_parametric(nldata.nlobj)
    nnz_hess += length(hI)

    function eval_h(
        x::Vector{Float64},         # Current solution
        mode,                       # Either :Structure or :Values
        rows::Vector{Int32},        # Sparsity structure - row indices
        cols::Vector{Int32},        # Sparsity structure - column indices
        obj_factor::Float64,        # Lagrangian multiplier for objective
        lambda::Vector{Float64},    # Multipliers for each constraint
        values::Vector{Float64})    # The values of the Hessian

        if mode == :Structure
            for i in 1:length(hI)
                rows[i] = nldata.nlobj.mapfromcanonical[hI[i]]
                cols[i] = nldata.nlobj.mapfromcanonical[hJ[i]]
            end
            idx = length(hI)+1
            for k in keys(nldata.nlconstr)
                lI, lJ, lhess = constr_hessian[k]
                for c in nldata.nlconstr[k]
                    for r in 1:length(lI)
                        rows[idx] = c.terms.mapfromcanonical[lI[r]]
                        cols[idx] = c.terms.mapfromcanonical[lJ[r]]
                        idx += 1
                    end
                end
            end
            @assert idx-1 == nnz_hess
        else
            tic()
            hfunc(x, sub(values, 1:length(hI)), nldata.nlobj)
            scale!(sub(values, 1:length(hI)), obj_factor)
            idx = length(hI)+1
            cnt = 0
            for k in keys(nldata.nlconstr)
                lI, lJ, lhess = constr_hessian[k]
                for c in nldata.nlconstr[k]
                    cnt += 1
                    lhess(x, sub(values, idx:(idx+length(lI)-1)),c.terms)
                    scale!(sub(values, idx:(idx+length(lI)-1)), lambda[cnt])
                    idx += length(lI)
                end
            end
            @assert idx-1 == nnz_hess
            th += toq()

        end
    end
    #print("LB: ");show([linrowlb,nlrowlb]);println()
    #print("UB: ");show([linrowub,nlrowub]);println()
    prob = createProblem(m.numCols, m.colLower, m.colUpper, length(m.linconstr)+n_nlconstr,
        [linrowlb,nlrowlb], [linrowub, nlrowub], nnz_jac, nnz_hess,
        eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h)

    if !any(isnan(m.colVal))
        prob.x = m.colVal
    else
        # solve LP to find feasible point
        # do we need an iterior point?
        lpsol = linprog(zeros(m.numCols), A, linrowlb, linrowub, m.colLower, m.colUpper)
        @assert lpsol.status == :Optimal
        prob.x = lpsol.sol
    end

    status = solveProblem(prob)
    println("STATUS: ", Ipopt.ApplicationReturnStatus[status])
    m.colVal = prob.x
    m.objVal = prob.obj_val

    #println("feval $tf\nfgrad $tgf\ngeval $tg\njaceval $tjg\nhess $th")

end




