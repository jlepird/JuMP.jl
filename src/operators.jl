#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################

# Overloads
#
# Different objects that must all interact:
# 1. Number
# 2. Variable
# 3. AffExpr
# 4. QuadExpr

# Number
# Number--Number obviously already taken care of!
# Number--Variable
(+)(lhs::Number, rhs::Variable) = AffExpr([rhs],[+1.],convert(Float64,lhs))
(-)(lhs::Number, rhs::Variable) = AffExpr([rhs],[-1.],convert(Float64,lhs))
(*)(lhs::Number, rhs::Variable) = AffExpr([rhs],[convert(Float64,lhs)], 0.)
(/)(lhs::Number, rhs::Variable) = error("Cannot divide by variable")
# Number--GenericAffExpr
(+)(lhs::Number, rhs::AffExpr)  = AffExpr(copy(rhs.vars),copy(rhs.coeffs),lhs+rhs.constant)
(-)(lhs::Number, rhs::AffExpr)  = AffExpr(copy(rhs.vars),    -rhs.coeffs ,lhs-rhs.constant)
(*)(lhs::Number, rhs::AffExpr)  = AffExpr(copy(rhs.vars), lhs*rhs.coeffs ,lhs*rhs.constant)
(/)(lhs::Number, rhs::AffExpr)  = error("Cannot divide by an affine expression")
# Number--QuadExpr
(+)(lhs::Number, rhs::QuadExpr) = QuadExpr(copy(rhs.qvars1),copy(rhs.qvars2),copy(rhs.qcoeffs),lhs+rhs.aff)
(-)(lhs::Number, rhs::QuadExpr) = QuadExpr(copy(rhs.qvars1),copy(rhs.qvars2),    -rhs.qcoeffs ,lhs-rhs.aff)
(*)(lhs::Number, rhs::QuadExpr) = QuadExpr(copy(rhs.qvars1),copy(rhs.qvars2), lhs*rhs.qcoeffs ,lhs*rhs.aff)
(/)(lhs::Number, rhs::QuadExpr) = error("Cannot divide by a quadratic expression")
# Number--MatrixVar
(+)(lhs::Number, rhs::MatrixVar) = error("Cannot add a scalar and a matrix variable")
(-)(lhs::Number, rhs::MatrixVar) = error("Cannot subtract a scalar and a matrix variable")
(*)(lhs::Number, rhs::MatrixVar) = MatrixExpr([rhs], lhs*Array[eye(rhs)], zero(rhs))
(/)(lhs::Number, rhs::MatrixVar) = error("Cannot divide a scalar by a matrix variable")
# Number--MatrixExpr
(+)(lhs::Number, rhs::MatrixExpr)  = error("Cannot add a scalar to a matrix expression")
(-)(lhs::Number, rhs::MatrixExpr)  = error("Cannot subtract a matrix expression from a number")
(*)(lhs::Number, rhs::MatrixExpr)  = MatrixExpr(copy(rhs.vars), lhs*rhs.coeffs ,lhs*rhs.constant)
(/)(lhs::Number, rhs::MatrixExpr)  = error("Cannot divide a scalar by a matrix expression")

# Variable
# Variable--Number
(+)(lhs::Variable, rhs::Number) = (+)( rhs,lhs)
(-)(lhs::Variable, rhs::Number) = (+)(-rhs,lhs)
(*)(lhs::Variable, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::Variable, rhs::Number) = (*)(1./rhs,lhs)
# Variable--Variable
(+)(lhs::Variable, rhs::Variable) = AffExpr([lhs,rhs], [1.,+1.], 0.)
(-)(lhs::Variable, rhs::Variable) = AffExpr([lhs,rhs], [1.,-1.], 0.)
(*)(lhs::Variable, rhs::Variable) = QuadExpr([lhs],[rhs],[1.],AffExpr(Variable[],Float64[],0.))
(/)(lhs::Variable, rhs::Variable) = error("Cannot divide a variable by a variable")
# Variable--AffExpr
(+)(lhs::Variable, rhs::AffExpr) = AffExpr(vcat(rhs.vars,lhs),vcat( rhs.coeffs,1.), rhs.constant)
(-)(lhs::Variable, rhs::AffExpr) = AffExpr(vcat(rhs.vars,lhs),vcat(-rhs.coeffs,1.),-rhs.constant)
function (*)(lhs::Variable, rhs::AffExpr)
    n = length(rhs.vars)
    if rhs.constant != 0.      
        ret = QuadExpr([lhs for i=1:n],copy(rhs.vars),copy(rhs.coeffs),AffExpr([lhs], [rhs.constant], 0.))
    else
        ret = QuadExpr([lhs for i=1:n],copy(rhs.vars),copy(rhs.coeffs),AffExpr())
    end
end
(/)(lhs::Variable, rhs::AffExpr) = error("Cannot divide a variable by an affine expression")
# Variable--QuadExpr
(+)(v::Variable, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),v+q.aff)
(-)(v::Variable, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),    -q.qcoeffs ,v-q.aff)
(*)(v::Variable, q::QuadExpr) = error("Cannot multiply a variable by a quadratic expression")
(/)(v::Variable, q::QuadExpr) = error("Cannot divide a variable by a quadratic expression")
# Variable--MatrixVar
(+)(lhs::Variable, rhs::MatrixVar) = error("Cannot add a scalar variable and a matrix variable") #TODO: make this work w/ 1x1 matrices
(-)(lhs::Variable, rhs::MatrixVar) = error("Cannot subtract a matrix variable from a scalar variable")
(*)(lhs::Variable, rhs::MatrixVar) = error("Cannot multiply a scalar variable and a matrix variable")
(/)(lhs::Variable, rhs::MatrixVar) = error("Cannot divide a scalar variable by a matrix variable")

# Variable--MatrixExpr
(+)(lhs::Variable, rhs::MatrixExpr) = error("Cannot add a scalar variable and a matrix expression") #TODO: make this work w/ 1x1 matrices
(-)(lhs::Variable, rhs::MatrixExpr) = error("Cannot subtract a matrix expression from a scalar variable")
(*)(lhs::Variable, rhs::MatrixExpr) = error("Cannot multiply a scalar variable and a matrix expression")
(/)(lhs::Variable, rhs::MatrixExpr) = error("Cannot divide a scalar variable by a matrix expression")

# AffExpr
# AffExpr--Number
(+)(lhs::GenericAffExpr, rhs::Number) = (+)(+rhs,lhs)
(-)(lhs::GenericAffExpr, rhs::Number) = (+)(-rhs,lhs)
(*)(lhs::GenericAffExpr, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::GenericAffExpr, rhs::Number) = (*)(1.0/rhs,lhs)
# AffExpr--Variable
(+)(lhs::AffExpr, rhs::Variable) = (+)(rhs,lhs)
(-)(lhs::AffExpr, rhs::Variable) = AffExpr(vcat(lhs.vars,rhs),vcat(+lhs.coeffs,-1.),lhs.constant)
(*)(lhs::AffExpr, rhs::Variable) = (*)(rhs,lhs)
(/)(lhs::AffExpr, rhs::Variable) = error("Cannot divide affine expression by a variable")
# AffExpr--AffExpr
(+)(lhs::GenericAffExpr, rhs::GenericAffExpr) = GenericAffExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs, rhs.coeffs),lhs.constant+rhs.constant)
(-)(lhs::GenericAffExpr, rhs::GenericAffExpr) = GenericAffExpr(vcat(lhs.vars,rhs.vars),vcat(lhs.coeffs,-rhs.coeffs),lhs.constant-rhs.constant)
function (*)(lhs::AffExpr, rhs::AffExpr)
    ret = QuadExpr(Variable[],Variable[],Float64[],AffExpr(Variable[],Float64[],0.))

    # Quadratic terms
    n = length(lhs.coeffs)
    m = length(rhs.coeffs)
    sizehint(ret.qvars1, n*m)
    sizehint(ret.qvars2, n*m)
    sizehint(ret.qcoeffs, n*m)
    for i = 1:n
        for j = 1:m
            push!(ret.qvars1,  lhs.vars[i])
            push!(ret.qvars2,  rhs.vars[j])
            push!(ret.qcoeffs, lhs.coeffs[i]*rhs.coeffs[j])
        end
    end
    
    # Try to preallocate space for aff
    if lhs.constant != 0 && rhs.constant != 0
        sizehint(ret.aff.vars, n+m)
        sizehint(ret.aff.coeffs, n+m)
    elseif lhs.constant != 0
        sizehint(ret.aff.vars, n)
        sizehint(ret.aff.coeffs, n)
    elseif rhs.constant != 0
        sizehint(ret.aff.vars, m)
        sizehint(ret.aff.coeffs, m)
    end

    # [LHS constant] * RHS
    if lhs.constant != 0
        c = lhs.constant
        for j = 1:m
            push!(ret.aff.vars,   rhs.vars[j])
            push!(ret.aff.coeffs, rhs.coeffs[j] * c)
        end
        ret.aff.constant += c * rhs.constant
    end
    
    # Expr 2 constant * Expr 1 terms
    if rhs.constant != 0
        c = rhs.constant
        for i = 1:m
            push!(ret.aff.vars,   lhs.vars[i])
            push!(ret.aff.coeffs, lhs.coeffs[i] * c)
        end
        # Don't do the following line
        #ret.aff.constant += c * lhs.constant
        # If lhs.constant is 0, its a waste of time
        # If lhs.constant is non-zero, its already done
    end
    
    return ret
end
(/)(lhs::AffExpr, rhs::AffExpr) = error("Cannot divide aff. expression by aff. expression")
# AffExpr--QuadExpr
(+)(a::AffExpr, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),a+q.aff)
(-)(a::AffExpr, q::QuadExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),    -q.qcoeffs ,a-q.aff)
(*)(a::AffExpr, q::QuadExpr) = error("Cannot multiply an aff. expression by a quadratic expression")
(/)(a::AffExpr, q::QuadExpr) = error("Cannot divide an aff. expression by a quadratic expression")
# AffExpr--MatrixVar
(+)(lhs::AffExpr, rhs::MatrixVar) = error("Cannot add a scalar expression and a matrix variable") #TODO: make this work w/ 1x1 matrices
(-)(lhs::AffExpr, rhs::MatrixVar) = error("Cannot subtract a matrix variable from a scalar expression")
(*)(lhs::AffExpr, rhs::MatrixVar) = error("Cannot multiply a scalar expression and a matrix variable")
(/)(lhs::AffExpr, rhs::MatrixVar) = error("Cannot divide a scalar expression by a matrix variable")
# AffExpr--MatrixExpr
(+)(lhs::AffExpr, rhs::MatrixExpr) = error("Cannot add a scalar expression and a matrix expression") #TODO: make this work w/ 1x1 matrices
(-)(lhs::AffExpr, rhs::MatrixExpr) = error("Cannot subtract a matrix expression from a scalar expression")
(*)(lhs::AffExpr, rhs::MatrixExpr) = error("Cannot multiply a scalar expression and a matrix expression")
(/)(lhs::AffExpr, rhs::MatrixExpr) = error("Cannot divide a scalar expression by a matrix expression")

# QuadExpr
# QuadExpr--Number
(+)(lhs::QuadExpr, rhs::Number) = (+)(+rhs,lhs)
(-)(lhs::QuadExpr, rhs::Number) = (+)(-rhs,lhs)
(*)(lhs::QuadExpr, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::QuadExpr, rhs::Number) = (*)(1.0/rhs,lhs)
# QuadExpr--Variable
(+)(q::QuadExpr, v::Variable) = (+)(v,q)
(-)(q::QuadExpr, v::Variable) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),q.aff-v)
(*)(q::QuadExpr, v::Variable) = error("Cannot multiply a quadratic expression by a variable")
(/)(q::QuadExpr, v::Variable) = error("Cannot divide a quadratic expression by a variable")
# QuadExpr--AffExpr
(+)(q::QuadExpr, a::AffExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),q.aff+a)
(-)(q::QuadExpr, a::AffExpr) = QuadExpr(copy(q.qvars1),copy(q.qvars2),copy(q.qcoeffs),q.aff-a)
(*)(q::QuadExpr, a::AffExpr) = error("Cannot multiply a quadratic expression by an aff. expression")
(/)(q::QuadExpr, a::AffExpr) = error("Cannot divide a quadratic expression by an aff. expression")
# QuadExpr--QuadExpr
(+)(q1::QuadExpr, q2::QuadExpr) = QuadExpr(vcat(q1.qvars1,   q2.qvars1),
                                                                                     vcat(q1.qvars2,   q2.qvars2),
                                                                                     vcat(q1.qcoeffs,  q2.qcoeffs),
                                                                                     q1.aff + q2.aff)
(-)(q1::QuadExpr, q2::QuadExpr) = QuadExpr(vcat(q1.qvars1,   q2.qvars1),
                                                                                     vcat(q1.qvars2,   q2.qvars2),
                                                                                     vcat(q1.qcoeffs, -q2.qcoeffs),
                                                                                     q1.aff - q2.aff)
(*)(q1::QuadExpr, q2::QuadExpr) = error("Cannot multiply two quadratic expressions")
(/)(q1::QuadExpr, q2::QuadExpr) = error("Cannot divide a quadratic expression by a quadratic expression")
# QuadExpr--MatrixVar
(+)(lhs::QuadExpr, rhs::MatrixVar) = error("Cannot add a scalar quadratic expression and a matrix variable") #TODO: make this work w/ 1x1 matrices
(-)(lhs::QuadExpr, rhs::MatrixVar) = error("Cannot subtract a matrix variable from a scalar quadratic expression")
(*)(lhs::QuadExpr, rhs::MatrixVar) = error("Cannot multiply a scalar quadratic expression and a matrix variable")
(/)(lhs::QuadExpr, rhs::MatrixVar) = error("Cannot divide a scalar quadratic expression by a matrix variable")
# QuadExpr--MatrixExpr
(+)(lhs::QuadExpr, rhs::MatrixExpr) = error("Cannot add a scalar quadratic expression and a matrix expression") #TODO: make this work w/ 1x1 matrices
(-)(lhs::QuadExpr, rhs::MatrixExpr) = error("Cannot subtract a matrix expression from a scalar quadratic expression")
(*)(lhs::QuadExpr, rhs::MatrixExpr) = error("Cannot multiply a scalar quadratic expression and a matrix expression")
(/)(lhs::QuadExpr, rhs::MatrixExpr) = error("Cannot divide a scalar quadratic expression by a matrix expression")

# Matrix
# Matrix--MatrixVar
function (+)(lhs::Matrix, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot add matrices of unequal size")
    MatrixExpr([rhs], Array[eye(rhs)], lhs)
end
function (-)(lhs::Matrix, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot subtract matrices of unequal size")
    MatrixExpr([rhs], -Array[eye(rhs)], lhs)
end
function (*)(lhs::Matrix, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot multiply matrices of incompatible sizes")
    MatrixEpxr([rhs], Array[eye(lhs)], zero(rhs))
end
(/)(lhs::MatrixVar, rhs::Matrix) = error("Cannot divide matrices")
# MatrixVar
# MatrixVar--Number
(+)(lhs::MatrixVar, rhs::Number) = error("Cannot add a matrix variable and a scalar")
(-)(lhs::MatrixVar, rhs::Number) = error("Cannot subtract a matrix variable and a scalar")
(*)(lhs::MatrixVar, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::MatrixVar, rhs::Number) = (*)(1/rhs,lhs)
# MatrixVar--Matrix
function (+)(lhs::MatrixVar, rhs::Matrix)
    (size(lhs) == size(rhs)) || error("Cannot add matrices of unequal size")
    MatrixExpr([lhs], Array[eye(lhs)], rhs)
end
function (-)(lhs::MatrixVar, rhs::Matrix)
    (size(lhs) == size(rhs)) || error("Cannot subtract matrices of unequal size")
    MatrixExpr([lhs], Array[eye(lhs)], -rhs)
end
(*)(lhs::MatrixVar, rhs::Matrix) = error("Post-multiplying a matrix variable with a matrix is not yet supported")
(/)(lhs::MatrixVar, rhs::Matrix) = error("Cannot divide matrices")
# MatrixVar--Variable
(+)(lhs::MatrixVar, rhs::Variable) = error("Cannot add a matrix variable and a variable")
(-)(lhs::MatrixVar, rhs::Variable) = error("Cannot subtract a matrix variable and a variable")
(*)(lhs::MatrixVar, rhs::Variable) = error("Cannot multiply a matrix variable and a variable")
(/)(lhs::MatrixVar, rhs::Variable) = error("Cannot divide a matrix variable and a variable")
# MatrixVar--AffExpr
(+)(lhs::MatrixVar, rhs::AffExpr) = error("Cannot add a matrix variable and an affine expression")
(-)(lhs::MatrixVar, rhs::AffExpr) = error("Cannot subtract a matrix variable and an affine expression")
(*)(lhs::MatrixVar, rhs::AffExpr) = error("Cannot multiply a matrix variable and an affine expression")
(/)(lhs::MatrixVar, rhs::AffExpr) = error("Cannot divide a matrix variable and an affine expression")
# MatrixVar--QuadExpr
(+)(lhs::MatrixVar, rhs::QuadExpr) = error("Cannot add a matrix variable and a quadratic expression")
(-)(lhs::MatrixVar, rhs::QuadExpr) = error("Cannot subtract a matrix variable and a quadratic expression")
(*)(lhs::MatrixVar, rhs::QuadExpr) = error("Cannot multiply a matrix variable and a quadratic expression")
(/)(lhs::MatrixVar, rhs::QuadExpr) = error("Cannot divide a matrix variable and a quadratic expression")
# MatrixVar--MatrixVar
function (+)(lhs::MatrixVar, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot add matrix variables of unequal size")
    MatrixExpr([lhs,rhs],Array[eye(lhs),eye(rhs)], zero(lhs))
end
function (-)(lhs::MatrixVar, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot subtract matrix variables of unequal size")
    MatrixExpr([lhs,rhs],Array[eye(lhs),-eye(rhs)], zero(lhs))
end
(*)(lhs::MatrixVar, rhs::MatrixVar) = error("Cannot multiply matrix variables")
(/)(lhs::MatrixVar, rhs::MatrixVar) = error("Cannot divide matrix variables")
# MatrixVar--MatrixExpr
(+)(lhs::MatrixVar, rhs::MatrixExpr) = MatrixExpr(vcat(rhs.vars,lhs),Array[ rhs.coeffs...,eye(lhs)], rhs.constant)
(-)(lhs::MatrixVar, rhs::MatrixExpr) = MatrixExpr(vcat(rhs.vars,lhs),Array[-rhs.coeffs...,eye(lhs)],-rhs.constant)
(*)(lhs::MatrixVar, rhs::MatrixExpr) = error("Cannot multiply a matrix variable and a matrix expression")
(/)(lhs::MatrixVar, rhs::MatrixExpr) = error("Cannot divide a matrix variable by a matrix expression")

# MatrixExpr
# MatrixExpr--Number
(+)(lhs::MatrixExpr, rhs::Number) = error("Cannot add a matrix expression and a scalar")
(-)(lhs::MatrixExpr, rhs::Number) = error("Cannot subtract a matrix expression and a scalar")
(*)(lhs::MatrixExpr, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::MatrixExpr, rhs::Number) = (*)(1/rhs,lhs)
# MatrixExpr--Matrix
function (+)(lhs::MatrixExpr, rhs::Matrix)
    (size(lhs.constant) == size(rhs)) || error("Cannot add matrices of unequal size")
    MatrixExpr(copy(lhs.vars), copy(lhs.coeffs), lhs.constant+rhs)
end
function (-)(lhs::MatrixExpr, rhs::Matrix)
    (size(lhs.constant) == size(rhs)) || error("Cannot subtract matrices of unequal size")
    MatrixExpr(copy(lhs.vars), copy(lhs.coeffs), lhs.constant-rhs)
end
(*)(lhs::MatrixExpr, rhs::Matrix) = error("Post-multiplying matrices is not yet supported")
(/)(lhs::MatrixExpr, rhs::Matrix) = error("Cannot divide matrices")
# MatrixExpr--Variable
(+)(lhs::MatrixExpr, rhs::Variable) = error("Cannot add a matrix expression and a variable")
(-)(lhs::MatrixExpr, rhs::Variable) = error("Cannot subtract a matrix expression and a variable")
(*)(lhs::MatrixExpr, rhs::Variable) = error("Cannot multiply a matrix expression and a variable")
(/)(lhs::MatrixExpr, rhs::Variable) = error("Cannot divide a matrix expression and a variable")
# MatrixExpr--AffExpr
(+)(lhs::MatrixExpr, rhs::AffExpr) = error("Cannot add a matrix expression and an affine expression")
(-)(lhs::MatrixExpr, rhs::AffExpr) = error("Cannot subtract a matrix expression and an affine expression")
(*)(lhs::MatrixExpr, rhs::AffExpr) = error("Cannot multiply a matrix expression and an affine expression")
(/)(lhs::MatrixExpr, rhs::AffExpr) = error("Cannot divide a matrix expression and an affine expression")
# MatrixExpr--QuadExpr
(+)(lhs::MatrixExpr, rhs::QuadExpr) = error("Cannot add a matrix expression and a quadratic expression")
(-)(lhs::MatrixExpr, rhs::QuadExpr) = error("Cannot subtract a matrix expression and a quadratic expression")
(*)(lhs::MatrixExpr, rhs::QuadExpr) = error("Cannot multiply a matrix expression and a quadratic expression")
(/)(lhs::MatrixExpr, rhs::QuadExpr) = error("Cannot divide a matrix expression and a quadratic expression")
# MatrixExpr--MatrixVar
function (+)(lhs::MatrixExpr, rhs::MatrixVar)
    (size(lhs.constant) == size(rhs)) || error("Cannot add matrix variables of unequal size")
    MatrixExpr(vcat(lhs.vars,rhs), Array[lhs.coeffs...,eye(rhs)], lhs.constant)
end
function (-)(lhs::MatrixExpr, rhs::MatrixVar)
    (size(lhs.constant) == size(rhs)) || error("Cannot subtract matrix variables of unequal size")
    MatrixExpr(vcat(lhs.vars,rhs), Array[lhs.coeffs...,-eye(rhs)], lhs.constant)
end
(*)(lhs::MatrixExpr, rhs::MatrixVar) = error("Cannot multiply matrix variables")
(/)(lhs::MatrixExpr, rhs::MatrixVar) = error("Cannot divide matrix variables")
# MatrixExpr--MatrixExpr
(+)(lhs::MatrixExpr, rhs::MatrixExpr) = MatrixExpr(vcat(lhs.vars,rhs.vars),Array[lhs.coeffs..., rhs.coeffs...],lhs.constant+rhs.constant)
(-)(lhs::MatrixExpr, rhs::MatrixExpr) = MatrixExpr(vcat(lhs.vars,rhs.vars),Array[lhs.coeffs...,-rhs.coeffs...],lhs.constant-rhs.constant)
(*)(lhs::MatrixExpr, rhs::MatrixExpr) = error("Cannot multiply a matrix variable and a matrix expression")
(/)(lhs::MatrixExpr, rhs::MatrixExpr) = error("Cannot divide a matrix variable by a matrix expression")

# LinearConstraint
# Number--???
(<=)(lhs::Number, rhs::Variable) = (>=)(rhs, lhs)
(==)(lhs::Number, rhs::Variable) = (==)(rhs, lhs)
(>=)(lhs::Number, rhs::Variable) = (<=)(rhs, lhs)

(<=)(lhs::Number, rhs::AffExpr) = (>=)(rhs, lhs)
(==)(lhs::Number, rhs::AffExpr) = (==)(rhs, lhs)
(>=)(lhs::Number, rhs::AffExpr) = (<=)(rhs, lhs)
(<=)(lhs::Number, rhs::QuadExpr) = (>=)(rhs, lhs)
(==)(lhs::Number, rhs::QuadExpr) = (==)(rhs, lhs)
(>=)(lhs::Number, rhs::QuadExpr) = (<=)(rhs, lhs)
# Variable--???
(<=)(lhs::Variable, rhs::Number) = (<=)(lhs - rhs, 0.0)
(==)(lhs::Variable, rhs::Number) = (==)(lhs - rhs, 0.0)
(>=)(lhs::Variable, rhs::Number) = (>=)(lhs - rhs, 0.0)

(<=)(lhs::Variable, rhs::Variable) = (<=)(lhs - rhs, 0.0)
(==)(lhs::Variable, rhs::Variable) = (==)(lhs - rhs, 0.0)
(>=)(lhs::Variable, rhs::Variable) = (>=)(lhs - rhs, 0.0)

(<=)(lhs::Variable, rhs::AffExpr) = (<=)(lhs - rhs, 0.0)
(==)(lhs::Variable, rhs::AffExpr) = (==)(lhs - rhs, 0.0)
(>=)(lhs::Variable, rhs::AffExpr) = (>=)(lhs - rhs, 0.0)
(<=)(lhs::Variable, rhs::QuadExpr) = (<=)(lhs - rhs, 0.0)
(==)(lhs::Variable, rhs::QuadExpr) = (==)(lhs - rhs, 0.0)
(>=)(lhs::Variable, rhs::QuadExpr) = (>=)(lhs - rhs, 0.0)
# AffExpr--???
(<=)(lhs::AffExpr, rhs::Number) = LinearConstraint(lhs,            -Inf,rhs-lhs.constant)
(==)(lhs::AffExpr, rhs::Number) = LinearConstraint(lhs,rhs-lhs.constant,rhs-lhs.constant)
(>=)(lhs::AffExpr, rhs::Number) = LinearConstraint(lhs,rhs-lhs.constant,             Inf)
(<=)(lhs::AffExpr, rhs::Variable) = (<=)(lhs-rhs, 0.0)
(==)(lhs::AffExpr, rhs::Variable) = (==)(lhs-rhs, 0.0)
(>=)(lhs::AffExpr, rhs::Variable) = (>=)(lhs-rhs, 0.0)
(<=)(lhs::AffExpr, rhs::AffExpr) = (<=)(lhs-rhs, 0.0)
(==)(lhs::AffExpr, rhs::AffExpr) = (==)(lhs-rhs, 0.0)
(>=)(lhs::AffExpr, rhs::AffExpr) = (>=)(lhs-rhs, 0.0)
(<=) (lhs::AffExpr, rhs::QuadExpr) = (<=)(lhs-rhs, 0)
(==) (lhs::AffExpr, rhs::QuadExpr) = (==)(lhs-rhs, 0)
(>=) (lhs::AffExpr, rhs::QuadExpr) = (>=)(lhs-rhs, 0)

# There's no easy way to allow operator overloads for range constraints.
# Use macros instead.

# QuadConstraint
# QuadConstraint--Number
(<=) (lhs::QuadExpr, rhs::Number)   = QuadConstraint( QuadExpr(copy(lhs.qvars1), copy(lhs.qvars2), lhs.qcoeffs,lhs.aff - rhs), :<=   )
(==) (lhs::QuadExpr, rhs::Number)   = QuadConstraint( QuadExpr(copy(lhs.qvars1), copy(lhs.qvars2), lhs.qcoeffs,lhs.aff - rhs), :(==) )
(>=) (lhs::QuadExpr, rhs::Number)   = QuadConstraint( QuadExpr(copy(lhs.qvars1), copy(lhs.qvars2), lhs.qcoeffs,lhs.aff - rhs), :>=   )
(<=) (lhs::QuadExpr, rhs::Variable) = (<=)(lhs-rhs, 0)
(==) (lhs::QuadExpr, rhs::Variable) = (==)(lhs-rhs, 0)
(>=) (lhs::QuadExpr, rhs::Variable) = (>=)(lhs-rhs, 0)
(<=) (lhs::QuadExpr, rhs::AffExpr)  = (<=)(lhs-rhs, 0)
(==) (lhs::QuadExpr, rhs::AffExpr)  = (==)(lhs-rhs, 0)
(>=) (lhs::QuadExpr, rhs::AffExpr)  = (>=)(lhs-rhs, 0)
(<=) (lhs::QuadExpr, rhs::QuadExpr) = (<=)(lhs-rhs, 0)
(==) (lhs::QuadExpr, rhs::QuadExpr) = (==)(lhs-rhs, 0)
(>=) (lhs::QuadExpr, rhs::QuadExpr) = (>=)(lhs-rhs, 0)

# Matrix
# Matrix--MatrixVar
function (<=)(lhs::Matrix, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr([rhs], -Array[eye(rhs)], lhs), :<=)
end
function (==)(lhs::Matrix, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr([rhs], -Array[eye(rhs)], lhs), :(==))
end
function (>=)(lhs::Matrix, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr([rhs], -Array[eye(rhs)], lhs), :>=)
end
# Matrix--MatrixExpr
function (<=)(lhs::Matrix, rhs::MatrixExpr)
    (size(lhs) == size(rhs.constant)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(copy(rhs.vars), -rhs.coeffs, lhs-rhs.constant), :<=)
end
function (==)(lhs::Matrix, rhs::MatrixExpr)
    (size(lhs) == size(rhs.constant)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(copy(rhs.vars), -rhs.coeffs, lhs-rhs.constant), :(==))
end
function (>=)(lhs::Matrix, rhs::MatrixExpr)
    (size(lhs) == size(rhs.constant)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(copy(rhs.vars), -rhs.coeffs, lhs-rhs.constant), :>=)
end

# MatrixVar
# MatrixVar--Matrix
function (<=)(lhs::MatrixVar, rhs::Matrix)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr([lhs], Array[eye(lhs)], -rhs), :<=)
end
function (==)(lhs::MatrixVar, rhs::Matrix)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr([lhs], Array[eye(lhs)], -rhs), :(==))
end
function (>=)(lhs::MatrixVar, rhs::Matrix)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr([lhs], Array[eye(lhs)], -rhs), :>=)
end
# MatrixVar--MatrixVar
function (<=)(lhs::MatrixVar, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr([lhs,rhs], vcat(eye(lhs),-eye(rhs)), zero(lhs)), :<=)
end
function (==)(lhs::MatrixVar, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr([lhs,rhs], vcat(eye(lhs),-eye(rhs)), zero(lhs)), :(==))
end
function (>=)(lhs::MatrixVar, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr([lhs,rhs], vcat(eye(lhs),-eye(rhs)), zero(lhs)), :>=)
end
# MatrixVar--MatrixExpr
function (<=)(lhs::MatrixVar, rhs::MatrixExpr)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(vcat(rhs.vars,lhs), Array[-rhs.coeffs...,eye(lhs)], -rhs.constant), :<=)
end
function (==)(lhs::MatrixVar, rhs::MatrixExpr)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(vcat(rhs.vars,lhs), Array[-rhs.coeffs...,eye(lhs)], -rhs.constant), :(==))
end
function (>=)(lhs::MatrixVar, rhs::MatrixExpr)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(vcat(rhs.vars,lhs), Array[-rhs.coeffs...,eye(lhs)], -rhs.constant), :>=)
end

# MatrixExpr
# MatrixExpr--Matrix
function (<=)(lhs::MatrixExpr, rhs::Matrix)
    (size(lhs.constant) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(copy(lhs.vars), copy(lhs.coeffs), lhs.constant-rhs), :<=)
end
function (==)(lhs::MatrixExpr, rhs::Matrix)
    (size(lhs.constant) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(copy(lhs.vars), copy(lhs.coeffs), lhs.constant-rhs), :(==))
end
function (>=)(lhs::MatrixExpr, rhs::Matrix)
    (size(lhs.constant) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(copy(lhs.vars), copy(lhs.coeffs), lhs.constant-rhs), :>=)
end
# MatrixExpr--MatrixVar
function (<=)(lhs::MatrixExpr, rhs::MatrixVar)
    (size(lhs.constant) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(copy(lhs.vars), copy(lhs.coeffs), lhs.constant-rhs), :<=)
end
function (==)(lhs::MatrixExpr, rhs::MatrixVar)
    (size(lhs.constant) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(copy(lhs.vars), copy(lhs.coeffs), lhs.constant-rhs), :(==))
end
function (>=)(lhs::MatrixExpr, rhs::MatrixVar)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(copy(lhs.vars), copy(lhs.coeffs), lhs.constant-rhs), :>=)
end
# MatrixExpr--MatrixExpr
function (<=)(lhs::MatrixExpr, rhs::MatrixExpr)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(vcat(lhs.vars,rhs.vars), Array[lhs.coeffs...,rhs.coeffs...], lhs.constant-rhs.constant), :<=)
end
function (==)(lhs::MatrixExpr, rhs::MatrixExpr)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(vcat(lhs.vars,rhs.vars), Array[lhs.coeffs...,rhs.coeffs...], lhs.constant-rhs.constant), :(==))
end
function (>=)(lhs::MatrixExpr, rhs::MatrixExpr)
    (size(lhs) == size(rhs)) || error("Cannot compare matrices of different sizes")
    MatrixConstraint(MatrixExpr(vcat(lhs.vars,rhs.vars), Array[lhs.coeffs...,rhs.coeffs...], lhs.constant-rhs.constant), :>=)
end
