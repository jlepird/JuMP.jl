##########################
# Matrix variable type
type MatrixVar
  innerArray::Array{Variable,2}
  name::String
end

getindex(d::MatrixVar) = getindex(d.innerArray)
setindex!(d::MatrixVar,val) = setindex!(d.innerArray,val)

typealias MatrixExpr GenericAffExpr{Matrix{Float64},MatrixVar}
MatrixExpr(n::Int) = MatrixExpr(MatrixVar[],Matrix{Float64}[],zeros(n,n))

type MatrixConstraint <: JuMPConstraint
    terms::MatrixExpr
    lb::Vector{Float64}
end

function addConstraint(m::Model, c::MatrixConstraint)
    push!(m.matrixconstr,c)
    return ConstraintRef{MatrixConstraint}(m,length(m.matrixconstr))
end

## operators we need:
## * inverse
## * trace
## * transpose
## * concatenation (e.g. [AX I; I X])
