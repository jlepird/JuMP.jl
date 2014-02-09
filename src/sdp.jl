##########################
# Matrix variable type
abstract MatrixVar

type MatrixVarExplicit <: MatrixVar
  innerArray::Array{Variable,2}
  name::String
end

type MatrixVarImplicit <: MatrixVar
  dim::Int
  name::String
end

getindex(d::MatrixVar) = getindex(d.innerArray)
setindex!(d::MatrixVar,val) = setindex!(d.innerArray,val)

size(d::MatrixVar) = size(d.innerArray)
one(d::MatrixVar) = one(d.innerArray)
zero(d::MatrixVar) = zero(d.innerArray)

typealias MatrixExpr GenericAffExpr{Matrix{Float64},MatrixVar}
MatrixExpr(n::Int) = MatrixExpr(MatrixVar[],Matrix{Float64}[],zeros(n,n))

type MatrixConstraint <: JuMPConstraint
    terms::MatrixExpr
    sense::Symbol
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
