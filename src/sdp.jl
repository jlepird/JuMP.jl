##########################
# Matrix variable type
abstract MatrixVar

type MatrixVarExplicit <: MatrixVar
    m::Model
    innerArray::Array{Variable,2}
end

type MatrixVarImplicit <: MatrixVar
    m::Model
    dim::Int
end

transpose(a::MatrixVar)  = a
ctranspose(a::MatrixVar) = a

# need to figure out how to get variable name here...
function show(io::IO,d::MatrixVar) 
    arr = Array(String, size(d))
    for i in 1:size(d)[1]
        for j in 1:size(d)[2]
            arr[i,j] = "X[$i,$j]"
        end
    end
    show(io, arr)
end

getindex(d::MatrixVarExplicit) = getindex(d.innerArray)
setindex!(d::MatrixVarExplicit,val) = setindex!(d.innerArray,val)

# Use this to work with elements of a MatrixVar maybe?
# Maybe extend to ranges so you can add constraints on blocks?
type MatrixVarRef
    var::MatrixVarImplicit
    x::Int
    y::Int
end

getindex(d::MatrixVarImplicit, x::Int, y::Int) = MatrixVarRef(d,x,y)

size(d::MatrixVarExplicit) = size(d.innerArray)
eye(d::MatrixVarExplicit) = eye(d.innerArray)
zero(d::MatrixVarExplicit) = zero(d.innerArray)

size(d::MatrixVarImplicit) = (d.dim,d.dim)
eye(d::MatrixVarImplicit) = eye(d.dim)
zero(d::MatrixVarImplicit) = zeros(d.dim,d.dim)

type MatrixSparseExpr
    vars::Vector{MatrixVar}
    pre::Vector{SparseMatrixCSC{Float64}}
    post::Vector{SparseMatrixCSC{Float64}}
    constant::SparseMatrixCSC{Float64}
end

type MatrixExpr
    vars::Vector{MatrixVar}
    pre::Vector{Matrix{Float64}}
    post::Vector{Matrix{Float64}}
    constant::Matrix{Float64}
end # TODO: check that all matrices are symmetric!

MatrixExpr(n::Int) = MatrixExpr(MatrixVar[],Matrix{Float64}[],Matrix{Float64}[],zeros(n,n))

transpose(d::MatrixExpr)  = MatrixExpr(copy(d.vars), copy(d.post), copy(d.pre), d.constant)
ctranspose(d::MatrixExpr) = MatrixExpr(copy(d.vars), copy(d.post), copy(d.pre), d.constant)

type NestedMatrixExpr
    arr::Matrix{NestedMatrixExpr}
end

type MatrixConstraint <: JuMPConstraint
    terms::MatrixExpr
    sense::Symbol
end

type MatrixFunctionExpr
    vars
    constant::Real
end

type MatrixFunctionConstraint <: JuMPConstraint
    expr::MatrixFunctionExpr
    sense::Symbol
end

# some operation on MatrixExpr (e.g. norm, trace, etc.). Must return a scalar.
abstract MatrixFunction

type MatrixTrace <: MatrixFunction
    expr::MatrixExpr
end

type MatrixNorm <: MatrixFunction
    expr::MatrixExpr
end

trace(c::MatrixExpr) = MatrixTrace(c)
norm(c::MatrixExpr) = MatrixNorm(c)

function addConstraint(m::Model, c::MatrixConstraint)
    push!(m.matrixconstr,c)
    return ConstraintRef{MatrixConstraint}(m,length(m.matrixconstr))
end

# function hcat(args::)

# end

# function vcat()

# end

# function hvcat()

# end

## operators we need:
## * trace
## * norm
## * concatenation (e.g. [AX I; I X])
## Need to handle linear constraints on elements of MatrixVar (see SDPA docs 9.1)
