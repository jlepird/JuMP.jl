#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# Hock-Schittkowski Nonlinear Test Suite
# These files are JuMP implementations of the model described in 
#  W. Hock, K. Schittkowski, Test Examples for Nonlinear Programming
#  Codes, Lecture Notes in Economics and Mathematical Systems, 
#  Springer, No, 187, 1981 
# More information, including original model description, at
# http://www.ai7.uni-bayreuth.de/downloads.htm
#
# This file will run all the HS test problems, allowing for easy testing
#############################################################################

if Pkg.installed("Ipopt") == nothing
    error("Ipopt must be installed to run the NLP tests")
end

HS_path = joinpath(Pkg.dir("JuMP"),"test","hockschittkowski")

files = readdir(HS_path)

println("Running HS Tests")
for file in files
    contains(file, "runhs") && continue
    try
        #readall(`julia `)
        readall(`julia $(joinpath(HS_path,file))`)
        println("  $file  PASSED")
    catch
        println("  $file  FAILED")
    end
end