#=
Generates a bunch of random matrices, any of which can be used to induce
boundary roughness in our model
=#

using FileIO

# Random matrix parameters
M = 16               #number of sine modes. complexity of edge roughness
matrix_size = 100     #size of the matrix used for the random coefficients

for matrixnumber in 1:3

  matrixname = "Uniform_M=$(M)_size=$(matrix_size)_$(matrixnumber).jld2"
  matrix = zeros(Float64, matrix_size, matrix_size, M)

  for k in 1:M, j in 1:matrix_size, i in 1:matrix_size
    matrix[i,j,k] = rand() - 0.5
  end
  save(matrixname, "matrix", matrix) # write file to disk

end
