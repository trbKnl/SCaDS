#####################################################
# Create a package for the sparce SCA cpp function
#####################################################

require(Rcpp)

Rcpp.package.skeleton(name = "sparseSCAcppFunction", list = character(), 
	environment = .GlobalEnv, path = ".", force = FALSE, 
	code_files = character(), cpp_files = "./sparseSCA.cpp",
	example_code = TRUE, attributes = TRUE, module = FALSE, 
	author = "Niek de Schipper", 
	maintainer = "Niek de Schipper", 
	email = "your@email.com", 
	license = "GPL (>= 2)"
	)



