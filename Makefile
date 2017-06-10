install:
	@echo 'roxygen2::roxygenise()' | R --vanilla --silent
	@cd docs; make; cd -
	@R --slave -e 'Rcpp::compileAttributes()'
	@R CMD build ./ --no-manual
	@((R CMD INSTALL mashr_*.tar.gz -l $(shell echo "cat(.libPaths()[1])" | R --slave) && rm -rf tmp.* mashr_*.tar.gz) || ($(ECHO) "Installation failed"))
