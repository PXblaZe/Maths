LDir=libs
SDir=src


build: src/*

	@set -e
	
	@if test ! -d ${LDir}; then \
		mkdir -p ${LDir}; \
	fi

	@for file in ${PWD}/${SDir}/* ; do \
		filename="$${file##*/}" ; \
		echo "Compiling $${file} file..." ; \
		g++ -c $${file} -I${PWD} ; \
		echo "Building $${PWD}/${LDir}/$${filename%.*}.o file..." ; \
		mv $${filename%.*}.o ${LDir}/ ; \
		echo ; \
	done

	@echo "Successfully done."

# install: build

clean: ${LDir}/

	@echo "Cleaning..."
	@rm -rf ${LDir}
	@echo "Successfully done."

	