LDir=libs
SDir=src


build: ${SDir}/*

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


install: ${LDir}/*

	@set -e

	@echo "Installing..."
	@ar rsv libmaths.so ${PWD}/${LDir}/*
	@sudo mv libmaths.so /usr/local/lib/
	@sudo cp -r ${PWD}/maths /usr/local/include/
	@echo "Done."


uninstall: /usr/local/lib/libmaths.so /usr/local/include/maths

	@echo "Uninstalling..."
	@sudo rm -rf /usr/local/include/maths
	@sudo rm -f /usr/local/lib/libmaths.so
	@echo "Done."


clean: ${LDir}/

	@if test -d ${LDir}; then \
		echo "Cleaning..." ; \
		rm -rf ${LDir} ; \
		echo "Successfully done." ; \
	fi

	