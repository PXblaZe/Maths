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
	@sudo mv libmaths.so /usr/lib/
	@sudo cp -r ${PWD}/maths /usr/include/
	@echo "Done."


uninstall: /usr/lib/libmaths.so /usr/include/maths

	@echo "Uninstalling..."
	@sudo rm -rf /usr/include/maths
	@sudo rm -f /usr/lib/libmaths.so
	@echo "Done."


clean: ${LDir}/

	@if test -d ${LDir}; then \
		echo "Cleaning..." ; \
		rm -rf ${LDir} ; \
		echo "Successfully done." ; \
	fi

	