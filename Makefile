LDir=libs
SDir=src


one_shot_install: 
	@make -s build 
	@make -s install 
	@make -s clean

build: ${SDir}/*

	@set -e
	
	@if test ! -d ${LDir}; then \
		mkdir -p ${LDir}; \
	fi

	@echo "\nBuilding..."

	@for file in ${PWD}/${SDir}/* ; do \
		filename="$${file##*/}" ; \
		echo "Compiling $${file} file..." ; \
		g++ -c $${file} -I${PWD} ; \
		echo "Building $${PWD}/${LDir}/$${filename%.*}.o file..." ; \
		mv $${filename%.*}.o ${LDir}/ ; \
		echo ; \
	done

	@echo "Successfully done.\n"


install: ${LDir}/*

	@set -e

	@echo "\nInstalling..."
	@sudo ar rsv libmaths.so ${PWD}/${LDir}/*
	@sudo mv libmaths.so /usr/local/lib/
	@sudo cp -r ${PWD}/maths /usr/local/include/
	@echo "Done.\n"


uninstall: /usr/local/lib/libmaths.so /usr/local/include/maths

	@echo "\nUninstalling..."
	@sudo rm -rf /usr/local/include/maths
	@sudo rm -f /usr/local/lib/libmaths.so
	@echo "Done.\n"


clean: ${LDir}/

	@if test -d ${LDir}; then \
		echo "\nCleaning..." ; \
		rm -rf ${LDir} ; \
		echo "Successfully done.\n" ; \
	fi

	