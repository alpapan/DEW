SHELL := /bin/bash

all: 
	if [ ! -d 3rd_party/bin ]; then mkdir 3rd_party/bin; fi
	cd 3rd_party/bin && ln ../biokanga ../express .
	cd 3rd_party/justpreprocessmyreads && $(MAKE) 3rd_party
	cd 3rd_party/samtools && $(MAKE) && cp samtools ../bin/ && $(MAKE) clean
	cd 3rd_party/bedtools && if [ ! -d bin ]; then mkdir bin; fi && $(MAKE) && cp bin/* ../bin/ && $(MAKE) clean
	chmod a+rx 3rd_party/bin/*
	echo "Installation complete."
clean:
	cd 3rd_party/samtools && $(MAKE) clean
	cd 3rd_party/bedtools && $(MAKE) clean
	cd 3rd_party/bin && rm -fr *


###################################################################


