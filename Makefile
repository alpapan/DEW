SHELL := /bin/bash

all:
	if [ ! -d 3rd_party/bin ]; then mkdir 3rd_party/bin; fi
	cd 3rd_party/bin && if [ ! -f biokanga ]; then ln ../biokanga; fi && if [ ! -f express ]; then ln ../express .;fi
	cd 3rd_party/justpreprocessmyreads && $(MAKE) 3rd_party
	cd 3rd_party/samtools && $(MAKE) && cp -f samtools ../bin/ && $(MAKE) clean
	cd 3rd_party/bwa && $(MAKE) && cp -f bwa qualfa2fq.pl xa2multi.pl ../bin/ && $(MAKE) clean
	cd 3rd_party && cp -f bedtools bwa bin/
	chmod a+rx 3rd_party/bin/*
	echo "Installation complete."
clean:
	cd 3rd_party/samtools && $(MAKE) clean
	cd 3rd_party/bedtools && $(MAKE) clean
	cd 3rd_party/bin && rm -fr *


###################################################################


