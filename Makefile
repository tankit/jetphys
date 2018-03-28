
clean:
	rm -f ./*_C*
	rm -f ./CondFormats/JetMETObjects/src/*_cc*

distclean:
	rm -f ./*_C*
	rm -f ./CondFormats/JetMETObjects/src/*_cc*
	rm -f reports/*.log
	rm -f reports/*.json
	rm -f ./*.txt
	rm -rf doc/latex/*
	rm -rf doc/html/*
	rm -f ./*.root
