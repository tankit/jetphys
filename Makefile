
clean:
	rm -f ./*_C*
	rm -f ./JetMETCorrections/Modules/src/*_cc*
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
	rm -f pypdf/*.pdf
	rm -f tex/*.toc
	rm -f tex/*.doc
	rm -f tex/*.log
	rm -f tex/*.nav
	rm -f tex/*.aux
	rm -f tex/*.out

