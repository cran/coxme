.PHONY: fun clean doc all

PARTS = main.Rnw \
	basic.Rnw \
	build.Rnw \
	formula.Rnw \
	varfun.Rnw  varfun2.Rnw \
	fit.Rnw \
	ranef.Rnw \
	lmekin.Rnw \
	bdsmatrix.Rnw
#	reconcile.Rnw

SFUN =  coxme.R \
	coxme.fit.R\
	expand.nested.R \
	formula.R \
	coxmeFull.R \
	coxmeMlist.R \
	lmekin.R\
	ranef.R\
	bdsmatrix.R \
#	matrix.reconcile.R

CFUN =  bds_dsc.c

RDIR = ../R
all: code.pdf fun

code.pdf: code.tex noweb.sty
	pdflatex code.tex

code.nw: $(PARTS)
	cat $(PARTS) > code.nw
	echo "\\end{document}" >> code.nw

code.tex: code.nw
	echo "require(noweb); noweave('code.nw', target='$(*F)', out='$@')" | R --slave


RFUN = $(SFUN:%=$(RDIR)/%)
CFUN2= $(CFUN:%=../src/%)

fun: $(RFUN) $(CFUN2)

$(RFUN) : code.nw 
$(SFUN) : code.nw
$(CFUN2): code.nw

%.R: 
	echo "# Automatically generated from the noweb directory" > $@
	echo "require(noweb); notangle('code.nw', target='$(*F)', out='zz')" | R --slave
	cat zz >> $@
	rm zz
#	notangle -R$(*F) code.nw >> $@

%.c: 
	echo "/* Automatically generated from the noweb directory */" > $@
	echo "require(noweb); notangle('code.nw', target='$(*F)', out='zz')" | R --slave
	cat zz >> $@
	rm zz

clean:
	-rm code.nw code.log code.aux code.toc code.tex code.pdf code.out
	-rm noweb.sty

noweb.sty:
	echo 'library(noweb); data(noweb); cat(noweb.sty, sep="\n", file="noweb.sty")' |  R --vanilla --slave
