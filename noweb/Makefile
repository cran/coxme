.PHONY: fun clean doc

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
all: all.pdf fun

all.pdf: all.tex
	pdflatex all

all.nw: $(PARTS)
	cat $(PARTS) > all.nw
	echo "\\end{document}" >> all.nw

all.tex: all.nw
	noweave -delay all.nw > all.tex

RFUN = $(SFUN:%=$(RDIR)/%)
CFUN2= $(CFUN:%=../src/%)

fun: $(RFUN) $(CFUN2)

$(RFUN) : all.nw 
$(SFUN) : all.nw
$(CFUN2): all.nw

%.R: 
	echo "# Automatically generated from all.nw using noweb" > $@
	echo "source('noweb.R'); print(notangle(nwread('all.nw'), target='$(*F)'), file='zz')" |  $(R_HOME)/bin$(R_ARCH_BIN)/R --vanilla --slave
	cat zz >> $@
	rm zz
#	notangle -R$(*F) all.nw >> $@

%.c: 
	echo "/* Automatically generated from all.nw using noweb */" > $@
	echo "source('noweb.R'); print(notangle(nwread('all.nw'), target='$(*F)'), file='zz')" |  $(R_HOME)/bin$(R_ARCH_BIN)/R  --vanilla --slave
	cat zz >> $@
	rm zz

clean:
	-rm all.nw all.log all.aux all.toc all.tex all.pdf
	-rm $(RFUN)

noweb.R: noweb.Rnw
	notangle noweb.Rnw > noweb.R
