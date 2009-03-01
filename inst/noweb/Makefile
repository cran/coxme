PARTS = main.Rnw \
	basic.Rnw \
	build.Rnw \
	formula.Rnw \
	varfun.Rnw  varfun2.Rnw \
	fit.Rnw 
#	reconcile.Rnw

SFUN =  coxme.R \
	coxme.fit.R\
	expand.nested.R \
	formula.R \
	coxmeFull.R \
	coxmeMlist.R 
#	matrix.reconcile.R

RDIR = ../../R

all.pdf: all.tex
	pdflatex all

all.nw: $(PARTS)
	cat $(PARTS) > all.nw
	echo "\\end{document}" >> all.nw

all.tex: all.nw
	noweave -delay all.nw > all.tex

RFUN = $(SFUN:%=$(RDIR)/%)

fun:  $(RFUN)

$(RFUN) : all.nw

%.R: 
	echo "# Automatically generated from all.nw using noweb" > $@
	notangle -R$(*F) all.nw >> $@

clean:
	-rm all.nw all.log all.aux all.toc all.tex all.pdf
	-rm $(RFUN)