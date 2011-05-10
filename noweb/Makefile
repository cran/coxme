PARTS = main.Rnw \
	basic.Rnw \
	build.Rnw \
	formula.Rnw \
	varfun.Rnw  varfun2.Rnw \
	fit.Rnw \
	ranef.Rnw \
#	reconcile.Rnw

SFUN =  coxme.R \
	coxme.fit.R\
	expand.nested.R \
	formula.R \
	coxmeFull.R \
	coxmeMlist.R \
	ranef.R
#	matrix.reconcile.R

RDIR = ../R

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
$(SFUN) : all.nw

%.R: 
	echo "# Automatically generated from all.nw using noweb" > $@
	echo "source('noweb.R'); print(notangle(nwread('all.nw'), target='$(*F)'), file='zz')" |  $(R_HOME)/bin$(R_ARCH_BIN)/R --vanilla --slave
	cat zz >> $@
	rm zz
#	notangle -R$(*F) all.nw >> $@


clean:
	-rm all.nw all.log all.aux all.toc all.tex all.pdf
	-rm $(RFUN)

noweb.R: noweb.Rnw
	notangle noweb.Rnw > noweb.R
