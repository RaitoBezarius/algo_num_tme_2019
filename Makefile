source := src_rapports
output := rendus
images := images

sources := $(wildcard $(source)/*.md)
objects := $(patsubst %.md,%.pdf,$(subst $(source),$(output),$(sources)))

all: $(objects)

$(output)/%.pdf: $(source)/%.md $(source)/%.bib
		pandoc $< \
			--filter pandoc-citeproc \
			--resource-path=$(images) \
			--pdf-engine=pdflatex \
			--bibliography $(word 2,$^) \
			--csl=$(source)/chicago-author-date.csl \
			-o $@


$(output)/%.pdf: $(source)/%.md
		pandoc $< \
			--resource-path=$(images) \
			--pdf-engine=pdflatex \
			-o $@


.PHONY: clean

clean:
	rm -f $(output)/*.pdf
