HTML_FILES := $(patsubst %.Rmd, %.html ,$(wildcard *.Rmd)) \
              $(patsubst %.md, %.html ,$(wildcard *.md))

all: html

html: $(HTML_FILES)

devel:
	touch restart.txt
	rsync -rav \
	. minge.ibt.lt:ShinyApps/VAI/PD2018/IndexPage

%.html: %.Rmd
	R --slave -e "set.seed(100);rmarkdown::render('$<', encoding = 'UTF-8')"
