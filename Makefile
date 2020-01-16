objects=\
	IndexPage/index.html \
	Appendix_PDvsControls_Padlock/index.html \
  Mice_CecalPatch_Padlock/index.html \
  Mice_DSS_Padlock/index.html \
  Appendix_PDvsControls_Padlock/index.html \
  Appendix_PDvsControls_RNAseq/index.html \
  Appendix_AgeAcceleration_Padlock/index.html \
  Appendix_AgeAcceleration_RNAseq/index.html \
  Brain_PFC_Padlock_CGonly/index.html \
  Brain_OFB_Padlock_CGonly/index.html \
  Brain_PFCRep_Padlock_withGLU/index.html \
  Brain_PFCRep_Padlock_withGLU_Braak/index.html \
  Brain_AgeAcceleration_Padlock/index.html \
  Discover_Pathways/index.html \
  Results_Appendix/index.html \
  Figure1/index.html \
  Results_Brain/index.html \
  Figure2/index.html \
  Results_Proteomics/index.html \
  Figure3/index.html \
  Results_Aging/index.html \
  Figure4/index.html \
  Results_Mice/index.html \
  Figure5/index.html \
  Results_etc/index.html

all: navbar $(objects)


navbar: include/before_body.html

# Prepare navigation bar
include/before_body.html: code/generateNavigationBar.R index.json
	Rscript code/generateNavigationBar.R

# Rendering a Rmd file
# has to be redone if navbar has changed
%.html: %.Rmd include/before_body.html
	Rscript -e 'rmarkdown::render("$<", knit_root_dir="./")'
	touch $(dir $<)restart.txt

runAmiGO:
	# run by hand
	mkdir -p AmiGO/index
	cd AmiGO/index
	curl -O http://release.geneontology.org/2019-07-01/products/solr/golr-index-contents.tgz
	tar -zxvf golr-index-contents.tgz
	cd ..
	docker run -p 8080:8080 -p 9999:9999 -v ${PWD}:/srv/solr/data -t geneontology/amigo-standalone

run:
	mkdir -p tmp/log
	docker run -p 3838:3838 \
	-v ${PWD}:/srv/shiny-server/myapp:ro \
	-v ${PWD}/tmp/log:/var/log/shiny-server:Z \
	vugene/shinyserver


exportenv:
	conda env export --from-history | grep -v "^prefix: " > env.yml

importenv:
	conda env create --force --file env.yml