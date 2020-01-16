
# Epigenetic inactivation of the autophagy–lysosomal system in the Parkinson’s disease appendix

- [Source code](http://git.vugene.eu/gordevicius/pdappendixalp)

- [Website of data analysis results](http://www.vugene.eu/VAI/pdappendixalp/)


# Running the analysis on local machine

- Install conda 
- Install docker
- [Download the tarball of code and data](www/archive_20200115.tar.gz)

Create execution environment:

<pre><code>
make importenv
conda activate ALP
</code></pre>

Repeat all analyses by running R and Rmd scripts:

<pre><code>
make all
</code></pre>

Build a shiny-server docker container:

<pre><code>
make -C shinyserver all
</code></pre>

Run the results website:

<pre><code>
make run
</code></pre>

The website will be available at [http://localhost:3838](http://localhost:3838).

# Contact

Send a message to juozas.gordevicius at vai.org