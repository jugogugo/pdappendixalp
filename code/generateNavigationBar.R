# generate before_body.html
require(glue)

local({
  before_body <- "include/before_body.html"


  index <- jsonlite::read_json("index.json")  



  cat('
<div class="navbar navbar-default navbar-fixed-top" role="navigation">
  <div class="container">
    <div class="navbar-header">
      <a class="navbar-brand" href="../IndexPage/">
        <img src="../include/gene.png" height="28" width="28">
      </a>
    </div>
    <div id="navbar" class="navbar-collapse collapse">
      <ul class="nav navbar-nav">
        <li><a href="../IndexPage/">Home</a></li>', file = before_body)

  for (page in index) {    
    if (!is.null(page$children)) {
      # collection of pages
      cat(glue('
        <li class="dropdown">
          <a href="authoring" class="dropdown-toggle" data-toggle="dropdown" role="button" aria-expanded="false">
            {page$title}<span class="caret"></span>
          </a>
          <ul class="dropdown-menu" role="menu">'), file = before_body, append=TRUE)
      for (child in page$children) {
        cat(glue(
            "<li><a href={child$url}>{child$title}</a></li>"),
          file = before_body, append = TRUE)
      }
      cat('
          </ul>
        </li>', file=before_body, append=TRUE)
    } else {
      # single page
      cat(glue("<li><a href={page$url}>{page$title}</a></li>"),
        file = before_body, append = TRUE)
    }
  }

  cat('
      </ul>
    </div><!--/.nav-collapse -->
  </div><!--/.container -->
</div><!--/.navbar -->
', file = before_body, append = TRUE)
})
