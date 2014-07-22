
Update the web page:

(1) install pkg2html from R-forge

(2) install jekyll from http://jekyllrb.com/

(3) run setup.R in www_src

(4) in the newly created dir html, run 

    html> jekyll serve --watch

    and point you web browser to localhost:4000

(5) if everything is OK, run publish.R in www_src
    and commit changes in www
