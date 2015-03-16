(TeX-add-style-hook
 "main"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "12pt" "letterpaper")))
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art12"
    "mystyle")))

