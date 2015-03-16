(TeX-add-style-hook
 "mystyle"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("geometry" "left=2cm" "right=2cm" "top=2cm" "bottom=2cm") ("matlab-prettifier" "numbered" "framed")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (TeX-run-style-hooks
    "inputenc"
    "amsmath"
    "amsfonts"
    "amssymb"
    "graphicx"
    "geometry"
    "subfiles"
    "listings"
    "xcolor"
    "cancel"
    "bm"
    "matlab-prettifier"
    "array")))

